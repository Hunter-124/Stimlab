#include "agent/Agent.h"

#include <utility>

namespace biocad::agent {
namespace {

constexpr int kMaxToolIterations = 8;

ChatMessage makeUser(const std::string& text) {
    ChatMessage m;
    m.role = ChatRole::User;
    m.text = text;
    return m;
}

}  // namespace

Agent::~Agent() {
    {
        std::lock_guard<std::mutex> lk(mu_);
        abort_ = true;
        approval_ = -1;  // unblock an AskFirst wait
    }
    cv_.notify_all();
    if (worker_.joinable()) worker_.join();
}

void Agent::configure(const ILlmProvider* provider, const IToolRegistry* registry,
                      std::string systemPrompt) {
    std::lock_guard<std::mutex> lk(mu_);
    provider_ = provider;
    registry_ = registry;
    system_ = std::move(systemPrompt);
}

void Agent::setProvider(const ILlmProvider* provider) {
    if (busy()) return;
    poll();
    std::lock_guard<std::mutex> lk(mu_);
    provider_ = provider;
}

void Agent::setModel(std::string model) {
    std::lock_guard<std::mutex> lk(mu_);
    model_ = std::move(model);
}

bool Agent::busy() const {
    std::lock_guard<std::mutex> lk(mu_);
    return status_ == AgentStatus::Running || status_ == AgentStatus::AwaitingApproval;
}

bool Agent::abortRequested() const {
    std::lock_guard<std::mutex> lk(mu_);
    return abort_;
}

void Agent::poll() {
    if (threadDone_.load() && worker_.joinable()) {
        worker_.join();
        threadDone_.store(false);
    }
}

void Agent::submit(const std::string& userText) {
    if (userText.empty() || busy()) return;
    poll();  // join a previously finished worker before reusing worker_
    {
        std::lock_guard<std::mutex> lk(mu_);
        if (!provider_) return;
        history_.push_back(makeUser(userText));
        transcript_.push_back({TranscriptEntry::Kind::User, userText});
        streaming_.clear();
        status_ = AgentStatus::Running;
        abort_ = false;
    }
    threadDone_.store(false);
    worker_ = std::thread([this] {
        processTurn(/*threaded=*/true);
        threadDone_.store(true);
    });
}

void Agent::runBlocking(const std::string& userText) {
    {
        std::lock_guard<std::mutex> lk(mu_);
        history_.push_back(makeUser(userText));
        transcript_.push_back({TranscriptEntry::Kind::User, userText});
        streaming_.clear();
        status_ = AgentStatus::Running;
        abort_ = false;
    }
    processTurn(/*threaded=*/false);
}

void Agent::approvePending() {
    {
        std::lock_guard<std::mutex> lk(mu_);
        approval_ = 1;
    }
    cv_.notify_all();
}

void Agent::denyPending() {
    {
        std::lock_guard<std::mutex> lk(mu_);
        approval_ = -1;
    }
    cv_.notify_all();
}

void Agent::reset() {
    if (busy()) return;
    poll();
    std::lock_guard<std::mutex> lk(mu_);
    history_.clear();
    transcript_.clear();
    pending_.clear();
    streaming_.clear();
    lastError_.clear();
    status_ = AgentStatus::Idle;
    turns_ = 0;
}

bool Agent::waitForApproval(const std::vector<ToolCall>& calls) {
    std::unique_lock<std::mutex> lk(mu_);
    pending_.clear();
    for (const auto& c : calls) pending_.push_back({c.id, c.name, c.arguments});
    status_ = AgentStatus::AwaitingApproval;
    approval_ = 0;
    cv_.wait(lk, [this] { return approval_ != 0 || abort_; });
    if (abort_) return false;
    const bool ok = (approval_ == 1);
    approval_ = 0;
    status_ = AgentStatus::Running;
    return ok;
}

AgentSnapshot Agent::snapshot() const {
    std::lock_guard<std::mutex> lk(mu_);
    AgentSnapshot s;
    s.status = status_;
    s.transcript = transcript_;
    s.streaming = streaming_;
    s.pending = pending_;
    s.lastError = lastError_;
    s.turns = turns_;
    return s;
}

void Agent::processTurn(bool threaded) {
    const ILlmProvider* provider = provider_;  // set before submit(); stable for this turn
    const IToolRegistry* registry = registry_;
    const bool askFirst = (mode_.load() == AgentMode::AskFirst);

    if (!provider) {
        std::lock_guard<std::mutex> lk(mu_);
        status_ = AgentStatus::Error;
        lastError_ = "no provider configured";
        transcript_.push_back({TranscriptEntry::Kind::Error, "No assistant provider configured."});
        return;
    }

    for (int iter = 0; iter < kMaxToolIterations; ++iter) {
        LlmRequest req;
        {
            std::lock_guard<std::mutex> lk(mu_);
            req.system = system_;
            req.model = model_;
            req.messages = history_;
        }
        req.maxTokens = 1024;
        if (registry) req.tools = registry->specs();

        // Stream text deltas into streaming_ (read each frame by the UI).
        LlmResponse resp = provider->send(req, [this](const std::string& delta) {
            std::lock_guard<std::mutex> lk(mu_);
            streaming_ += delta;
        });

        if (resp.stop == StopReason::Error) {
            std::lock_guard<std::mutex> lk(mu_);
            status_ = AgentStatus::Error;
            lastError_ = resp.error;
            transcript_.push_back({TranscriptEntry::Kind::Error, "Error: " + resp.error});
            streaming_.clear();
            return;
        }

        // Commit the assistant message.
        {
            std::lock_guard<std::mutex> lk(mu_);
            history_.push_back(resp.message);
            if (!resp.message.text.empty())
                transcript_.push_back({TranscriptEntry::Kind::Assistant, resp.message.text});
            streaming_.clear();
            ++turns_;
        }

        if (resp.stop != StopReason::ToolUse || resp.message.toolCalls.empty()) {
            std::lock_guard<std::mutex> lk(mu_);
            status_ = AgentStatus::Done;
            return;
        }

        // The model asked to run one or more tools.
        const std::vector<ToolCall>& calls = resp.message.toolCalls;
        bool approved = true;
        if (askFirst && threaded) {
            approved = waitForApproval(calls);
            if (abortRequested()) return;
        }

        ChatMessage toolMsg;
        toolMsg.role = ChatRole::Tool;
        for (const auto& c : calls) {
            ToolResult out;
            if (!approved) {
                out.content = "User declined to run this tool.";
            } else if (registry && registry->has(c.name)) {
                out = registry->dispatch(c.name, c.arguments);
            } else {
                out.content = "Unknown tool: " + c.name;
                out.isError = true;
            }
            toolMsg.toolResults.push_back({c.id, out.content, out.isError});
            std::lock_guard<std::mutex> lk(mu_);
            transcript_.push_back(
                {TranscriptEntry::Kind::Tool, (approved ? "\xE2\x86\x92 " : "declined ") + c.name});
        }
        {
            std::lock_guard<std::mutex> lk(mu_);
            history_.push_back(std::move(toolMsg));
            pending_.clear();
            status_ = AgentStatus::Running;
        }
        // loop: let the model react to the tool results
    }

    std::lock_guard<std::mutex> lk(mu_);
    status_ = AgentStatus::Done;
    transcript_.push_back(
        {TranscriptEntry::Kind::System, "(stopped: reached the tool-call limit for one turn)"});
}

}  // namespace biocad::agent
