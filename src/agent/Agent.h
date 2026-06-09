// agent/Agent.h - the tool-calling loop / orchestrator (WP-J).
//
// Drives an ILlmProvider against an IToolRegistry: append the user turn, call
// the provider (streaming text), run any requested tools (autopilot), feed the
// results back, repeat until the model stops asking for tools. Supports an
// ask-first mode that pauses for user approval before running a tool batch.
//
// Threading: each submit() runs the loop on a short-lived worker thread so the
// UI never blocks on the network. All observable state is mutex-guarded and read
// via snapshot(). Tools that mutate UI state must marshal back to the UI thread
// themselves (the highlight/navigate tools post to a thread-safe inbox on
// AppShell) - the Agent never touches ImGui. runBlocking() is a synchronous
// entry for tests (no thread, auto-approves tools).
#pragma once

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include <nlohmann/json.hpp>

#include "contracts/IAgentTools.h"
#include "contracts/ILlmProvider.h"

namespace stimlab::agent {

enum class AgentMode { Autopilot, AskFirst };
enum class AgentStatus { Idle, Running, AwaitingApproval, Done, Error };

struct TranscriptEntry {
    enum class Kind { User, Assistant, Tool, System, Error };
    Kind        kind = Kind::Assistant;
    std::string text;
};

struct PendingCall {
    std::string    id;
    std::string    name;
    nlohmann::json arguments;
};

// An immutable copy of the renderable state, taken under the lock.
struct AgentSnapshot {
    AgentStatus                  status = AgentStatus::Idle;
    std::vector<TranscriptEntry> transcript;
    std::string                  streaming;   // partial assistant text mid-stream
    std::vector<PendingCall>     pending;     // AskFirst: awaiting approval
    std::string                  lastError;
    int                          turns = 0;   // assistant turns completed
};

class Agent {
public:
    Agent() = default;
    ~Agent();
    Agent(const Agent&) = delete;
    Agent& operator=(const Agent&) = delete;

    // Configuration. Provider + registry must outlive the Agent. Call when idle.
    void configure(const ILlmProvider* provider, const IToolRegistry* registry,
                   std::string systemPrompt);
    void setProvider(const ILlmProvider* provider);  // hot-swap when idle (key entered)
    void setModel(std::string model);
    void setMode(AgentMode mode) { mode_.store(mode); }
    [[nodiscard]] AgentMode mode() const { return mode_.load(); }

    // Drive (threaded). submit() is ignored while a turn is in flight.
    void submit(const std::string& userText);
    void approvePending();   // AskFirst: run the pending tool batch
    void denyPending();      // AskFirst: decline it (model is told)
    void reset();            // clear the conversation (when idle)
    void poll();             // call each frame: joins a finished worker

    [[nodiscard]] AgentSnapshot snapshot() const;
    [[nodiscard]] bool          busy() const;

    // Synchronous turn for tests: no thread, tools auto-approved.
    void runBlocking(const std::string& userText);

private:
    void processTurn(bool threaded);
    bool waitForApproval(const std::vector<ToolCall>& calls);  // threaded AskFirst
    [[nodiscard]] bool abortRequested() const;

    // Config (set when idle; read on the worker after a happens-before submit()).
    const ILlmProvider*     provider_ = nullptr;
    const IToolRegistry*    registry_ = nullptr;
    std::atomic<AgentMode>  mode_{AgentMode::Autopilot};

    mutable std::mutex          mu_;
    std::condition_variable     cv_;
    std::string                 system_;
    std::string                 model_ = "claude-opus-4-8";
    std::vector<ChatMessage>    history_;
    std::vector<TranscriptEntry> transcript_;
    std::string                 streaming_;
    std::vector<PendingCall>    pending_;
    AgentStatus                 status_ = AgentStatus::Idle;
    std::string                 lastError_;
    int                         turns_ = 0;
    int                         approval_ = 0;   // 0 none, 1 approve, -1 deny
    bool                        abort_ = false;

    std::thread        worker_;
    std::atomic<bool>  threadDone_{false};
};

}  // namespace stimlab::agent
