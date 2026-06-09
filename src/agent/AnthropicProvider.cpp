#include "agent/AnthropicProvider.h"

#include <utility>

#ifdef STIMLAB_HAVE_SCIENCE
#  include <mutex>
#  include <vector>
#  include <curl/curl.h>
#  include <nlohmann/json.hpp>
#endif

namespace stimlab::agent {

void AnthropicProvider::setApiKey(std::string key) {
    std::lock_guard<std::mutex> lk(mu_);
    apiKey_ = std::move(key);
}

bool AnthropicProvider::hasKey() const {
    std::lock_guard<std::mutex> lk(mu_);
    return !apiKey_.empty();
}

bool AnthropicProvider::transportAvailable() {
#ifdef STIMLAB_HAVE_SCIENCE
    return true;
#else
    return false;
#endif
}

bool AnthropicProvider::ready() const {
#ifdef STIMLAB_HAVE_SCIENCE
    return hasKey();
#else
    return false;
#endif
}

#ifdef STIMLAB_HAVE_SCIENCE

namespace {

using nlohmann::json;

StopReason mapStop(const std::string& s) {
    if (s == "tool_use") return StopReason::ToolUse;
    if (s == "max_tokens") return StopReason::MaxTokens;
    if (s == "refusal") return StopReason::Refusal;
    return StopReason::EndTurn;  // end_turn / stop_sequence / unknown
}

// Serialize the conversation to the Anthropic messages[] shape.
json buildMessages(const LlmRequest& req) {
    json msgs = json::array();
    for (const auto& m : req.messages) {
        switch (m.role) {
            case ChatRole::User: {
                msgs.push_back({{"role", "user"},
                                {"content", json::array({{{"type", "text"}, {"text", m.text}}})}});
                break;
            }
            case ChatRole::Assistant: {
                json content = json::array();
                if (!m.text.empty()) content.push_back({{"type", "text"}, {"text", m.text}});
                for (const auto& c : m.toolCalls)
                    content.push_back({{"type", "tool_use"},
                                       {"id", c.id},
                                       {"name", c.name},
                                       {"input", c.arguments}});
                if (content.empty()) content.push_back({{"type", "text"}, {"text", ""}});
                msgs.push_back({{"role", "assistant"}, {"content", content}});
                break;
            }
            case ChatRole::Tool: {
                json content = json::array();
                for (const auto& tr : m.toolResults)
                    content.push_back({{"type", "tool_result"},
                                       {"tool_use_id", tr.callId},
                                       {"content", tr.content},
                                       {"is_error", tr.isError}});
                msgs.push_back({{"role", "user"}, {"content", content}});
                break;
            }
            case ChatRole::System:
                break;  // system is top-level, not a message
        }
    }
    return msgs;
}

std::string buildBody(const LlmRequest& req) {
    json b;
    b["model"] = req.model;
    b["max_tokens"] = req.maxTokens;
    b["stream"] = true;
    // NOTE: deliberately NO temperature/top_p/top_k - Opus 4.x rejects sampling
    // params (400). Thinking is left off (this is a snappy UI navigator).
    if (!req.system.empty()) b["system"] = req.system;
    if (!req.tools.empty()) {
        json tools = json::array();
        for (const auto& t : req.tools)
            tools.push_back(
                {{"name", t.name}, {"description", t.description}, {"input_schema", t.inputSchema}});
        b["tools"] = tools;
    }
    b["messages"] = buildMessages(req);
    return b.dump();
}

// Accumulates the streamed response as SSE lines arrive.
struct SseCtx {
    std::string lineBuf;       // partial trailing line between write callbacks
    std::string raw;           // full raw body (used for the non-200 error path)

    std::string text;          // assembled assistant text
    std::vector<ToolCall> calls;

    // current content block
    std::string curType;       // "text" / "tool_use" / ""
    std::string curToolId, curToolName, curPartial;

    StopReason stop = StopReason::EndTurn;
    int inputTokens = 0, outputTokens = 0;
    std::string errMsg;

    const StreamCallback* onText = nullptr;
};

void handleEvent(SseCtx& ctx, const std::string& payload) {
    json j = json::parse(payload, nullptr, /*allow_exceptions=*/false);
    if (j.is_discarded() || !j.is_object()) return;
    const std::string type = j.value("type", "");

    if (type == "message_start") {
        if (j.contains("message") && j["message"].contains("usage"))
            ctx.inputTokens = j["message"]["usage"].value("input_tokens", 0);
    } else if (type == "content_block_start") {
        const json& cb = j["content_block"];
        ctx.curType = cb.value("type", "");
        if (ctx.curType == "tool_use") {
            ctx.curToolId = cb.value("id", "");
            ctx.curToolName = cb.value("name", "");
            ctx.curPartial.clear();
        }
    } else if (type == "content_block_delta") {
        const json& d = j["delta"];
        const std::string dt = d.value("type", "");
        if (dt == "text_delta") {
            const std::string t = d.value("text", "");
            ctx.text += t;
            if (ctx.onText && *ctx.onText) (*ctx.onText)(t);
        } else if (dt == "input_json_delta") {
            ctx.curPartial += d.value("partial_json", "");
        }
    } else if (type == "content_block_stop") {
        if (ctx.curType == "tool_use") {
            ToolCall c;
            c.id = ctx.curToolId;
            c.name = ctx.curToolName;
            json args = ctx.curPartial.empty()
                            ? json::object()
                            : json::parse(ctx.curPartial, nullptr, /*allow_exceptions=*/false);
            c.arguments = args.is_discarded() ? json::object() : args;
            ctx.calls.push_back(std::move(c));
        }
        ctx.curType.clear();
    } else if (type == "message_delta") {
        if (j.contains("delta") && j["delta"].contains("stop_reason") &&
            j["delta"]["stop_reason"].is_string())
            ctx.stop = mapStop(j["delta"]["stop_reason"].get<std::string>());
        if (j.contains("usage")) ctx.outputTokens = j["usage"].value("output_tokens", ctx.outputTokens);
    } else if (type == "error") {
        ctx.errMsg = j.contains("error") ? j["error"].value("message", "stream error") : "stream error";
        ctx.stop = StopReason::Error;
    }
}

void handleLine(SseCtx& ctx, std::string line) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.rfind("data:", 0) != 0) return;  // ignore "event:" lines + blanks
    std::string payload = line.substr(5);
    if (!payload.empty() && payload.front() == ' ') payload.erase(0, 1);
    if (payload.empty() || payload == "[DONE]") return;
    handleEvent(ctx, payload);
}

std::size_t writeCb(char* ptr, std::size_t size, std::size_t nmemb, void* userdata) {
    const std::size_t n = size * nmemb;
    auto* ctx = static_cast<SseCtx*>(userdata);
    ctx->raw.append(ptr, n);
    ctx->lineBuf.append(ptr, n);
    std::size_t pos;
    while ((pos = ctx->lineBuf.find('\n')) != std::string::npos) {
        handleLine(*ctx, ctx->lineBuf.substr(0, pos));
        ctx->lineBuf.erase(0, pos + 1);
    }
    return n;
}

void ensureCurlInit() {
    static std::once_flag once;
    std::call_once(once, [] { curl_global_init(CURL_GLOBAL_DEFAULT); });
}

}  // namespace

LlmResponse AnthropicProvider::send(const LlmRequest& req, const StreamCallback& onText) const {
    LlmResponse r;
    std::string key;
    {
        std::lock_guard<std::mutex> lk(mu_);
        key = apiKey_;
    }
    if (key.empty()) {
        r.stop = StopReason::Error;
        r.error = "no Anthropic API key configured (enter one in Settings).";
        return r;
    }

    ensureCurlInit();
    CURL* curl = curl_easy_init();
    if (!curl) {
        r.stop = StopReason::Error;
        r.error = "could not initialize libcurl.";
        return r;
    }

    const std::string body = buildBody(req);

    SseCtx ctx;
    ctx.onText = &onText;

    curl_slist* headers = nullptr;
    headers = curl_slist_append(headers, "content-type: application/json");
    headers = curl_slist_append(headers, "anthropic-version: 2023-06-01");
    const std::string keyHeader = "x-api-key: " + key;
    headers = curl_slist_append(headers, keyHeader.c_str());

    curl_easy_setopt(curl, CURLOPT_URL, "https://api.anthropic.com/v1/messages");
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
    curl_easy_setopt(curl, CURLOPT_POST, 1L);
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, body.c_str());
    curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, static_cast<long>(body.size()));
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writeCb);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &ctx);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, 120L);
    curl_easy_setopt(curl, CURLOPT_CONNECTTIMEOUT, 20L);
    curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 1L);
    curl_easy_setopt(curl, CURLOPT_USERAGENT, "StimLab/0.1");
    curl_easy_setopt(curl, CURLOPT_ACCEPT_ENCODING, "");  // allow gzip if offered

    const CURLcode rc = curl_easy_perform(curl);
    long httpCode = 0;
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &httpCode);
    curl_slist_free_all(headers);
    curl_easy_cleanup(curl);

    if (rc != CURLE_OK) {
        r.stop = StopReason::Error;
        r.error = std::string("network error: ") + curl_easy_strerror(rc);
        return r;
    }
    if (httpCode != 200) {
        std::string msg = "HTTP " + std::to_string(httpCode);
        json j = json::parse(ctx.raw, nullptr, false);
        if (!j.is_discarded() && j.contains("error") && j["error"].contains("message"))
            msg += ": " + j["error"]["message"].get<std::string>();
        r.stop = StopReason::Error;
        r.error = msg;
        return r;
    }
    if (!ctx.errMsg.empty()) {
        r.stop = StopReason::Error;
        r.error = ctx.errMsg;
        return r;
    }

    r.message.role = ChatRole::Assistant;
    r.message.text = ctx.text;
    r.message.toolCalls = std::move(ctx.calls);
    r.stop = ctx.stop;
    if (!r.message.toolCalls.empty() && r.stop == StopReason::EndTurn) r.stop = StopReason::ToolUse;
    r.usage.inputTokens = ctx.inputTokens;
    r.usage.outputTokens = ctx.outputTokens;
    return r;
}

#else  // !STIMLAB_HAVE_SCIENCE - no transport compiled in.

LlmResponse AnthropicProvider::send(const LlmRequest&, const StreamCallback&) const {
    LlmResponse r;
    r.stop = StopReason::Error;
    r.error =
        "Anthropic provider unavailable: this build has no networking. Rebuild with the "
        "windows-science preset (libcurl) to enable the live assistant.";
    return r;
}

#endif

}  // namespace stimlab::agent
