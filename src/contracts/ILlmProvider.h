// contracts/ILlmProvider.h - FROZEN LLM-provider seam (Phase D, WP-J).
//
// An ILlmProvider turns a request (system prompt + conversation + tool specs)
// into one assistant turn, streaming text deltas as they arrive and returning a
// structured message that may carry tool-use requests. The Anthropic provider
// (libcurl + SSE Messages API) and a deterministic MockProvider both implement
// it; the tool-calling loop (agent::Agent) drives whichever is configured.
//
// SAFETY SCOPE: the provider is a transport. The safety boundary (no synthesis /
// route / precursor / manufacturability content) is enforced by the system
// prompt the loop supplies and by the absence of any such tool - see
// IAgentTools.h. Docking "binding affinity" is a pharmacology signal only.
#pragma once

#include <functional>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "contracts/IAgentTools.h"  // ToolSpec

namespace stimlab {

enum class ChatRole { System, User, Assistant, Tool };

// An assistant request to call a tool (Anthropic tool_use block).
struct ToolCall {
    std::string    id;        // provider-assigned tool_use id
    std::string    name;
    nlohmann::json arguments = nlohmann::json::object();
};

// A tool result fed back to the model (Anthropic tool_result block).
struct ToolResultBlock {
    std::string callId;       // the ToolCall::id this answers
    std::string content;
    bool        isError = false;
};

// One message in the conversation. A `Tool` message bundles the results of the
// previous assistant turn's tool calls (serialized as a user/tool_result turn).
struct ChatMessage {
    ChatRole                     role = ChatRole::User;
    std::string                  text;         // user/assistant free text
    std::vector<ToolCall>        toolCalls;    // assistant: requested calls
    std::vector<ToolResultBlock> toolResults;  // role==Tool: results
};

enum class StopReason { EndTurn, ToolUse, MaxTokens, Refusal, Error };

struct LlmUsage {
    int inputTokens = 0;
    int outputTokens = 0;
};

// One assembled assistant turn.
struct LlmResponse {
    ChatMessage message;                       // text + any tool calls
    StopReason  stop = StopReason::EndTurn;
    std::string error;                          // non-empty iff stop==Error
    LlmUsage    usage;
};

struct LlmRequest {
    std::string              system;
    std::vector<ChatMessage> messages;
    std::vector<ToolSpec>    tools;
    std::string              model = "claude-opus-4-8";
    int                      maxTokens = 1024;
};

// Streamed text-delta sink. NOTE: may be invoked from a worker thread - the
// callback must be thread-safe with respect to whatever it touches.
using StreamCallback = std::function<void(const std::string& textDelta)>;

class ILlmProvider {
public:
    virtual ~ILlmProvider() = default;

    virtual std::string id() const = 0;            // "anthropic" / "mock"
    virtual std::string displayName() const = 0;   // human label
    virtual bool        ready() const = 0;         // configured (key present, transport built)

    // Run one turn. Implementations must NOT throw across this boundary: on any
    // failure return a response with stop==Error and a populated `error`.
    virtual LlmResponse send(const LlmRequest& req, const StreamCallback& onText) const = 0;
};

}  // namespace stimlab
