// contracts/IAgentTools.h - FROZEN agent tool seam (Phase D).
//
// An ITool is a single capability the LLM assistant may invoke (name + JSON
// schema + an execute fn); an IToolRegistry is the set the agent loop dispatches
// against. Tools describe themselves to the provider via ToolSpec. Concrete
// tools (highlight_panel, navigate_ui, get_active_compound, ...) are built in the
// UI layer so they can bind to AppShell without the agent library depending on UI.
//
// SAFETY SCOPE: the tool surface is navigation + explanation + read-only
// analysis of what a compound IS and DOES. There is intentionally NO tool that
// produces synthesis routes, reaction conditions, precursors, or
// manufacturability guidance - that is out of scope by design and the agent's
// system prompt enforces it too.
#pragma once

#include <string>
#include <vector>

#include <nlohmann/json.hpp>

namespace stimlab {

// How a tool advertises itself to the LLM provider (Anthropic `tools[]` entry).
struct ToolSpec {
    std::string    name;          // unique, snake_case (e.g. "highlight_panel")
    std::string    description;   // tells the model when/why to call it
    nlohmann::json inputSchema = nlohmann::json::object();  // JSON Schema (object)
};

// The outcome of running a tool: the text handed back to the model as the
// tool_result content, plus whether it should be flagged is_error.
struct ToolResult {
    std::string content;
    bool        isError = false;
};

// One invokable capability.
class ITool {
public:
    virtual ~ITool() = default;
    virtual ToolSpec   spec() const = 0;
    virtual ToolResult invoke(const nlohmann::json& args) const = 0;
};

// The set of tools the agent loop can dispatch against.
class IToolRegistry {
public:
    virtual ~IToolRegistry() = default;
    virtual std::vector<ToolSpec> specs() const = 0;
    virtual bool        has(const std::string& name) const = 0;
    virtual ToolResult  dispatch(const std::string& name, const nlohmann::json& args) const = 0;
};

}  // namespace stimlab
