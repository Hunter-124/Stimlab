#include "agent/Tools.h"

#include <utility>

namespace stimlab::agent {

FunctionTool::FunctionTool(std::string name, std::string description, nlohmann::json schema,
                           Handler handler)
    : handler_(std::move(handler)) {
    spec_.name = std::move(name);
    spec_.description = std::move(description);
    spec_.inputSchema = std::move(schema);
}

ToolResult FunctionTool::invoke(const nlohmann::json& args) const {
    if (!handler_) return {"tool has no handler", true};
    try {
        return handler_(args);
    } catch (const std::exception& e) {
        return {std::string("tool error: ") + e.what(), true};
    }
}

void ToolRegistry::add(std::unique_ptr<ITool> tool) {
    if (tool) tools_.push_back(std::move(tool));
}

std::vector<ToolSpec> ToolRegistry::specs() const {
    std::vector<ToolSpec> out;
    out.reserve(tools_.size());
    for (const auto& t : tools_) out.push_back(t->spec());
    return out;
}

bool ToolRegistry::has(const std::string& name) const {
    for (const auto& t : tools_)
        if (t->spec().name == name) return true;
    return false;
}

ToolResult ToolRegistry::dispatch(const std::string& name, const nlohmann::json& args) const {
    for (const auto& t : tools_)
        if (t->spec().name == name) return t->invoke(args);
    return {"unknown tool: " + name, true};
}

}  // namespace stimlab::agent
