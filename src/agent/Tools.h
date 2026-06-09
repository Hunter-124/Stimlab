// agent/Tools.h - concrete IToolRegistry + a std::function-backed ITool.
//
// FunctionTool lets the UI layer bind tools to AppShell (highlight_panel, etc.)
// via a closure, so the agent library never depends on the UI. ToolRegistry owns
// the tools and dispatches by name.
#pragma once

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "contracts/IAgentTools.h"

namespace stimlab::agent {

class FunctionTool : public ITool {
public:
    using Handler = std::function<ToolResult(const nlohmann::json&)>;

    FunctionTool(std::string name, std::string description, nlohmann::json schema,
                 Handler handler);

    ToolSpec   spec() const override { return spec_; }
    ToolResult invoke(const nlohmann::json& args) const override;

private:
    ToolSpec spec_;
    Handler  handler_;
};

class ToolRegistry : public IToolRegistry {
public:
    void add(std::unique_ptr<ITool> tool);

    std::vector<ToolSpec> specs() const override;
    bool                  has(const std::string& name) const override;
    ToolResult            dispatch(const std::string& name, const nlohmann::json& args) const override;

    [[nodiscard]] std::size_t size() const { return tools_.size(); }

private:
    std::vector<std::unique_ptr<ITool>> tools_;
};

}  // namespace stimlab::agent
