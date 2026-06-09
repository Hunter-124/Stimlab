// ui/AppShell.h - top-level UI: host dockspace, left navigator, central content
// router, right assistant. Holds the shared UiState, the assistant->UI highlight
// mechanism, and (Phase D) the real LLM agent: an ILlmProvider + tool registry +
// tool-calling loop whose highlight_panel / navigate_ui tools drive this shell.
#pragma once

#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "contracts/Services.h"
#include "data/Domain.h"

struct ID3D11Device;
struct ID3D11DeviceContext;

namespace stimlab {

class Config;
namespace render { class MolViewport; }
namespace agent {
class Agent;
class ToolRegistry;
class MockProvider;
class AnthropicProvider;
}  // namespace agent


struct PanelInfo {
    std::string id;     // stable key used by the content router + highlight
    std::string label;  // display label
    std::string help;   // one-line description shown by the assistant
};

struct UiState {
    std::string selectedMolecule = "amphetamine";
    Molecule    customMolecule;             // user-entered SMILES result (id "__custom__")
    bool        hasCustom = false;          // customMolecule analyzed + selectable as active
    std::string activePanel = "Dashboard";
    std::string dockTarget;                 // chosen docking target
    std::string highlight;                  // panel id currently pulsing ("" = none)
    double      highlightStart = -1000.0;   // ImGui time when highlight began
    bool        showAssistant = true;
    bool        showAbout = false;
    bool        quitRequested = false;
    std::vector<std::string> assistantLog;  // legacy canned-button log (kept as fallback)
};

class AppShell {
public:
    explicit AppShell(Services services);
    ~AppShell();  // defined in the .cpp (MolViewport + agent types are incomplete here)

    void draw();  // call once per frame between NewFrame() and Render()

    UiState&  state()    { return state_; }
    Services& services() { return svc_; }
    [[nodiscard]] Molecule currentMolecule() const;

    // 3D viewer wiring (WP-2). WinMain hands us the live DX11 device so panels can
    // render the off-screen molecular viewport via viewer().
    void setRenderDevice(ID3D11Device* device, ID3D11DeviceContext* context);
    [[nodiscard]] render::MolViewport* viewer();

    // Assistant -> UI bridge.
    void requestHighlight(const std::string& panelId, const std::string& explanation);
    void frameHighlightCurrentWindow(const std::string& panelId);  // call inside a panel
    [[nodiscard]] bool isHighlighted(const std::string& panelId) const;

    [[nodiscard]] const std::vector<PanelInfo>& panels() const { return panels_; }

    // ---- Phase D agent --------------------------------------------------------
    // WinMain hands us the persisted Config so the agent can read the API key /
    // provider / model and Settings can write them back.
    void setConfig(Config* config);
    [[nodiscard]] Config* config() { return config_; }
    [[nodiscard]] agent::Agent* agent() { return agent_.get(); }

    // Settings helpers (encapsulate Config + DPAPI + agent reconfigure).
    [[nodiscard]] bool hasApiKey() const;
    void saveApiKey(const std::string& plaintext);   // encrypts via DPAPI, persists, reconfigures
    void clearApiKey();
    void setAgentModel(const std::string& model);
    [[nodiscard]] std::string agentModel() const;
    void setAgentProviderIndex(int idx);             // 0 = Anthropic, 1 = Offline mock
    [[nodiscard]] int  agentProviderIndex() const;
    void setAutopilot(bool on);
    [[nodiscard]] bool autopilot() const;
    [[nodiscard]] bool anthropicReady() const;       // key present AND transport built
    [[nodiscard]] bool anthropicTransport() const;   // built with libcurl (science feature)
    [[nodiscard]] std::string activeProviderLabel() const;

    // Tool bridge (called from the agent worker thread - thread-safe, enqueues a
    // UI action applied on the UI thread at the top of draw()).
    void postHighlightAction(const std::string& panelId, const std::string& explanation);
    void postNavigateAction(const std::string& panelId);
    [[nodiscard]] bool isValidPanel(const std::string& panelId) const;
    [[nodiscard]] Molecule agentCompoundSnapshot() const;  // thread-safe copy for tools

private:
    void drawMainMenuBar();
    void drawNavigator();
    void drawContent();
    void drawAssistant();
    void drawAboutModal();

    void buildAgent();             // construct registry/tools/providers/agent (ctor)
    void reconfigureAgent();       // pick active provider + mode from Config
    void drainAgentActions();      // apply queued highlight/navigate on the UI thread
    [[nodiscard]] std::string buildSystemPrompt() const;

    Services svc_;
    UiState  state_;
    std::vector<PanelInfo> panels_;

    ID3D11Device*        renderDev_ = nullptr;
    ID3D11DeviceContext* renderCtx_ = nullptr;
    std::unique_ptr<render::MolViewport> viewer_;

    // ---- agent ----
    // DECLARATION ORDER IS LOAD-BEARING: every member the worker thread touches
    // (registry, providers, and the UI inbox below) is declared BEFORE agent_, so
    // ~Agent - which joins the worker - runs while they are still alive. agent_ is
    // therefore the LAST member.
    Config* config_ = nullptr;
    bool    agentUsingAnthropic_ = false;
    char    chatBuf_[1024] = {0};

    std::unique_ptr<agent::ToolRegistry>      registry_;
    std::unique_ptr<agent::MockProvider>      mock_;
    std::unique_ptr<agent::AnthropicProvider> anthropic_;

    // Thread-safe inbox between the agent worker and the UI thread.
    struct AgentUiAction {
        bool        navigate = false;  // true = navigate_ui, false = highlight_panel
        std::string panel;
        std::string explanation;
    };
    mutable std::mutex         agentMu_;
    std::vector<AgentUiAction> agentInbox_;
    Molecule                   agentCompound_;  // refreshed each frame on the UI thread

    std::unique_ptr<agent::Agent> agent_;  // LAST: ~Agent joins its worker first
};

}  // namespace stimlab
