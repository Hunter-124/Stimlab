// ui/AppShell.h - top-level UI: host dockspace, left navigator, central content
// router, right assistant. Holds the shared UiState, the assistant->UI highlight
// mechanism, and (Phase D) the real LLM agent: an ILlmProvider + tool registry +
// tool-calling loop whose highlight_panel / navigate_ui tools drive this shell.
#pragma once

#include <atomic>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <thread>
#include <vector>

#include "contracts/Services.h"
#include "data/Domain.h"
#include "workflow/Dag.h"

struct ID3D11Device;
struct ID3D11DeviceContext;

namespace stimlab {

class Config;
class IToolRegistry;  // contracts/IAgentTools.h
namespace render { class MolViewport; }
namespace docking { class Provisioner; }
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
    std::string group;  // navigator section, e.g. "Compound" / "Analysis" / "Discovery"
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

    // Docking engine + receptor provisioning (WP-F/WP-G). The Provisioner runs the
    // slow best-effort work off the UI thread; the panel polls its status. The first
    // access kicks a cheap locate-only probe (no network) to seed initial status.
    [[nodiscard]] docking::Provisioner& provisioner();
    void provisionDocking();  // start a full provision (download engine + headline receptors)

    // On-demand provisioning of a single selected target (any of the 29 CNS presets,
    // not just the 4 headlines). `target` may be a preset id or display name; returns
    // false if it resolves to no preset. A no-op while a provision is already running.
    bool provisionTarget(const std::string& target);
    // Cache-only readiness probe (no network): is `target`'s receptor prepared on disk?
    [[nodiscard]] bool receptorReady(const std::string& target) const;

    // Provision the Vina-GPU (OpenCL) engine on demand: download its binaries + compile a
    // GPU-matched kernel off the UI thread. This is heavier than the Vina/receptor
    // provision (a kernel build), so it runs on its own worker. No-op while one is already
    // running. Poll vinaGpuProvisioning()/vinaGpuStatus(); vinaGpuReady() is a cheap
    // cache-only readiness check (no network).
    void                      provisionVinaGpu();
    [[nodiscard]] bool        vinaGpuProvisioning() const { return vinaGpuProvisioning_.load(); }
    [[nodiscard]] std::string vinaGpuStatus() const;
    [[nodiscard]] bool        vinaGpuReady() const;

    // Async docking: returns the latest completed result for (molecule, target),
    // computed ONCE on a background thread (a real engine dock must never run on the
    // UI thread). `computing` is set while a fresh result is being produced. Safe and
    // cheap to call every frame; results are cached by molecule+target.
    DockJobResult dockFor(const Molecule& m, const std::string& target, bool& computing);

    // Drop the cached dock so the next dockFor recomputes - call after provisioning
    // completes so a result computed while the engine was absent is re-docked for real.
    void invalidateDock();

    // ---- Workflows / DAG (WP-D) ----------------------------------------------
    // Run the prep->dock pipeline as a live DAG off the UI thread. The Workflows
    // panel polls workflowSnapshot() each frame for a thread-safe view of node states.
    struct WfNode {
        std::string              id;
        std::string              module;
        workflow::NodeStatus     status = workflow::NodeStatus::Pending;
        std::string              detail;
        std::vector<std::string> deps;
    };
    struct WorkflowSnapshot {
        bool                running = false;
        bool                hasResult = false;
        std::string         label;
        std::vector<WfNode> nodes;
        int                 ran = 0, cached = 0, failed = 0, cancelled = 0;
    };
    void runWorkflow(const std::string& smiles, const std::string& targetId,
                     const std::string& label);
    void cancelWorkflow();
    [[nodiscard]] WorkflowSnapshot workflowSnapshot() const;
    [[nodiscard]] bool workflowRunning() const;

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
    // The agent's tool registry (enumeration / direct dispatch, e.g. tests).
    [[nodiscard]] IToolRegistry* toolRegistry();

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
    // Docking compute mode (0=Auto, 1=GPU, 2=CPU): persisted via Config and pushed to
    // the docking module's engine selection. Changing it re-docks under the new engine.
    void setComputeMode(int mode);
    [[nodiscard]] int  computeMode() const;
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
    void drawCommandPalette();

    void buildAgent();             // construct registry/tools/providers/agent (ctor)
    void registerAgentServiceTools();  // bind dock/admet/analyze/run/list tools to Services
    void registerAgentWebTools();      // bind web_search/web_fetch (science builds only)
    // Resolve a tool's compound argument: "" -> active; else a library name/id; else a
    // raw SMILES analyzed on the fly. nullopt if it is none of these.
    [[nodiscard]] std::optional<Molecule> resolveAgentCompound(const std::string& arg) const;
    void reconfigureAgent();       // pick active provider + mode from Config
    void drainAgentActions();      // apply queued highlight/navigate on the UI thread
    [[nodiscard]] std::string buildSystemPrompt() const;

    Services svc_;
    UiState  state_;
    std::vector<PanelInfo> panels_;

    ID3D11Device*        renderDev_ = nullptr;
    ID3D11DeviceContext* renderCtx_ = nullptr;
    std::unique_ptr<render::MolViewport> viewer_;

    std::unique_ptr<docking::Provisioner> provisioner_;
    bool provisionProbed_ = false;  // locate-only probe kicked once on first access

    // Vina-GPU (OpenCL) provisioning runs on its OWN worker (download + kernel compile is
    // heavier than the Vina/receptor provision); joined in the destructor. The Settings
    // panel polls vinaGpuStatus_ each frame.
    std::thread        vinaGpuThread_;
    std::atomic<bool>  vinaGpuProvisioning_{false};
    mutable std::mutex vinaGpuMu_;
    std::string        vinaGpuStatus_{"Vina-GPU (OpenCL) engine not provisioned."};

    // Async docking cache (single-slot). The worker computes dockDetailed() off the
    // UI thread; the panel reads the cached copy. Joined in the destructor.
    std::thread           dockThread_;
    mutable std::mutex    dockMu_;
    std::string           dockKey_;        // molecule|target currently cached/computing
    DockJobResult         dockResult_;     // last completed result
    bool                  dockComputing_ = false;
    bool                  dockReady_ = false;

    // Workflow run state (WP-D). The DAG runs on wfThread_ (off the UI thread);
    // onProgress + the final result update wfNodes_/counters under wfMu_. Joined in
    // the destructor. wfJobs_/wfCache_ are created lazily and outlive the worker.
    std::unique_ptr<workflow::JobSystem>     wfJobs_;
    std::unique_ptr<workflow::DiskNodeCache> wfCache_;
    std::thread                              wfThread_;
    mutable std::mutex                       wfMu_;
    std::vector<WfNode>                      wfNodes_;
    std::string                              wfLabel_;
    bool                                     wfRunning_ = false;
    bool                                     wfHasResult_ = false;
    int                                      wfRan_ = 0, wfCached_ = 0, wfFailed_ = 0,
                                             wfCancelled_ = 0;
    workflow::CancelToken                    wfCancel_;

    // ---- agent ----
    // DECLARATION ORDER IS LOAD-BEARING: every member the worker thread touches
    // (registry, providers, and the UI inbox below) is declared BEFORE agent_, so
    // ~Agent - which joins the worker - runs while they are still alive. agent_ is
    // therefore the LAST member.
    Config* config_ = nullptr;
    bool    agentUsingAnthropic_ = false;
    char    chatBuf_[1024] = {0};

    // Command palette state.
    bool cmdPaletteOpen_ = false;
    char cmdPaletteBuf_[128] = {};

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
