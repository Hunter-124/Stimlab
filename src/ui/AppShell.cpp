#include "ui/AppShell.h"

#include <cmath>
#include <utility>

#include <imgui.h>
#include <nlohmann/json.hpp>

#include "agent/Agent.h"
#include "agent/AnthropicProvider.h"
#include "agent/MockProvider.h"
#include "agent/SystemPrompt.h"
#include "agent/Tools.h"
#include "core/Config.h"
#include "core/Secrets.h"
#include "modules/Pipelines.h"
#include "modules/docking/Presets.h"
#include "modules/docking/Provisioning.h"
#include "render/MolViewport.h"
#include "ui/Panels.h"
#include "ui/Theme.h"
#include "workflow/Dag.h"
#include "workflow/JobSystem.h"

namespace stimlab {

namespace {
// Config keys for the persisted agent settings.
constexpr const char* kCfgProvider = "agent.provider";  // 0 = Anthropic, 1 = Offline mock
constexpr const char* kCfgApiKey   = "agent.apiKey";    // DPAPI-encrypted base64 blob
constexpr const char* kCfgModel    = "agent.model";
constexpr const char* kCfgMode     = "agent.mode";      // "autopilot" | "askfirst"
constexpr const char* kDefaultModel = "claude-opus-4-8";
}  // namespace

namespace {

constexpr double kHighlightSeconds = 5.0;

constexpr ImGuiWindowFlags kPaneFlags =
    ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse |
    ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoDocking;

void pulseBorder(const ImVec2& min, const ImVec2& max, double elapsed) {
    const float a = 0.45f + 0.45f * static_cast<float>(std::sin(elapsed * 6.0));
    const ImU32 col = IM_COL32(255, 191, 64, static_cast<int>(a * 255.0f));
    ImDrawList* dl = ImGui::GetForegroundDrawList();
    dl->AddRect(ImVec2(min.x - 2, min.y - 2), ImVec2(max.x + 2, max.y + 2), col,
                ImGui::GetStyle().WindowRounding, 0, 3.0f);
}

}  // namespace

AppShell::AppShell(Services services) : svc_(services) {
    panels_ = {
        {"Dashboard", "Dashboard",
         "Overview: library size, recent runs, and a snapshot of the selected compound."},
        {"Structure", "Structure Workbench",
         "Identity + physicochemical properties of the selected molecule with a live 3D viewer."},
        {"Input", "Molecule Input",
         "Enter any SMILES (or load a library compound) to analyze an arbitrary structure."},
        {"Analog", "Analog Explorer",
         "Model a candidate derivative's profile and check it vs existing samples + predicted byproducts."},
        {"Compare", "Compare",
         "Side-by-side comparison of up to three compounds across stability, absorption and ADMET."},
        {"Stability", "Stability",
         "Degradation liabilities (hydrolysis/oxidation/photolysis/thermal/pH) and a shelf-life estimate."},
        {"Absorption", "Absorption / PK",
         "Permeability, oral bioavailability, BBB penetration and P-gp efflux - to narrow candidates."},
        {"Metabolism", "Metabolism (ADMET)",
         "Metabolic routes, harmful metabolites, drug-drug interactions, hERG and safety flags."},
        {"Similarity", "Similarity",
         "Structural + pharmacophore similarity vs the curated known-substance reference set."},
        {"Legal", "Legal Analog",
         "Substantial-similarity scorecard vs controlled references (illustrative, not legal advice)."},
        {"Docking", "Docking",
         "Predicted ligand->target BINDING AFFINITY (pharmacology) at DAT/NET/SERT/TAAR1."},
        {"Workflows", "Workflows",
         "Re-runnable prep->dock pipeline as a live, content-cached, cancellable DAG."},
        {"Library", "Library",
         "Browse, search and select compounds from the default + imported library."},
        {"Runs", "Runs",
         "History of analyses with status and summaries."},
        {"Presets", "Presets / Targets",
         "CNS target presets and reusable analysis panels."},
        {"Settings", "Settings",
         "AI provider/keys, GPU mode, storage paths."},
    };
    provisioner_ = std::make_unique<docking::Provisioner>();
    buildAgent();
}

AppShell::~AppShell() {
    // Join workers BEFORE members tear down (they read svc_ / wfJobs_ / wfCache_).
    wfCancel_.cancel();
    if (wfThread_.joinable()) wfThread_.join();
    if (dockThread_.joinable()) dockThread_.join();
}

void AppShell::runWorkflow(const std::string& smiles, const std::string& targetId,
                           const std::string& label) {
    {
        std::lock_guard<std::mutex> lk(wfMu_);
        if (wfRunning_) return;  // one run at a time
    }
    if (wfThread_.joinable()) wfThread_.join();  // reap a previously-finished worker

    // Build the DAG on the UI thread (cheap: SMILES + filesystem probes, no network /
    // no engine subprocess - those run inside the node functions on the worker).
    workflow::Dag dag = workflow::buildDockingPipeline(smiles, targetId, svc_);
    {
        std::lock_guard<std::mutex> lk(wfMu_);
        wfNodes_.clear();
        for (const auto& n : dag.nodes())
            wfNodes_.push_back({n.id, n.module, workflow::NodeStatus::Pending, "", n.deps});
        wfLabel_ = label;
        wfRunning_ = true;
        wfHasResult_ = false;
        wfRan_ = wfCached_ = wfFailed_ = wfCancelled_ = 0;
    }
    wfCancel_ = workflow::CancelToken{};  // fresh token for this run
    if (!wfJobs_) wfJobs_ = std::make_unique<workflow::JobSystem>();
    if (!wfCache_) wfCache_ = std::make_unique<workflow::DiskNodeCache>("workflow");
    const workflow::CancelToken tok = wfCancel_;

    wfThread_ = std::thread([this, dag = std::move(dag), tok]() mutable {
        workflow::DagExecutor exec(*wfJobs_, *wfCache_);
        auto onProgress = [this](const workflow::NodeProgress& p) {
            std::lock_guard<std::mutex> lk(wfMu_);
            for (auto& n : wfNodes_)
                if (n.id == p.id) {
                    n.status = p.status;
                    n.detail = p.detail;
                    break;
                }
        };
        const workflow::DagRunResult r = exec.run(dag, tok, onProgress);
        std::lock_guard<std::mutex> lk(wfMu_);
        for (auto& n : wfNodes_) {
            n.status = r.statusOf(n.id);
            auto it = r.output.find(n.id);  // show the node's actual output (JSON)
            if (it != r.output.end() && !it->second.empty()) n.detail = it->second;
        }
        wfRan_ = r.ran;
        wfCached_ = r.cached;
        wfFailed_ = r.failed;
        wfCancelled_ = r.cancelled;
        wfRunning_ = false;
        wfHasResult_ = true;
    });
}

void AppShell::cancelWorkflow() { wfCancel_.cancel(); }

bool AppShell::workflowRunning() const {
    std::lock_guard<std::mutex> lk(wfMu_);
    return wfRunning_;
}

AppShell::WorkflowSnapshot AppShell::workflowSnapshot() const {
    std::lock_guard<std::mutex> lk(wfMu_);
    WorkflowSnapshot s;
    s.running = wfRunning_;
    s.hasResult = wfHasResult_;
    s.label = wfLabel_;
    s.nodes = wfNodes_;
    s.ran = wfRan_;
    s.cached = wfCached_;
    s.failed = wfFailed_;
    s.cancelled = wfCancelled_;
    return s;
}

DockJobResult AppShell::dockFor(const Molecule& m, const std::string& target, bool& computing) {
    const std::string key = m.id + "|" + m.smiles + "|" + target;
    std::lock_guard<std::mutex> lk(dockMu_);

    if (key == dockKey_) {
        computing = dockComputing_;
        return dockReady_ ? dockResult_ : DockJobResult{};
    }

    // New (molecule,target): start a fresh background compute. Reap the prior worker
    // first (it has finished for any superseded key once dockComputing_ is cleared).
    if (!dockComputing_ && dockThread_.joinable()) dockThread_.join();
    if (dockComputing_) {
        // A previous compute is still running for an older key; let it finish (it will
        // store under its own key and we will relaunch next frame). Avoid piling up.
        computing = true;
        return DockJobResult{};
    }

    dockKey_ = key;
    dockReady_ = false;
    dockComputing_ = true;
    IDockingModule* dock = svc_.docking;
    dockThread_ = std::thread([this, dock, m, target, key]() {
        DockJobResult r = dock ? dock->dockDetailed(m, target) : DockJobResult{};
        std::lock_guard<std::mutex> lk(dockMu_);
        if (dockKey_ == key) {  // still the active selection
            dockResult_ = std::move(r);
            dockReady_ = true;
        }
        dockComputing_ = false;
    });
    computing = true;
    return DockJobResult{};
}

void AppShell::invalidateDock() {
    std::lock_guard<std::mutex> lk(dockMu_);
    // Clearing the key invalidates any in-flight worker (its key check will fail and
    // discard the stale result) and forces a fresh compute next frame.
    dockKey_.clear();
    dockReady_ = false;
}

docking::Provisioner& AppShell::provisioner() {
    if (!provisionProbed_) {
        provisionProbed_ = true;
        // Cheap locate-only probe (no network) to seed the initial status/flags.
        provisioner_->start(/*allowDownload=*/false, docking::headlinePresets());
    }
    return *provisioner_;
}

void AppShell::provisionDocking() {
    provisionProbed_ = true;  // an explicit provision supersedes the locate-only probe
    provisioner_->start(/*allowDownload=*/true, docking::headlinePresets());
}

void AppShell::setRenderDevice(ID3D11Device* device, ID3D11DeviceContext* context) {
    renderDev_ = device;
    renderCtx_ = context;
    // Drop any viewport built against an old device; it is rebuilt lazily.
    viewer_.reset();
}

render::MolViewport* AppShell::viewer() {
    if (!renderDev_ || !renderCtx_) return nullptr;
    if (!viewer_) {
        viewer_ = std::make_unique<render::MolViewport>(renderDev_, renderCtx_);
        if (!viewer_->valid()) {
            viewer_.reset();
            return nullptr;
        }
    }
    return viewer_.get();
}

Molecule AppShell::currentMolecule() const {
    if (state_.hasCustom && state_.selectedMolecule == "__custom__") return state_.customMolecule;
    if (svc_.library) {
        if (auto m = svc_.library->byId(state_.selectedMolecule)) return *m;
        auto all = svc_.library->all();
        if (!all.empty()) return all.front();
    }
    Molecule fallback;
    fallback.name = "(no compound)";
    return fallback;
}

void AppShell::requestHighlight(const std::string& panelId, const std::string& explanation) {
    state_.highlight = panelId;
    state_.highlightStart = ImGui::GetTime();
    for (const auto& p : panels_) {
        if (p.id == panelId) { state_.activePanel = panelId; break; }
    }
    if (!explanation.empty()) state_.assistantLog.push_back(explanation);
}

bool AppShell::isHighlighted(const std::string& panelId) const {
    if (state_.highlight != panelId) return false;
    return (ImGui::GetTime() - state_.highlightStart) < kHighlightSeconds;
}

void AppShell::frameHighlightCurrentWindow(const std::string& panelId) {
    if (!isHighlighted(panelId)) return;
    const ImVec2 pos = ImGui::GetWindowPos();
    const ImVec2 size = ImGui::GetWindowSize();
    pulseBorder(pos, ImVec2(pos.x + size.x, pos.y + size.y),
                ImGui::GetTime() - state_.highlightStart);
}

// ---------------------------------------------------------------- agent wiring
std::string AppShell::buildSystemPrompt() const {
    std::string s = agent::safetySystemPrompt();
    s += "\n\nAvailable panels (use these exact ids with navigate_ui / highlight_panel):\n";
    for (const auto& p : panels_) {
        s += "- " + p.id + " (" + p.label + "): " + p.help + "\n";
    }
    s += "\nThe currently selected compound's real, structure-derived properties are available via "
         "the get_active_compound tool - use it instead of guessing values.";
    return s;
}

bool AppShell::isValidPanel(const std::string& panelId) const {
    for (const auto& p : panels_)
        if (p.id == panelId) return true;
    return false;
}

void AppShell::postHighlightAction(const std::string& panelId, const std::string& explanation) {
    std::lock_guard<std::mutex> lk(agentMu_);
    agentInbox_.push_back({false, panelId, explanation});
}

void AppShell::postNavigateAction(const std::string& panelId) {
    std::lock_guard<std::mutex> lk(agentMu_);
    agentInbox_.push_back({true, panelId, ""});
}

Molecule AppShell::agentCompoundSnapshot() const {
    std::lock_guard<std::mutex> lk(agentMu_);
    return agentCompound_;
}

void AppShell::drainAgentActions() {
    std::vector<AgentUiAction> actions;
    {
        std::lock_guard<std::mutex> lk(agentMu_);
        actions.swap(agentInbox_);
        agentCompound_ = currentMolecule();  // publish for the worker-thread tools
    }
    for (const auto& a : actions) {
        if (!isValidPanel(a.panel)) continue;
        if (a.navigate) {
            state_.activePanel = a.panel;
        } else {
            // Pulse + focus; the assistant's own text carries the explanation, so
            // we pass an empty string to avoid double-logging into assistantLog.
            requestHighlight(a.panel, "");
        }
    }
}

void AppShell::buildAgent() {
    registry_ = std::make_unique<agent::ToolRegistry>();
    mock_ = std::make_unique<agent::MockProvider>();
    anthropic_ = std::make_unique<agent::AnthropicProvider>();
    agent_ = std::make_unique<agent::Agent>();

    // Enum of valid panel ids for the navigation tool schemas.
    nlohmann::json panelEnum = nlohmann::json::array();
    for (const auto& p : panels_) panelEnum.push_back(p.id);

    auto labelOf = [this](const std::string& id) -> std::string {
        for (const auto& p : panels_)
            if (p.id == id) return p.label;
        return id;
    };

    using agent::FunctionTool;
    using nlohmann::json;

    // highlight_panel: pulse + focus a panel and explain why.
    {
        json schema = {
            {"type", "object"},
            {"properties",
             {{"panel", {{"type", "string"}, {"enum", panelEnum},
                         {"description", "Panel id to focus and pulse."}}},
              {"explanation",
               {{"type", "string"}, {"description", "One line on what the user will find there."}}}}},
            {"required", json::array({"panel"})}};
        registry_->add(std::make_unique<FunctionTool>(
            "highlight_panel",
            "Focus a panel and pulse a highlight around it so the user sees where to look. Use this "
            "when answering 'where / how do I ...' questions about the UI.",
            schema, [this, labelOf](const json& args) -> ToolResult {
                const std::string panel = args.value("panel", "");
                if (!isValidPanel(panel)) return {"Unknown panel id: '" + panel + "'.", true};
                postHighlightAction(panel, args.value("explanation", ""));
                return {"Focused and highlighted the " + labelOf(panel) + " panel.", false};
            }));
    }

    // navigate_ui: switch panel without the pulse.
    {
        json schema = {{"type", "object"},
                       {"properties",
                        {{"panel", {{"type", "string"}, {"enum", panelEnum},
                                    {"description", "Panel id to switch to."}}}}},
                       {"required", json::array({"panel"})}};
        registry_->add(std::make_unique<FunctionTool>(
            "navigate_ui",
            "Switch the workspace to a panel (no pulse). Use when the user wants to go somewhere.",
            schema, [this, labelOf](const json& args) -> ToolResult {
                const std::string panel = args.value("panel", "");
                if (!isValidPanel(panel)) return {"Unknown panel id: '" + panel + "'.", true};
                postNavigateAction(panel);
                return {"Switched to the " + labelOf(panel) + " panel.", false};
            }));
    }

    // list_panels: ground the model in the available panels.
    {
        json schema = {{"type", "object"}, {"properties", json::object()}};
        registry_->add(std::make_unique<FunctionTool>(
            "list_panels", "List every StimLab panel with its id and what it shows.", schema,
            [this](const json&) -> ToolResult {
                json arr = json::array();
                for (const auto& p : panels_)
                    arr.push_back({{"id", p.id}, {"label", p.label}, {"shows", p.help}});
                return {arr.dump(), false};
            }));
    }

    // get_active_compound: read the selected compound's real properties (in scope:
    // identity + physicochemical descriptors; NO synthesis data exists to return).
    {
        json schema = {{"type", "object"}, {"properties", json::object()}};
        registry_->add(std::make_unique<FunctionTool>(
            "get_active_compound",
            "Get the currently selected compound's real, structure-derived properties (name, SMILES, "
            "formula, MW, logP, TPSA, H-bond donors/acceptors, drug class, legal status).",
            schema, [this](const json&) -> ToolResult {
                const Molecule m = agentCompoundSnapshot();
                json j = {{"name", m.name},          {"smiles", m.smiles},
                          {"formula", m.formula},    {"molWeight", m.molWeight},
                          {"logP", m.logP},          {"tpsa", m.tpsa},
                          {"hbd", m.hbd},            {"hba", m.hba},
                          {"drugClass", m.drugClass},{"legalStatus", m.legalStatus}};
                return {j.dump(), false};
            }));
    }

    // what_can_stimlab_do: the capability / scope summary.
    {
        json schema = {{"type", "object"}, {"properties", json::object()}};
        registry_->add(std::make_unique<FunctionTool>(
            "what_can_stimlab_do",
            "Summarize what StimLab can and cannot do (its capabilities and its safety scope).",
            schema, [](const json&) -> ToolResult {
                return {"StimLab predicts what a CNS-stimulant compound IS and DOES: structure and "
                        "physicochemical properties, molecular stability, absorption/PK, "
                        "ADMET/metabolism, target binding affinity (docking), similarity to known "
                        "substances, and legal-analog scoring. It does NOT and will not provide "
                        "synthesis routes, reaction conditions, precursors, or manufacturability "
                        "guidance - that is out of scope by design.",
                        false};
            }));
    }

    agent_->configure(mock_.get(), registry_.get(), buildSystemPrompt());
    agent_->setModel(kDefaultModel);
    agentUsingAnthropic_ = false;
}

void AppShell::reconfigureAgent() {
    if (!agent_) return;
    std::string model = kDefaultModel;
    int providerIdx = 0;
    std::string modeStr = "autopilot";
    std::string keyBlob;
    if (config_) {
        model = config_->get<std::string>(kCfgModel, kDefaultModel);
        providerIdx = config_->get<int>(kCfgProvider, 0);
        modeStr = config_->get<std::string>(kCfgMode, "autopilot");
        keyBlob = config_->get<std::string>(kCfgApiKey, "");
    }
    if (model.empty()) model = kDefaultModel;

    // Decrypt the stored key (DPAPI) into the Anthropic provider.
    if (anthropic_) {
        if (!keyBlob.empty()) {
            if (auto plain = Secrets::unprotect(keyBlob); plain.ok())
                anthropic_->setApiKey(plain.value());
            else
                anthropic_->setApiKey("");
        } else {
            anthropic_->setApiKey("");
        }
    }

    const ILlmProvider* active = mock_.get();
    if (providerIdx == 0 && anthropic_ && anthropic_->ready()) {
        active = anthropic_.get();
        agentUsingAnthropic_ = true;
    } else {
        agentUsingAnthropic_ = false;
    }

    agent_->setProvider(active);
    agent_->setModel(model);
    agent_->setMode(modeStr == "askfirst" ? agent::AgentMode::AskFirst
                                          : agent::AgentMode::Autopilot);
}

void AppShell::setConfig(Config* config) {
    config_ = config;
    reconfigureAgent();
}

// ---- Settings helpers ----
bool AppShell::hasApiKey() const {
    return config_ && config_->has(kCfgApiKey) &&
           !config_->get<std::string>(kCfgApiKey, "").empty();
}

void AppShell::saveApiKey(const std::string& plaintext) {
    if (!config_) return;
    if (plaintext.empty()) { clearApiKey(); return; }
    if (auto enc = Secrets::protect(plaintext); enc.ok()) {
        config_->set(kCfgApiKey, enc.value());
        config_->set(kCfgProvider, 0);  // entering a key means "use Anthropic"
        config_->save();
        reconfigureAgent();
    }
}

void AppShell::clearApiKey() {
    if (!config_) return;
    config_->set(kCfgApiKey, std::string(""));
    config_->save();
    reconfigureAgent();
}

void AppShell::setAgentModel(const std::string& model) {
    if (config_) { config_->set(kCfgModel, model); config_->save(); }
    if (agent_) agent_->setModel(model.empty() ? kDefaultModel : model);
}

std::string AppShell::agentModel() const {
    return config_ ? config_->get<std::string>(kCfgModel, kDefaultModel) : std::string(kDefaultModel);
}

void AppShell::setAgentProviderIndex(int idx) {
    if (config_) { config_->set(kCfgProvider, idx); config_->save(); }
    reconfigureAgent();
}

int AppShell::agentProviderIndex() const {
    return config_ ? config_->get<int>(kCfgProvider, 0) : 0;
}

void AppShell::setAutopilot(bool on) {
    if (config_) { config_->set(kCfgMode, std::string(on ? "autopilot" : "askfirst")); config_->save(); }
    if (agent_) agent_->setMode(on ? agent::AgentMode::Autopilot : agent::AgentMode::AskFirst);
}

bool AppShell::autopilot() const {
    return !config_ || config_->get<std::string>(kCfgMode, "autopilot") != "askfirst";
}

bool AppShell::anthropicReady() const { return anthropic_ && anthropic_->ready(); }
bool AppShell::anthropicTransport() const { return agent::AnthropicProvider::transportAvailable(); }

std::string AppShell::activeProviderLabel() const {
    return agentUsingAnthropic_ ? std::string("Anthropic (Claude)")
                                : (mock_ ? mock_->displayName() : std::string("Offline"));
}

void AppShell::draw() {
    drainAgentActions();
    if (agent_) agent_->poll();
    drawMainMenuBar();

    const ImGuiViewport* vp = ImGui::GetMainViewport();
    const ImVec2 wp = vp->WorkPos;   // already below the main menu bar
    const ImVec2 ws = vp->WorkSize;
    const float navW = 248.0f;
    const float asstW = state_.showAssistant ? 372.0f : 0.0f;
    const float contentW = ws.x - navW - asstW;

    ImGui::SetNextWindowPos(wp, ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(navW, ws.y), ImGuiCond_Always);
    drawNavigator();

    ImGui::SetNextWindowPos(ImVec2(wp.x + navW, wp.y), ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(contentW, ws.y), ImGuiCond_Always);
    drawContent();

    if (state_.showAssistant) {
        ImGui::SetNextWindowPos(ImVec2(wp.x + navW + contentW, wp.y), ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2(asstW, ws.y), ImGuiCond_Always);
        drawAssistant();
    }

    drawAboutModal();
}

void AppShell::drawMainMenuBar() {
    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("Exit")) state_.quitRequested = true;
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("View")) {
            ImGui::MenuItem("Assistant panel", nullptr, &state_.showAssistant);
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Help")) {
            if (ImGui::MenuItem("About StimLab")) state_.showAbout = true;
            ImGui::EndMenu();
        }
        ImGui::Separator();
        ImGui::TextDisabled("StimLab  -  CNS-stimulant analysis suite  -  real chem engine + 3D + docking");
        ImGui::EndMainMenuBar();
    }
}

void AppShell::drawNavigator() {
    if (ImGui::Begin("Navigator", nullptr, kPaneFlags)) {
        ImGui::PushStyleColor(ImGuiCol_Text, theme::verdictColor(1));
        ImGui::SetWindowFontScale(1.25f);
        ImGui::TextUnformatted("StimLab");
        ImGui::SetWindowFontScale(1.0f);
        ImGui::PopStyleColor();
        ImGui::TextDisabled("Stimulant Laboratory");
        ImGui::Spacing();

        ImGui::TextDisabled("ACTIVE COMPOUND");
        if (svc_.library) {
            const Molecule cur = currentMolecule();
            if (ImGui::BeginCombo("##compound", cur.name.c_str())) {
                for (const auto& m : svc_.library->all()) {
                    const bool sel = (m.id == state_.selectedMolecule);
                    if (ImGui::Selectable(m.name.c_str(), sel)) state_.selectedMolecule = m.id;
                    if (sel) ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }
        }
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::TextDisabled("WORKSPACE");
        ImGui::Spacing();

        for (const auto& p : panels_) {
            const bool selected = (state_.activePanel == p.id);
            if (ImGui::Selectable(p.label.c_str(), selected, 0, ImVec2(0, 26))) {
                state_.activePanel = p.id;
            }
            if (isHighlighted(p.id)) {
                pulseBorder(ImGui::GetItemRectMin(), ImGui::GetItemRectMax(),
                            ImGui::GetTime() - state_.highlightStart);
            }
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", p.help.c_str());
        }

        ImGui::Spacing();
        ImGui::Separator();
        if (!state_.showAssistant && ImGui::Button("Open Assistant", ImVec2(-1, 0))) {
            state_.showAssistant = true;
        }
    }
    frameHighlightCurrentWindow("Navigator");
    ImGui::End();
}

void AppShell::drawContent() {
    if (ImGui::Begin("Workspace", nullptr, kPaneFlags)) {
        const Molecule cur = currentMolecule();
        std::string label = "Dashboard";
        std::string help;
        for (const auto& p : panels_) {
            if (p.id == state_.activePanel) { label = p.label; help = p.help; break; }
        }

        ImGui::SetWindowFontScale(1.35f);
        ImGui::TextUnformatted(label.c_str());
        ImGui::SetWindowFontScale(1.0f);
        ImGui::SameLine();
        ImGui::TextDisabled("   *  %s", cur.name.c_str());
        if (!help.empty()) ImGui::TextDisabled("%s", help.c_str());
        ImGui::Separator();
        ImGui::Spacing();

        const std::string& panel = state_.activePanel;
        if (panel == "Dashboard")        panels::dashboard(*this);
        else if (panel == "Structure")   panels::structureWorkbench(*this);
        else if (panel == "Input")       panels::moleculeInput(*this);
        else if (panel == "Analog")      panels::analogExplorer(*this);
        else if (panel == "Compare")     panels::compare(*this);
        else if (panel == "Stability")   panels::stability(*this);
        else if (panel == "Absorption")  panels::absorption(*this);
        else if (panel == "Metabolism")  panels::metabolism(*this);
        else if (panel == "Similarity")  panels::similarity(*this);
        else if (panel == "Legal")       panels::legal(*this);
        else if (panel == "Docking")     panels::docking(*this);
        else if (panel == "Workflows")   panels::workflows(*this);
        else if (panel == "Library")     panels::library(*this);
        else if (panel == "Runs")        panels::runs(*this);
        else if (panel == "Presets")     panels::presets(*this);
        else if (panel == "Settings")    panels::settings(*this);
        else                              panels::dashboard(*this);
    }
    frameHighlightCurrentWindow(state_.activePanel);
    ImGui::End();
}

void AppShell::drawAssistant() {
    if (ImGui::Begin("Assistant", nullptr, kPaneFlags)) {
        ImGui::PushStyleColor(ImGuiCol_Text, theme::verdictColor(0));
        ImGui::TextUnformatted("StimLab Assistant");
        ImGui::PopStyleColor();

        const agent::AgentSnapshot snap = agent_ ? agent_->snapshot() : agent::AgentSnapshot{};
        const bool busy = snap.status == agent::AgentStatus::Running ||
                          snap.status == agent::AgentStatus::AwaitingApproval;

        // Provider status line.
        if (agentUsingAnthropic_) {
            ImGui::TextDisabled("Live: %s  -  model %s", activeProviderLabel().c_str(),
                                agentModel().c_str());
        } else if (anthropicTransport()) {
            ImGui::TextDisabled("Offline assistant - add an API key in Settings for live chat.");
        } else {
            ImGui::TextDisabled("Offline assistant (this build has no networking).");
        }

        // Mode + reset row.
        bool autop = autopilot();
        if (ImGui::Checkbox("Autopilot", &autop)) setAutopilot(autop);
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("On: tools (navigate/highlight) run automatically.\n"
                              "Off (ask-first): you approve each tool batch.");
        ImGui::SameLine();
        ImGui::BeginDisabled(busy);
        if (ImGui::SmallButton("New chat") && agent_) agent_->reset();
        ImGui::EndDisabled();

        // Quick prompts (canned fallbacks - now routed through the real loop).
        ImGui::TextDisabled("TRY ASKING");
        ImGui::BeginDisabled(busy);
        struct QP { const char* label; const char* prompt; };
        static const QP kQuick[] = {
            {"How do I change the docking target?", "How do I change the docking target?"},
            {"Where is absorption / bioavailability?", "Where do I see absorption and bioavailability?"},
            {"What can StimLab do?", "What can StimLab do?"},
            {"Tell me about the selected compound", "Tell me about the currently selected compound."},
        };
        for (const auto& q : kQuick) {
            if (ImGui::Button(q.label, ImVec2(-1, 0)) && agent_) agent_->submit(q.prompt);
        }
        ImGui::EndDisabled();

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::TextDisabled("CONVERSATION");

        const float reserve = (snap.status == agent::AgentStatus::AwaitingApproval) ? 132.0f : 40.0f;
        ImGui::BeginChild("##log", ImVec2(0, -reserve), ImGuiChildFlags_Borders);
        if (snap.transcript.empty() && snap.streaming.empty()) {
            ImGui::TextDisabled("Ask a question or tap one above. The assistant can navigate and "
                                "highlight the UI, and read the selected compound's real properties.");
        }
        for (const auto& e : snap.transcript) {
            switch (e.kind) {
                case agent::TranscriptEntry::Kind::User:
                    ImGui::PushStyleColor(ImGuiCol_Text, theme::verdictColor(1));
                    ImGui::TextWrapped("You: %s", e.text.c_str());
                    ImGui::PopStyleColor();
                    break;
                case agent::TranscriptEntry::Kind::Assistant:
                    ImGui::TextWrapped("%s", e.text.c_str());
                    break;
                case agent::TranscriptEntry::Kind::Tool:
                    ImGui::TextDisabled("  %s", e.text.c_str());
                    break;
                case agent::TranscriptEntry::Kind::Error:
                    ImGui::PushStyleColor(ImGuiCol_Text, theme::verdictColor(3));
                    ImGui::TextWrapped("%s", e.text.c_str());
                    ImGui::PopStyleColor();
                    break;
                case agent::TranscriptEntry::Kind::System:
                    ImGui::TextDisabled("%s", e.text.c_str());
                    break;
            }
            ImGui::Spacing();
        }
        if (!snap.streaming.empty()) ImGui::TextWrapped("%s", snap.streaming.c_str());
        if (busy) {
            ImGui::TextDisabled("...");
            ImGui::SetScrollHereY(1.0f);
        }
        ImGui::EndChild();

        // Ask-first approval gate.
        if (snap.status == agent::AgentStatus::AwaitingApproval) {
            ImGui::PushStyleColor(ImGuiCol_Text, theme::verdictColor(2));
            ImGui::TextWrapped("The assistant wants to run:");
            ImGui::PopStyleColor();
            for (const auto& p : snap.pending)
                ImGui::BulletText("%s %s", p.name.c_str(), p.arguments.dump().c_str());
            if (ImGui::Button("Approve", ImVec2(120, 0)) && agent_) agent_->approvePending();
            ImGui::SameLine();
            if (ImGui::Button("Deny", ImVec2(120, 0)) && agent_) agent_->denyPending();
        }

        // Input row.
        ImGui::BeginDisabled(busy);
        ImGui::SetNextItemWidth(-1);
        const bool submitted = ImGui::InputTextWithHint(
            "##chat", "Ask about a panel, the compound, or how to do something...", chatBuf_,
            sizeof(chatBuf_), ImGuiInputTextFlags_EnterReturnsTrue);
        ImGui::EndDisabled();
        if (submitted && !busy && agent_ && chatBuf_[0] != '\0') {
            agent_->submit(chatBuf_);
            chatBuf_[0] = '\0';
            ImGui::SetKeyboardFocusHere(-1);
        }
    }
    frameHighlightCurrentWindow("Assistant");
    ImGui::End();
}

void AppShell::drawAboutModal() {
    if (state_.showAbout) {
        ImGui::OpenPopup("About StimLab");
        state_.showAbout = false;
    }
    ImGui::SetNextWindowSize(ImVec2(520, 0), ImGuiCond_Appearing);
    if (ImGui::BeginPopupModal("About StimLab", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::TextUnformatted("StimLab - native CNS-stimulant analysis suite");
        ImGui::TextDisabled("In-house cheminformatics engine - 3D viewer, real docking backends, SQLite.");
        ImGui::Separator();
        ImGui::TextWrapped(
            "IN SCOPE: structure/properties, molecular stability, absorption/PK, ADMET/metabolism, "
            "target binding affinity (docking), similarity to known substances, legal-analog scoring, "
            "and an assistant that explains and highlights the UI.");
        ImGui::Spacing();
        ImGui::PushStyleColor(ImGuiCol_Text, theme::verdictColor(2));
        ImGui::TextWrapped(
            "OUT OF SCOPE (by design): synthesis routes/steps, reaction conditions, precursor "
            "selection, and any manufacturability / ease-of-manufacture scoring.");
        ImGui::PopStyleColor();
        ImGui::Spacing();
        if (ImGui::Button("Close", ImVec2(120, 0))) ImGui::CloseCurrentPopup();
        ImGui::EndPopup();
    }
}

}  // namespace stimlab
