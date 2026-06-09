#include "ui/AppShell.h"

#include <cmath>

#include <imgui.h>

#include "render/MolViewport.h"
#include "ui/Panels.h"
#include "ui/Theme.h"

namespace stimlab {

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
         "Identity + physicochemical properties of the selected molecule (2D depiction lands in Phase C)."},
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
        {"Library", "Library",
         "Browse, search and select compounds from the default + imported library."},
        {"Runs", "Runs",
         "History of analyses with status and summaries."},
        {"Presets", "Presets / Targets",
         "CNS target presets and reusable analysis panels."},
        {"Settings", "Settings",
         "AI provider/keys, GPU mode, storage paths."},
    };
}

AppShell::~AppShell() = default;  // here MolViewport is a complete type

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

void AppShell::draw() {
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
        ImGui::TextDisabled("StimLab  -  CNS-stimulant analysis suite  -  Phase B skeleton");
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
        else if (panel == "Analog")      panels::analogExplorer(*this);
        else if (panel == "Compare")     panels::compare(*this);
        else if (panel == "Stability")   panels::stability(*this);
        else if (panel == "Absorption")  panels::absorption(*this);
        else if (panel == "Metabolism")  panels::metabolism(*this);
        else if (panel == "Similarity")  panels::similarity(*this);
        else if (panel == "Legal")       panels::legal(*this);
        else if (panel == "Docking")     panels::docking(*this);
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
        ImGui::TextWrapped(
            "I can explain any panel and highlight where to click. Wire a provider in "
            "Settings (Phase D) for full conversational control.");
        ImGui::Spacing();
        ImGui::TextDisabled("TRY ASKING");

        if (ImGui::Button("How do I change the docking target?", ImVec2(-1, 0))) {
            requestHighlight("Docking",
                "Open the Docking panel - the target dropdown at the top selects DAT/NET/SERT/TAAR1. "
                "I've highlighted it for you.");
        }
        if (ImGui::Button("Where do I see absorption / bioavailability?", ImVec2(-1, 0))) {
            requestHighlight("Absorption",
                "The Absorption / PK panel shows HIA, oral F%, Caco-2 permeability, BBB partition and "
                "P-gp efflux. Highlighted now.");
        }
        if (ImGui::Button("How is stability scored (not manufacturability)?", ImVec2(-1, 0))) {
            requestHighlight("Stability",
                "Stability scores resistance to hydrolysis/oxidation/photolysis/thermal/pH stress and "
                "estimates shelf-life. It deliberately REPLACES any manufacturability score - out of scope.");
        }
        if (ImGui::Button("How do I pick a compound?", ImVec2(-1, 0))) {
            requestHighlight("Navigator",
                "Use the ACTIVE COMPOUND dropdown at the top of the Navigator (highlighted) to switch "
                "the molecule every panel analyzes.");
        }
        if (ImGui::Button("How do I check a derivative vs known samples?", ImVec2(-1, 0))) {
            requestHighlight("Analog",
                "Open Analog Explorer (highlighted): tune a candidate's properties + functional groups "
                "and instantly see its nearest existing sample, legal-analog score, and predicted "
                "byproducts/interactions - all analysis, no synthesis guidance.");
        }
        if (ImGui::Button("Can I compare compounds side by side?", ImVec2(-1, 0))) {
            requestHighlight("Compare",
                "Yes - the Compare panel (highlighted) charts stability, oral F%, HIA and oxidation "
                "across up to three compounds, with a full metric matrix below.");
        }
        if (ImGui::Button("What can StimLab do?", ImVec2(-1, 0))) {
            state_.assistantLog.push_back(
                "StimLab predicts what a compound IS and DOES: structure/properties, stability, "
                "absorption/PK, ADMET/metabolism, target binding (docking), similarity to known "
                "substances, and legal-analog scoring. It does NOT provide synthesis routes.");
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::TextDisabled("CONVERSATION");
        ImGui::BeginChild("##log", ImVec2(0, -38), ImGuiChildFlags_Borders);
        if (state_.assistantLog.empty()) {
            ImGui::TextDisabled("(Tap a question above. Responses appear here.)");
        }
        for (const auto& line : state_.assistantLog) {
            ImGui::TextWrapped("- %s", line.c_str());
            ImGui::Spacing();
        }
        ImGui::EndChild();

        ImGui::BeginDisabled(true);
        char buf[8] = {0};
        ImGui::SetNextItemWidth(-1);
        ImGui::InputTextWithHint("##chat", "Connect a provider in Settings to chat...", buf,
                                 sizeof(buf));
        ImGui::EndDisabled();
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
        ImGui::TextDisabled("Phase B skeleton - thick fakes behind a clean DX11/ImGui GUI.");
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
