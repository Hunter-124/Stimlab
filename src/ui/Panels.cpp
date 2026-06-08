#include "ui/Panels.h"

#include <algorithm>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include <imgui.h>
#include <implot.h>

#include "core/AppPaths.h"
#include "ui/AppShell.h"
#include "ui/Theme.h"

namespace stimlab::panels {
namespace {

std::string f2(double v) { char b[40]; std::snprintf(b, sizeof b, "%.2f", v); return b; }
std::string f0(double v) { char b[40]; std::snprintf(b, sizeof b, "%.0f", v); return b; }

void verdictText(Verdict v) {
    ImGui::TextColored(theme::verdictColor(static_cast<int>(v)), "%s", verdictLabel(v));
}

void statCard(const char* title, const std::string& value, const char* sub, float w = 168.0f) {
    ImGui::BeginChild(title, ImVec2(w, 84), ImGuiChildFlags_Borders);
    ImGui::TextDisabled("%s", title);
    ImGui::SetWindowFontScale(1.6f);
    ImGui::TextUnformatted(value.c_str());
    ImGui::SetWindowFontScale(1.0f);
    if (sub) ImGui::TextDisabled("%s", sub);
    ImGui::EndChild();
}

// A small decorative 2D schematic (real RDKit depiction arrives in Phase C).
void moleculeSchematic(const Molecule& m, ImVec2 size) {
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const ImVec2 o = ImGui::GetCursorScreenPos();
    ImGui::Dummy(size);

    const ImU32 bond = IM_COL32(150, 160, 175, 255);
    const ImU32 carbon = IM_COL32(70, 80, 96, 255);
    const ImU32 nitrogen = IM_COL32(96, 165, 250, 255);
    const ImU32 oxygen = IM_COL32(248, 113, 113, 255);

    const ImVec2 c(o.x + size.x * 0.32f, o.y + size.y * 0.5f);
    const float r = std::min(size.x, size.y) * 0.22f;

    ImVec2 hex[6];
    for (int i = 0; i < 6; ++i) {
        const float a = static_cast<float>(i) / 6.0f * 6.2831853f;
        hex[i] = ImVec2(c.x + r * std::cos(a), c.y + r * std::sin(a));
    }
    for (int i = 0; i < 6; ++i) dl->AddLine(hex[i], hex[(i + 1) % 6], bond, 2.0f);
    dl->AddCircle(c, r * 0.55f, bond, 24, 1.4f);  // aromatic ring marker
    for (int i = 0; i < 6; ++i) dl->AddCircleFilled(hex[i], 4.0f, carbon);

    // Side chain off the rightmost vertex -> amine (phenethylamine motif).
    ImVec2 a1 = ImVec2(hex[0].x + r * 0.9f, hex[0].y - r * 0.4f);
    ImVec2 a2 = ImVec2(a1.x + r * 0.8f, a1.y + r * 0.5f);
    ImVec2 nAtom = ImVec2(a2.x + r * 0.8f, a2.y - r * 0.4f);
    dl->AddLine(hex[0], a1, bond, 2.0f);
    dl->AddLine(a1, a2, bond, 2.0f);
    dl->AddLine(a2, nAtom, bond, 2.0f);
    dl->AddCircleFilled(a1, 4.0f, carbon);
    dl->AddCircleFilled(a2, 4.0f, carbon);
    dl->AddCircleFilled(nAtom, 6.0f, nitrogen);
    dl->AddText(ImVec2(nAtom.x + 6, nAtom.y - 7), nitrogen, "N");

    if (m.smiles.find("O") != std::string::npos) {
        ImVec2 oAtom = ImVec2(hex[3].x - r * 0.8f, hex[3].y + r * 0.4f);
        dl->AddLine(hex[3], oAtom, bond, 2.0f);
        dl->AddCircleFilled(oAtom, 6.0f, oxygen);
        dl->AddText(ImVec2(oAtom.x - 16, oAtom.y - 7), oxygen, "O");
    }
}

// Build a full Markdown analysis report for one compound (identity -> legal).
std::string buildReportMarkdown(Services& s, const Molecule& m) {
    std::string o;
    auto line = [&](const std::string& x) { o += x; o += "\n"; };
    line("# StimLab report - " + m.name);
    line("");
    line("## Identity");
    line("- SMILES: " + m.smiles);
    line("- Formula: " + m.formula);
    line("- Class: " + m.drugClass);
    line("- Legal status: " + m.legalStatus);
    line("");
    line("## Physicochemical");
    line("- MW " + f2(m.molWeight) + " | logP " + f2(m.logP) + " | TPSA " + f2(m.tpsa));
    line("- HBD " + std::to_string(m.hbd) + " | HBA " + std::to_string(m.hba) +
         " | RotB " + std::to_string(m.rotatableBonds));
    if (s.stability) {
        const auto r = s.stability->analyze(m);
        line("");
        line("## Stability: " + f0(r.overallScore) + "/100 (" + r.shelfLifeEstimate + ")");
        for (const auto& f : r.factors) line("- " + f.name + ": " + f0(f.score) + " - " + f.rationale);
        for (const auto& d : r.degradants) line("- Degradant: " + d.name + " (" + d.pathway + ")");
    }
    if (s.absorption) {
        const auto r = s.absorption->predict(m);
        line("");
        line("## Absorption / PK");
        for (const auto& mt : r.metrics) line("- " + mt.name + ": " + f2(mt.value) + " " + mt.unit);
    }
    if (s.admet) {
        const auto r = s.admet->screen(m);
        line("");
        line(std::string("## ADMET: ") + verdictLabel(r.overall));
        for (const auto& e : r.endpoints)
            line("- [" + std::string(verdictLabel(e.verdict)) + "] " + e.name + " - " + e.detail);
    }
    if (s.similarity) {
        const auto r = s.similarity->search(m);
        line("");
        line("## Similarity (nearest: " + r.nearestName + ")");
        for (size_t i = 0; i < r.hits.size() && i < 5; ++i)
            line("- " + r.hits[i].referenceName + ": Tanimoto " + f2(r.hits[i].tanimoto));
    }
    if (s.legal) {
        const auto r = s.legal->score(m);
        line("");
        line("## Legal-analog: " + f0(r.substantialSimilarity) + "/100 - " + r.classification);
    }
    line("");
    line("_Analysis only. No synthesis / manufacturability content (out of scope by design)._");
    return o;
}

}  // namespace

// ------------------------------------------------------------------ Dashboard
void dashboard(AppShell& shell) {
    Services& s = shell.services();
    const Molecule m = shell.currentMolecule();

    statCard("LIBRARY", std::to_string(s.library ? s.library->count() : 0), "compounds");
    ImGui::SameLine();
    statCard("RECENT RUNS", std::to_string(s.runs ? s.runs->recent().size() : 0), "this session");
    ImGui::SameLine();
    statCard("ACTIVE", m.name, m.drugClass.c_str(), 240.0f);
    ImGui::Spacing();

    ImGui::TextDisabled("SNAPSHOT  -  %s", m.name.c_str());
    ImGui::Separator();

    if (s.stability && s.absorption && s.admet) {
        const auto stab = s.stability->analyze(m);
        const auto abs = s.absorption->predict(m);
        const auto adm = s.admet->screen(m);

        statCard("STABILITY", f0(stab.overallScore) + "/100", stab.shelfLifeEstimate.c_str(), 200.0f);
        ImGui::SameLine();
        statCard("ORAL F", f0(abs.bioavailabilityPct) + "%",
                 abs.cnsPenetrant ? "CNS-penetrant" : "low CNS", 168.0f);
        ImGui::SameLine();
        statCard("HIA", f0(abs.hiaPct) + "%", "intestinal abs.", 168.0f);
        ImGui::SameLine();
        {
            ImGui::BeginChild("ADMETcard", ImVec2(200, 84), ImGuiChildFlags_Borders);
            ImGui::TextDisabled("ADMET");
            ImGui::SetWindowFontScale(1.4f);
            verdictText(adm.overall);
            ImGui::SetWindowFontScale(1.0f);
            ImGui::TextDisabled("%zu endpoints", adm.endpoints.size());
            ImGui::EndChild();
        }
        ImGui::Spacing();
        ImGui::TextWrapped("%s", stab.summary.c_str());
        ImGui::TextWrapped("%s", abs.summary.c_str());
    }
}

// ---------------------------------------------------------- Structure Workbench
void structureWorkbench(AppShell& shell) {
    const Molecule m = shell.currentMolecule();

    if (ImGui::BeginTable("idprop", 2, ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_Resizable)) {
        ImGui::TableSetupColumn("Identity", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("2D schematic", ImGuiTableColumnFlags_WidthFixed, 280.0f);
        ImGui::TableNextRow();

        ImGui::TableSetColumnIndex(0);
        if (ImGui::BeginTable("ident", 2, ImGuiTableFlags_BordersInnerH)) {
            auto row = [](const char* k, const std::string& v) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::TextDisabled("%s", k);
                ImGui::TableSetColumnIndex(1); ImGui::TextWrapped("%s", v.c_str());
            };
            row("Name", m.name);
            row("Formula", m.formula);
            row("SMILES", m.smiles);
            row("Class", m.drugClass);
            row("Legal status", m.legalStatus);
            row("Notes", m.notes);
            ImGui::EndTable();
        }
        ImGui::Spacing();
        ImGui::TextDisabled("PHYSICOCHEMICAL");
        if (ImGui::BeginTable("props", 6, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
            for (const char* h : {"MW", "logP", "TPSA", "HBD", "HBA", "RotB"})
                ImGui::TableSetupColumn(h);
            ImGui::TableHeadersRow();
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(f2(m.molWeight).c_str());
            ImGui::TableSetColumnIndex(1); ImGui::TextUnformatted(f2(m.logP).c_str());
            ImGui::TableSetColumnIndex(2); ImGui::TextUnformatted(f2(m.tpsa).c_str());
            ImGui::TableSetColumnIndex(3); ImGui::Text("%d", m.hbd);
            ImGui::TableSetColumnIndex(4); ImGui::Text("%d", m.hba);
            ImGui::TableSetColumnIndex(5); ImGui::Text("%d", m.rotatableBonds);
            ImGui::EndTable();
        }

        ImGui::TableSetColumnIndex(1);
        ImGui::TextDisabled("Schematic (illustrative)");
        moleculeSchematic(m, ImVec2(260, 220));
        ImGui::TextDisabled("Real 2D/3D depiction lands with RDKit (Phase C).");
        ImGui::EndTable();
    }

    ImGui::Spacing();
    static std::string exportStatus;
    if (ImGui::Button("Export full report (.md)")) {
        const auto path = AppPaths::instance().root() / ("report-" + m.id + ".md");
        std::ofstream out(path);
        if (out) {
            out << buildReportMarkdown(shell.services(), m);
            exportStatus = "Saved: " + path.string();
        } else {
            exportStatus = "Failed to write report.";
        }
    }
    if (!exportStatus.empty()) {
        ImGui::SameLine();
        ImGui::TextDisabled("%s", exportStatus.c_str());
    }
}

// ------------------------------------------------------------------- Stability
void stability(AppShell& shell) {
    const Molecule m = shell.currentMolecule();
    if (!shell.services().stability) return;
    const auto r = shell.services().stability->analyze(m);

    statCard("OVERALL", f0(r.overallScore) + "/100", "higher = more stable", 200.0f);
    ImGui::SameLine();
    statCard("SHELF LIFE", r.shelfLifeEstimate, "estimated", 300.0f);
    ImGui::Spacing();

    std::vector<double> vals;
    std::vector<const char*> labels;
    for (const auto& f : r.factors) { vals.push_back(f.score); labels.push_back(f.name.c_str()); }

    if (ImPlot::BeginPlot("##stabchart", ImVec2(-1, 200), ImPlotFlags_NoMouseText | ImPlotFlags_NoLegend)) {
        std::vector<double> pos;
        for (size_t i = 0; i < labels.size(); ++i) pos.push_back(static_cast<double>(i));
        ImPlot::SetupAxes(nullptr, "score (0-100)", ImPlotAxisFlags_NoGridLines, 0);
        ImPlot::SetupAxesLimits(-0.6, static_cast<double>(labels.size()) - 0.4, 0, 100, ImPlotCond_Always);
        ImPlot::SetupAxisTicks(ImAxis_X1, pos.data(), static_cast<int>(pos.size()), labels.data());
        ImPlot::PlotBars("score", vals.data(), static_cast<int>(vals.size()), 0.6);
        ImPlot::EndPlot();
    }

    if (ImGui::BeginTable("facs", 2, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Factor", ImGuiTableColumnFlags_WidthFixed, 180.0f);
        ImGui::TableSetupColumn("Rationale");
        ImGui::TableHeadersRow();
        for (const auto& f : r.factors) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0); ImGui::Text("%s  (%s)", f.name.c_str(), f0(f.score).c_str());
            ImGui::TableSetColumnIndex(1); ImGui::TextWrapped("%s", f.rationale.c_str());
        }
        ImGui::EndTable();
    }
    ImGui::Spacing();
    ImGui::TextDisabled("LIKELY DEGRADANTS");
    for (const auto& d : r.degradants)
        ImGui::BulletText("%s  -  %s (%s)", d.name.c_str(), d.pathway.c_str(), d.note.c_str());
}

// ------------------------------------------------------------------ Absorption
void absorption(AppShell& shell) {
    const Molecule m = shell.currentMolecule();
    if (!shell.services().absorption) return;
    const auto r = shell.services().absorption->predict(m);

    statCard("ORAL F", f0(r.bioavailabilityPct) + "%", "predicted bioavailability", 190.0f);
    ImGui::SameLine();
    statCard("HIA", f0(r.hiaPct) + "%", "intestinal absorption", 170.0f);
    ImGui::SameLine();
    statCard("logBB", f2(r.logBB), r.cnsPenetrant ? "CNS-penetrant" : "peripheral", 150.0f);
    ImGui::SameLine();
    statCard("P-gp", r.pgpSubstrate ? "substrate" : "no", "efflux", 150.0f);
    ImGui::Spacing();

    // Percent-scale metrics charted together (F and HIA), others in the table.
    if (ImPlot::BeginPlot("##pkchart", ImVec2(-1, 180), ImPlotFlags_NoMouseText | ImPlotFlags_NoLegend)) {
        const char* labels[] = {"Oral F %", "HIA %"};
        const double pos[] = {0, 1};
        const double vals[] = {r.bioavailabilityPct, r.hiaPct};
        ImPlot::SetupAxes(nullptr, "%", ImPlotAxisFlags_NoGridLines, 0);
        ImPlot::SetupAxesLimits(-0.6, 1.6, 0, 100, ImPlotCond_Always);
        ImPlot::SetupAxisTicks(ImAxis_X1, pos, 2, labels);
        ImPlot::PlotBars("pct", vals, 2, 0.55);
        ImPlot::EndPlot();
    }

    if (ImGui::BeginTable("pk", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Metric", ImGuiTableColumnFlags_WidthFixed, 200.0f);
        ImGui::TableSetupColumn("Value", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Unit", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Note");
        ImGui::TableHeadersRow();
        for (const auto& mt : r.metrics) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextColored(theme::verdictColor(static_cast<int>(mt.band)), "%s", mt.name.c_str());
            ImGui::TableSetColumnIndex(1); ImGui::TextUnformatted(f2(mt.value).c_str());
            ImGui::TableSetColumnIndex(2); ImGui::TextUnformatted(mt.unit.c_str());
            ImGui::TableSetColumnIndex(3); ImGui::TextWrapped("%s", mt.rationale.c_str());
        }
        ImGui::EndTable();
    }
    ImGui::Spacing();
    ImGui::TextWrapped("%s", r.summary.c_str());
}

// ------------------------------------------------------------------ Metabolism
void metabolism(AppShell& shell) {
    const Molecule m = shell.currentMolecule();
    if (!shell.services().admet) return;
    const auto r = shell.services().admet->screen(m);

    ImGui::TextUnformatted("Overall verdict:");
    ImGui::SameLine();
    verdictText(r.overall);
    ImGui::TextWrapped("%s", r.summary.c_str());
    ImGui::Separator();

    if (ImGui::BeginTable("endpoints", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Endpoint", ImGuiTableColumnFlags_WidthFixed, 230.0f);
        ImGui::TableSetupColumn("Verdict", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Detail");
        ImGui::TableHeadersRow();
        for (const auto& e : r.endpoints) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(e.name.c_str());
            ImGui::TableSetColumnIndex(1); verdictText(e.verdict);
            ImGui::TableSetColumnIndex(2); ImGui::TextWrapped("%s", e.detail.c_str());
        }
        ImGui::EndTable();
    }
}

// ------------------------------------------------------------------ Similarity
void similarity(AppShell& shell) {
    const Molecule m = shell.currentMolecule();
    if (!shell.services().similarity) return;
    const auto r = shell.services().similarity->search(m);

    statCard("NEAREST", r.nearestName.empty() ? "-" : r.nearestName,
             ("Tanimoto " + f2(r.nearestScore)).c_str(), 280.0f);
    ImGui::Spacing();

    if (!r.hits.empty() &&
        ImPlot::BeginPlot("##simchart", ImVec2(-1, 200), ImPlotFlags_NoMouseText | ImPlotFlags_NoLegend)) {
        std::vector<double> vals, pos;
        std::vector<const char*> labels;
        for (size_t i = 0; i < r.hits.size(); ++i) {
            vals.push_back(r.hits[i].tanimoto);
            pos.push_back(static_cast<double>(i));
            labels.push_back(r.hits[i].referenceName.c_str());
        }
        ImPlot::SetupAxes(nullptr, "Tanimoto", ImPlotAxisFlags_NoGridLines, 0);
        ImPlot::SetupAxesLimits(-0.6, static_cast<double>(r.hits.size()) - 0.4, 0, 1, ImPlotCond_Always);
        ImPlot::SetupAxisTicks(ImAxis_X1, pos.data(), static_cast<int>(pos.size()), labels.data());
        ImPlot::PlotBars("tani", vals.data(), static_cast<int>(vals.size()), 0.6);
        ImPlot::EndPlot();
    }

    if (ImGui::BeginTable("hits", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Reference");
        ImGui::TableSetupColumn("Class");
        ImGui::TableSetupColumn("Legal");
        ImGui::TableSetupColumn("Tanimoto", ImGuiTableColumnFlags_WidthFixed, 80.0f);
        ImGui::TableSetupColumn("Pharmacophore", ImGuiTableColumnFlags_WidthFixed, 110.0f);
        ImGui::TableHeadersRow();
        for (const auto& h : r.hits) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(h.referenceName.c_str());
            ImGui::TableSetColumnIndex(1); ImGui::TextUnformatted(h.referenceClass.c_str());
            ImGui::TableSetColumnIndex(2); ImGui::TextUnformatted(h.legalStatus.c_str());
            ImGui::TableSetColumnIndex(3); ImGui::TextUnformatted(f2(h.tanimoto).c_str());
            ImGui::TableSetColumnIndex(4); ImGui::TextUnformatted(f2(h.pharmacophore).c_str());
        }
        ImGui::EndTable();
    }
}

// ----------------------------------------------------------------------- Legal
void legal(AppShell& shell) {
    const Molecule m = shell.currentMolecule();
    if (!shell.services().legal) return;
    const auto r = shell.services().legal->score(m);

    ImGui::TextDisabled("JURISDICTION");
    ImGui::TextUnformatted(r.jurisdiction.c_str());
    ImGui::Spacing();
    ImGui::TextDisabled("SUBSTANTIAL SIMILARITY");
    ImGui::ProgressBar(static_cast<float>(r.substantialSimilarity / 100.0), ImVec2(-1, 0),
                       (f0(r.substantialSimilarity) + " / 100").c_str());
    ImGui::Spacing();
    ImGui::TextUnformatted("Classification:");
    ImGui::SameLine();
    const Verdict band = r.substantialSimilarity >= 75 ? Verdict::Danger
                         : (r.substantialSimilarity >= 50 ? Verdict::Warn : Verdict::Good);
    ImGui::TextColored(theme::verdictColor(static_cast<int>(band)), "%s", r.classification.c_str());
    ImGui::Separator();
    for (const auto& line : r.rationale) ImGui::BulletText("%s", line.c_str());
}

// --------------------------------------------------------------------- Docking
void docking(AppShell& shell) {
    const Molecule m = shell.currentMolecule();
    Services& s = shell.services();
    if (!s.docking) return;

    auto targets = s.docking->targets();
    UiState& st = shell.state();
    if (st.dockTarget.empty() && !targets.empty()) st.dockTarget = targets.front();

    ImGui::TextDisabled("TARGET (binding/pharmacology only)");
    ImGui::SetNextItemWidth(360);
    if (ImGui::BeginCombo("##target", st.dockTarget.c_str())) {
        for (const auto& t : targets) {
            const bool sel = (t == st.dockTarget);
            if (ImGui::Selectable(t.c_str(), sel)) st.dockTarget = t;
            if (sel) ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }
    ImGui::Spacing();

    const auto r = s.docking->dock(m, st.dockTarget);
    statCard("BEST AFFINITY", f2(r.bestAffinity), "kcal/mol (more negative = stronger)", 320.0f);
    ImGui::Spacing();

    if (!r.poses.empty() &&
        ImPlot::BeginPlot("##dock", ImVec2(-1, 190), ImPlotFlags_NoMouseText | ImPlotFlags_NoLegend)) {
        std::vector<double> xs, ys;
        for (const auto& p : r.poses) { xs.push_back(p.rank); ys.push_back(p.affinityKcalPerMol); }
        ImPlot::SetupAxes("pose rank", "affinity (kcal/mol)", 0, 0);
        ImPlot::PlotBars("affinity", ys.data(), static_cast<int>(ys.size()), 0.5, 1.0);
        ImPlot::EndPlot();
    }

    if (ImGui::BeginTable("poses", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Rank");
        ImGui::TableSetupColumn("Affinity (kcal/mol)");
        ImGui::TableSetupColumn("RMSD");
        ImGui::TableHeadersRow();
        for (const auto& p : r.poses) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0); ImGui::Text("%d", p.rank);
            ImGui::TableSetColumnIndex(1); ImGui::TextUnformatted(f2(p.affinityKcalPerMol).c_str());
            ImGui::TableSetColumnIndex(2); ImGui::TextUnformatted(f2(p.rmsd).c_str());
        }
        ImGui::EndTable();
    }
    ImGui::Spacing();
    ImGui::TextDisabled("Note: binding affinity is a target-engagement signal, never a make-it signal.");
}

// --------------------------------------------------------------------- Library
void library(AppShell& shell) {
    Services& s = shell.services();
    if (!s.library) return;
    static char filter[64] = {0};
    ImGui::SetNextItemWidth(320);
    ImGui::InputTextWithHint("##filter", "Filter by name or class...", filter, sizeof(filter));
    ImGui::SameLine();
    ImGui::TextDisabled("%zu compounds", s.library->count());

    std::string needle(filter);
    std::transform(needle.begin(), needle.end(), needle.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (ImGui::BeginTable("lib", 5,
                          ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                          ImGuiTableFlags_ScrollY | ImGuiTableFlags_Resizable,
                          ImVec2(0, 380))) {
        ImGui::TableSetupColumn("Name");
        ImGui::TableSetupColumn("Formula", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("MW", ImGuiTableColumnFlags_WidthFixed, 70.0f);
        ImGui::TableSetupColumn("Class");
        ImGui::TableSetupColumn("Legal");
        ImGui::TableHeadersRow();
        for (const auto& mol : s.library->all()) {
            std::string hay = mol.name + " " + mol.drugClass;
            std::transform(hay.begin(), hay.end(), hay.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (!needle.empty() && hay.find(needle) == std::string::npos) continue;

            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            const bool selected = (mol.id == shell.state().selectedMolecule);
            if (ImGui::Selectable(mol.name.c_str(), selected, ImGuiSelectableFlags_SpanAllColumns))
                shell.state().selectedMolecule = mol.id;
            ImGui::TableSetColumnIndex(1); ImGui::TextUnformatted(mol.formula.c_str());
            ImGui::TableSetColumnIndex(2); ImGui::TextUnformatted(f0(mol.molWeight).c_str());
            ImGui::TableSetColumnIndex(3); ImGui::TextUnformatted(mol.drugClass.c_str());
            ImGui::TableSetColumnIndex(4); ImGui::TextUnformatted(mol.legalStatus.c_str());
        }
        ImGui::EndTable();
    }
}

// ------------------------------------------------------------------------ Runs
void runs(AppShell& shell) {
    Services& s = shell.services();
    if (!s.runs) return;
    if (ImGui::BeginTable("runs", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("ID", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Kind", ImGuiTableColumnFlags_WidthFixed, 110.0f);
        ImGui::TableSetupColumn("Subject");
        ImGui::TableSetupColumn("Status", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Summary");
        ImGui::TableHeadersRow();
        for (const auto& r : s.runs->recent()) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(r.id.c_str());
            ImGui::TableSetColumnIndex(1); ImGui::TextUnformatted(r.kind.c_str());
            ImGui::TableSetColumnIndex(2); ImGui::TextUnformatted(r.subject.c_str());
            ImGui::TableSetColumnIndex(3);
            ImGui::TextColored(theme::verdictColor(r.status == "complete" ? 1 : 0), "%s",
                               r.status.c_str());
            ImGui::TableSetColumnIndex(4); ImGui::TextWrapped("%s", r.summary.c_str());
        }
        ImGui::EndTable();
    }
}

// --------------------------------------------------------------------- Presets
void presets(AppShell& shell) {
    Services& s = shell.services();
    ImGui::TextWrapped(
        "CNS target presets used for docking/pharmacology. The full set of 29 curated presets "
        "(with PDB references and binding-site boxes) loads from YAML in Phase C.");
    ImGui::Separator();
    if (s.docking) {
        for (const auto& t : s.docking->targets()) ImGui::BulletText("%s", t.c_str());
    }
}

// -------------------------------------------------------------------- Settings
void settings(AppShell& shell) {
    (void)shell;
    ImGui::TextDisabled("AI PROVIDER");
    static int provider = 0;
    const char* providers[] = {"Anthropic", "OpenAI-compatible", "DeepSeek", "Ollama (local)"};
    ImGui::SetNextItemWidth(280);
    ImGui::Combo("##prov", &provider, providers, IM_ARRAYSIZE(providers));
    static char key[96] = {0};
    ImGui::SetNextItemWidth(360);
    ImGui::InputTextWithHint("##key", "API key (stored via Windows DPAPI)...", key, sizeof(key),
                             ImGuiInputTextFlags_Password);
    ImGui::TextDisabled("Keys are encrypted at rest; the live agent loop arrives in Phase D.");
    ImGui::Spacing();

    ImGui::TextDisabled("COMPUTE");
    static int gpu = 0;
    ImGui::RadioButton("Auto", &gpu, 0); ImGui::SameLine();
    ImGui::RadioButton("GPU (CUDA)", &gpu, 1); ImGui::SameLine();
    ImGui::RadioButton("CPU", &gpu, 2);
    ImGui::Spacing();

    ImGui::TextDisabled("STORAGE");
    ImGui::TextWrapped("All state lives under %%APPDATA%%/StimLab (db, artifacts, runtime, presets, logs).");
}

// ------------------------------------------------------------------- Compare
void compare(AppShell& shell) {
    Services& s = shell.services();
    if (!s.library || !s.stability || !s.absorption || !s.admet) return;
    const auto lib = s.library->all();
    if (lib.size() < 2) return;

    static int slot[3] = {1, 2, 13};  // amphetamine vs methamphetamine vs nicotine-ish defaults
    ImGui::TextDisabled("Pick up to three compounds to compare an analog against existing samples.");
    for (int i = 0; i < 3; ++i) {
        ImGui::PushID(i);
        ImGui::SetNextItemWidth(220);
        const std::string cur =
            (slot[i] >= 0 && slot[i] < static_cast<int>(lib.size())) ? lib[slot[i]].name : "(none)";
        if (ImGui::BeginCombo("##cmp", cur.c_str())) {
            if (ImGui::Selectable("(none)", slot[i] < 0)) slot[i] = -1;
            for (int k = 0; k < static_cast<int>(lib.size()); ++k)
                if (ImGui::Selectable(lib[k].name.c_str(), slot[i] == k)) slot[i] = k;
            ImGui::EndCombo();
        }
        if (i < 2) ImGui::SameLine();
        ImGui::PopID();
    }

    struct Row {
        Molecule m;
        StabilityReport st;
        AbsorptionReport ab;
        AdmetReport ad;
    };
    std::vector<Row> rows;
    for (int i = 0; i < 3; ++i) {
        if (slot[i] < 0 || slot[i] >= static_cast<int>(lib.size())) continue;
        const Molecule& m = lib[slot[i]];
        rows.push_back({m, s.stability->analyze(m), s.absorption->predict(m), s.admet->screen(m)});
    }
    if (rows.empty()) { ImGui::TextDisabled("Select at least one compound."); return; }

    if (ImPlot::BeginPlot("##cmpchart", ImVec2(-1, 230))) {
        constexpr int G = 4;
        const int I = static_cast<int>(rows.size());
        const char* groups[G] = {"Stability", "Oral F%", "HIA%", "Oxidation"};
        const double gpos[G] = {0, 1, 2, 3};
        std::vector<double> data(static_cast<size_t>(I) * G);
        std::vector<const char*> items;
        for (int i = 0; i < I; ++i) {
            items.push_back(rows[i].m.name.c_str());
            double ox = 0;
            for (const auto& f : rows[i].st.factors)
                if (f.name.find("Oxidation") != std::string::npos) ox = f.score;
            data[static_cast<size_t>(i) * G + 0] = rows[i].st.overallScore;
            data[static_cast<size_t>(i) * G + 1] = rows[i].ab.bioavailabilityPct;
            data[static_cast<size_t>(i) * G + 2] = rows[i].ab.hiaPct;
            data[static_cast<size_t>(i) * G + 3] = ox;
        }
        ImPlot::SetupAxes(nullptr, "score / %");
        ImPlot::SetupAxesLimits(-0.5, 3.5, 0, 100, ImPlotCond_Always);
        ImPlot::SetupAxisTicks(ImAxis_X1, gpos, G, groups);
        ImPlot::PlotBarGroups(items.data(), data.data(), I, G, 0.6);
        ImPlot::EndPlot();
    }

    if (ImGui::BeginTable("cmptab", 1 + static_cast<int>(rows.size()),
                          ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Metric", ImGuiTableColumnFlags_WidthFixed, 130.0f);
        for (const auto& r : rows) ImGui::TableSetupColumn(r.m.name.c_str());
        ImGui::TableHeadersRow();
        auto rowD = [&](const char* label, const std::function<std::string(const Row&)>& fn) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0); ImGui::TextDisabled("%s", label);
            for (int i = 0; i < static_cast<int>(rows.size()); ++i) {
                ImGui::TableSetColumnIndex(i + 1);
                ImGui::TextUnformatted(fn(rows[i]).c_str());
            }
        };
        rowD("Formula", [](const Row& r) { return r.m.formula; });
        rowD("MW", [](const Row& r) { return f0(r.m.molWeight); });
        rowD("logP", [](const Row& r) { return f2(r.m.logP); });
        rowD("TPSA", [](const Row& r) { return f2(r.m.tpsa); });
        rowD("Stability", [](const Row& r) { return f0(r.st.overallScore) + "/100"; });
        rowD("Shelf life", [](const Row& r) { return r.st.shelfLifeEstimate; });
        rowD("Oral F%", [](const Row& r) { return f0(r.ab.bioavailabilityPct) + "%"; });
        rowD("HIA%", [](const Row& r) { return f0(r.ab.hiaPct) + "%"; });
        rowD("logBB", [](const Row& r) { return f2(r.ab.logBB) + (r.ab.cnsPenetrant ? " (CNS)" : ""); });
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0); ImGui::TextDisabled("ADMET");
        for (int i = 0; i < static_cast<int>(rows.size()); ++i) {
            ImGui::TableSetColumnIndex(i + 1);
            verdictText(rows[i].ad.overall);
        }
        rowD("Legal status", [](const Row& r) { return r.m.legalStatus; });
        ImGui::EndTable();
    }
}

// ------------------------------------------------------------ Analog Explorer
void analogExplorer(AppShell& shell) {
    Services& s = shell.services();
    if (!s.library || !s.stability || !s.absorption || !s.admet || !s.similarity || !s.legal) return;
    const auto lib = s.library->all();

    static int parentIdx = 1;  // amphetamine
    static float logP = 1.76f, tpsa = 26.0f, mw = 135.0f;
    static int hbd = 1, hba = 1, rot = 2;
    static bool ester = false, catechol = false, arylKetone = false, mdoxy = false;
    static int lastParent = -999;

    ImGui::TextWrapped("Model a candidate analog's property profile, then check it against existing "
                       "samples and screen predicted byproducts/interactions. (Analysis only - no "
                       "synthesis guidance.)");
    ImGui::SetNextItemWidth(260);
    const std::string pcur =
        (parentIdx >= 0 && parentIdx < static_cast<int>(lib.size())) ? lib[parentIdx].name : "(scratch)";
    if (ImGui::BeginCombo("Seed from parent", pcur.c_str())) {
        if (ImGui::Selectable("(scratch)", parentIdx < 0)) parentIdx = -1;
        for (int k = 0; k < static_cast<int>(lib.size()); ++k)
            if (ImGui::Selectable(lib[k].name.c_str(), parentIdx == k)) parentIdx = k;
        ImGui::EndCombo();
    }
    if (parentIdx != lastParent) {
        lastParent = parentIdx;
        if (parentIdx >= 0 && parentIdx < static_cast<int>(lib.size())) {
            const Molecule& p = lib[parentIdx];
            logP = static_cast<float>(p.logP); tpsa = static_cast<float>(p.tpsa);
            mw = static_cast<float>(p.molWeight);
            hbd = p.hbd; hba = p.hba; rot = p.rotatableBonds;
            ester = p.smiles.find("OC(=O)") != std::string::npos ||
                    p.smiles.find("C(=O)OC") != std::string::npos;
            catechol = p.smiles.find("c(O)c(O)") != std::string::npos;
            arylKetone = p.smiles.find("C(=O)c") != std::string::npos;
            mdoxy = p.smiles.find("OCOc") != std::string::npos;
        }
    }

    ImGui::Separator();
    if (ImGui::BeginTable("aeknobs", 2, ImGuiTableFlags_BordersInnerV)) {
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::TextDisabled("PROPERTIES");
        ImGui::SetNextItemWidth(220);
        ImGui::SliderFloat("logP", &logP, -3.0f, 5.0f, "%.2f");
        ImGui::SetNextItemWidth(220);
        ImGui::SliderFloat("TPSA", &tpsa, 0.0f, 150.0f, "%.0f");
        ImGui::SetNextItemWidth(220);
        ImGui::SliderFloat("MW", &mw, 80.0f, 400.0f, "%.0f");
        ImGui::SetNextItemWidth(220);
        ImGui::SliderInt("H-bond donors", &hbd, 0, 6);
        ImGui::SetNextItemWidth(220);
        ImGui::SliderInt("H-bond acceptors", &hba, 0, 10);
        ImGui::SetNextItemWidth(220);
        ImGui::SliderInt("Rotatable bonds", &rot, 0, 12);
        ImGui::TableSetColumnIndex(1);
        ImGui::TextDisabled("FUNCTIONAL GROUPS");
        ImGui::Checkbox("Ester (hydrolysis-labile)", &ester);
        ImGui::Checkbox("Catechol (oxidation / COMT)", &catechol);
        ImGui::Checkbox("Aryl ketone (beta-keto)", &arylKetone);
        ImGui::Checkbox("Methylenedioxy ring", &mdoxy);
        ImGui::EndTable();
    }

    Molecule c;
    c.id = "candidate";
    c.name = "Candidate analog";
    c.logP = logP; c.tpsa = tpsa; c.molWeight = mw;
    c.hbd = hbd; c.hba = hba; c.rotatableBonds = rot;
    c.formula = "(modeled)";
    c.drugClass = (parentIdx >= 0 && parentIdx < static_cast<int>(lib.size()))
                      ? lib[parentIdx].drugClass : "Phenethylamine (candidate)";
    const std::string ring = catechol ? "c1ccc(O)c(O)c1" : (mdoxy ? "c1ccc2OCOc2c1" : "c1ccccc1");
    c.smiles = ring + "CC(N)C" + (ester ? "OC(=O)C" : "") + (arylKetone ? "C(=O)c1ccccc1" : "");

    const auto st = s.stability->analyze(c);
    const auto ab = s.absorption->predict(c);
    const auto ad = s.admet->screen(c);
    const auto si = s.similarity->search(c);
    const auto lg = s.legal->score(c);

    ImGui::Separator();
    statCard("STABILITY", f0(st.overallScore) + "/100", st.shelfLifeEstimate.c_str(), 190.0f);
    ImGui::SameLine();
    statCard("ORAL F", f0(ab.bioavailabilityPct) + "%",
             ab.cnsPenetrant ? "CNS-penetrant" : "peripheral", 160.0f);
    ImGui::SameLine();
    statCard("HIA", f0(ab.hiaPct) + "%", "intestinal abs.", 150.0f);
    ImGui::SameLine();
    {
        ImGui::BeginChild("aeADMET", ImVec2(170, 84), ImGuiChildFlags_Borders);
        ImGui::TextDisabled("ADMET");
        ImGui::SetWindowFontScale(1.3f);
        verdictText(ad.overall);
        ImGui::SetWindowFontScale(1.0f);
        ImGui::TextDisabled("%zu endpoints", ad.endpoints.size());
        ImGui::EndChild();
    }
    ImGui::Spacing();

    ImGui::TextDisabled("NEAREST EXISTING SAMPLE");
    if (!si.hits.empty()) {
        ImGui::Text("%s", si.nearestName.c_str());
        ImGui::SameLine();
        ImGui::TextDisabled("(Tanimoto %.2f - %s, %s)", si.nearestScore,
                            si.hits.front().referenceClass.c_str(),
                            si.hits.front().legalStatus.c_str());
    }
    ImGui::TextDisabled("LEGAL-ANALOG SCORE");
    const Verdict lband = lg.substantialSimilarity >= 75 ? Verdict::Danger
                          : (lg.substantialSimilarity >= 50 ? Verdict::Warn : Verdict::Good);
    ImGui::PushStyleColor(ImGuiCol_PlotHistogram, theme::verdictColor(static_cast<int>(lband)));
    ImGui::ProgressBar(static_cast<float>(lg.substantialSimilarity / 100.0), ImVec2(330, 0),
                       (f0(lg.substantialSimilarity) + " / 100  -  " + lg.classification).c_str());
    ImGui::PopStyleColor();

    ImGui::Spacing();
    ImGui::TextDisabled("PREDICTED BYPRODUCTS / INTERACTIONS");
    for (const auto& d : st.degradants)
        ImGui::BulletText("Degradant: %s  (%s)", d.name.c_str(), d.pathway.c_str());
    for (const auto& e : ad.endpoints) {
        if (e.verdict == Verdict::Warn || e.verdict == Verdict::Danger) {
            ImGui::Bullet();
            verdictText(e.verdict);
            ImGui::SameLine();
            ImGui::TextWrapped("%s - %s", e.name.c_str(), e.detail.c_str());
        }
    }
}

}  // namespace stimlab::panels
