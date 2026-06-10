#include "ui/Panels.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <functional>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include <imgui.h>
#include <implot.h>

#include "chem/Descriptors.h"
#include "chem/Embed3D.h"
#include "chem/Smiles.h"
#include "core/AppPaths.h"
#include "core/Manifest.h"
#include "modules/docking/Presets.h"
#include "modules/docking/Provisioning.h"
#include "modules/docking/ReceptorPrep.h"
#include "render/MolViewport.h"
#include "ui/AppShell.h"
#include "ui/Theme.h"
#include "workflow/Dag.h"

namespace stimlab::panels {
namespace {

std::string f2(double v) { char b[40]; std::snprintf(b, sizeof b, "%.2f", v); return b; }
std::string f0(double v) { char b[40]; std::snprintf(b, sizeof b, "%.0f", v); return b; }

void verdictText(Verdict v) {
    ImGui::TextColored(theme::verdictColor(static_cast<int>(v)), "%s", verdictLabel(v));
}

void statCard(const char* title, const std::string& value, const char* sub, float w = 168.0f,
              float h = 84.0f, float valueScale = 1.6f) {
    ImGui::BeginChild(title, ImVec2(w, h), ImGuiChildFlags_Borders);
    ImGui::TextDisabled("%s", title);
    ImGui::SetWindowFontScale(valueScale);
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

// Persistent per-panel state for an embedded 3D molecular viewport.
struct ViewerUiState {
    bool spacefill = false;
    bool showH = false;
    bool showReceptor = true;  // overlay the binding-pocket receptor (docking only)
    std::string lastKey;       // cache key of the molecule currently framed
};

// CPK legend swatch (color comes from the same table the renderer uses).
void cpkSwatch(int z, const char* label) {
    const std::uint32_t c = render::cpkColor(z);  // 0xAABBGGRR
    const ImU32 col = IM_COL32(c & 0xFF, (c >> 8) & 0xFF, (c >> 16) & 0xFF, 0xFF);
    ImGui::PushStyleColor(ImGuiCol_Button, col);
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, col);
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, col);
    ImGui::Button("##sw", ImVec2(14, 14));
    ImGui::PopStyleColor(3);
    ImGui::SameLine(0, 4);
    ImGui::TextUnformatted(label);
}

// Parse a prepared receptor PDBQT into the BINDING-POCKET subset: heavy receptor
// atoms within `radius` A of any ligand heavy atom (the pose and the receptor share
// one coordinate frame, so they overlay directly). Bonds are inferred by covalent
// distance so the pocket draws as a wireframe context around the docked ligand.
chem::Conformer loadReceptorPocket(const std::string& path, const chem::Conformer& ligand,
                                   double radius) {
    chem::Conformer out;
    if (path.empty() || ligand.empty()) return out;
    std::ifstream in(path, std::ios::binary);
    if (!in) return out;

    std::vector<chem::Vec3> lig;
    for (int i = 0; i < ligand.size(); ++i)
        if (ligand.z[i] != 1) lig.push_back(ligand.pos[i]);
    if (lig.empty())
        for (int i = 0; i < ligand.size(); ++i) lig.push_back(ligand.pos[i]);
    const double r2 = radius * radius;

    auto adTypeToZ = [](std::string t) -> int {
        for (auto& ch : t) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        if (t == "C" || t == "A") return 6;
        if (t == "N" || t == "NA" || t == "NS") return 7;
        if (t == "O" || t == "OA" || t == "OS") return 8;
        if (t == "S" || t == "SA") return 16;
        if (t == "H" || t == "HD" || t == "HS") return 1;
        if (t == "P") return 15;   if (t == "F") return 9;
        if (t == "CL") return 17;  if (t == "BR") return 35;  if (t == "I") return 53;
        if (t == "FE") return 26;  if (t == "ZN") return 30;  if (t == "MG") return 12;
        if (t == "MN") return 25;  if (t == "CA") return 20;  if (t == "K") return 19;
        return 6;
    };
    auto lastToken = [](const std::string& s) -> std::string {
        const size_t e = s.find_last_not_of(" \t\r\n");
        if (e == std::string::npos) return {};
        const size_t b = s.find_last_of(" \t", e);
        return s.substr(b == std::string::npos ? 0 : b + 1, e - (b == std::string::npos ? 0 : b));
    };

    std::string line;
    while (std::getline(in, line)) {
        if (!(line.rfind("ATOM", 0) == 0 || line.rfind("HETATM", 0) == 0)) continue;
        if (line.size() < 54) continue;
        double x, y, z;
        try {
            x = std::stod(line.substr(30, 8));
            y = std::stod(line.substr(38, 8));
            z = std::stod(line.substr(46, 8));
        } catch (...) { continue; }
        if (!(std::isfinite(x) && std::isfinite(y) && std::isfinite(z))) continue;
        const int zz = adTypeToZ(lastToken(line));
        if (zz == 1) continue;  // skip receptor hydrogens (clutter)
        bool near = false;
        for (const auto& L : lig) {
            const double dx = L.x - x, dy = L.y - y, dz = L.z - z;
            if (dx * dx + dy * dy + dz * dz <= r2) { near = true; break; }
        }
        if (!near) continue;
        out.pos.push_back({x, y, z});
        out.z.push_back(zz);
        if (static_cast<int>(out.pos.size()) >= 1200) break;  // cap for performance
    }
    out.heavyCount = static_cast<int>(out.pos.size());
    const int N = static_cast<int>(out.pos.size());
    for (int i = 0; i < N; ++i)
        for (int j = i + 1; j < N; ++j) {
            const double dx = out.pos[i].x - out.pos[j].x, dy = out.pos[i].y - out.pos[j].y,
                         dz = out.pos[i].z - out.pos[j].z;
            const double d2 = dx * dx + dy * dy + dz * dz;
            const double cut = chem::covalentRadius(out.z[i]) + chem::covalentRadius(out.z[j]) + 0.45;
            if (d2 > 0.16 && d2 <= cut * cut) out.bonds.push_back({i, j});
        }
    return out;
}

// Overlay the receptor pocket into an existing ligand scene as muted thin sticks +
// small spheres, and grow the bounding radius so the pocket frames around the ligand.
void appendReceptorToScene(render::MolScene& scene, const chem::Conformer& rec) {
    if (rec.empty()) return;
    auto dim = [](std::uint32_t c) -> std::uint32_t {  // c = 0xAABBGGRR
        unsigned r = c & 0xFF, g = (c >> 8) & 0xFF, b = (c >> 16) & 0xFF;
        auto mix = [](unsigned x, unsigned t) { return static_cast<unsigned>(x * 0.42 + t * 0.58); };
        r = mix(r, 0x5a); g = mix(g, 0x64); b = mix(b, 0x72);
        return 0xFF000000u | (b << 16) | (g << 8) | r;
    };
    for (int i = 0; i < rec.size(); ++i) {
        render::AtomInst a;
        a.x = static_cast<float>(rec.pos[i].x);
        a.y = static_cast<float>(rec.pos[i].y);
        a.z = static_cast<float>(rec.pos[i].z);
        a.r = 0.26f;
        a.rgba = dim(render::cpkColor(rec.z[i]));
        scene.atoms.push_back(a);
    }
    for (const auto& bp : rec.bonds) {
        const int ia = bp.first, ib = bp.second;
        if (ia < 0 || ib < 0 || ia >= rec.size() || ib >= rec.size()) continue;
        render::BondInst b;
        b.ax = static_cast<float>(rec.pos[ia].x); b.ay = static_cast<float>(rec.pos[ia].y);
        b.az = static_cast<float>(rec.pos[ia].z);
        b.bx = static_cast<float>(rec.pos[ib].x); b.by = static_cast<float>(rec.pos[ib].y);
        b.bz = static_cast<float>(rec.pos[ib].z);
        b.rgbaA = dim(render::cpkColor(rec.z[ia]));
        b.rgbaB = dim(render::cpkColor(rec.z[ib]));
        b.radius = 0.07f;
        scene.bonds.push_back(b);
    }
    double maxR = scene.radius;
    for (int i = 0; i < rec.size(); ++i) {
        const double dx = rec.pos[i].x - scene.center.x, dy = rec.pos[i].y - scene.center.y,
                     dz = rec.pos[i].z - scene.center.z;
        maxR = std::max(maxR, std::sqrt(dx * dx + dy * dy + dz * dz) + 0.3);
    }
    scene.radius = static_cast<float>(maxR);
}

// Draw a chem::Conformer in an interactive 3D viewport (off-screen RT -> ImGui
// image) with toggles + a CPK legend. `cacheKey` identifies the molecule so the
// camera only auto-fits when the displayed structure changes. Returns false (and
// draws a fallback note) if no GPU viewer is available. When `receptor` is given
// (docking), its binding pocket is overlaid behind the ligand.
bool molViewer3D(AppShell& shell, const chem::Conformer& conf, const std::string& cacheKey,
                 ViewerUiState& ui, float height, const chem::Conformer* receptor = nullptr) {
    render::MolViewport* vp = shell.viewer();
    if (!vp) {
        ImGui::TextDisabled("3D viewer unavailable (no DirectX device).");
        return false;
    }

    // Controls row.
    ImGui::Checkbox("Spacefill", &ui.spacefill);
    ImGui::SameLine();
    ImGui::Checkbox("Show H", &ui.showH);
    if (receptor && !receptor->empty()) {
        ImGui::SameLine();
        ImGui::Checkbox("Receptor", &ui.showReceptor);
    }
    ImGui::SameLine();
    if (ImGui::SmallButton("Reset view")) ui.lastKey.clear();  // force a re-fit
    ImGui::SameLine();
    // Explicit zoom controls so zoom is discoverable without the wheel.
    if (ImGui::SmallButton("+##zoomin")) vp->zoom(1.0f);
    ImGui::SameLine();
    if (ImGui::SmallButton("-##zoomout")) vp->zoom(-1.0f);
    ImGui::SameLine();
    ImGui::TextDisabled("(drag: orbit  -  wheel/+/-: zoom  -  right-drag: pan)");

    const bool overlayReceptor = receptor && !receptor->empty() && ui.showReceptor;
    render::MolScene scene = render::buildMolScene(conf, ui.spacefill, ui.showH);
    if (overlayReceptor) appendReceptorToScene(scene, *receptor);

    // Auto-fit when the molecule (or toggle that changes framing) changes.
    const std::string key = cacheKey + (ui.showH ? "|H" : "") + (ui.spacefill ? "|S" : "") +
                            (overlayReceptor ? "|R" : "");
    if (key != ui.lastKey) {
        vp->resetCamera(scene);
        ui.lastKey = key;
    }

    const float avail = ImGui::GetContentRegionAvail().x;
    const int w = std::max(64, static_cast<int>(avail));
    const int h = std::max(64, static_cast<int>(height));

    if (conf.empty() || scene.atoms.empty()) {
        ImGui::Dummy(ImVec2(static_cast<float>(w), static_cast<float>(h)));
        ImGui::TextDisabled("(no atoms to display)");
        return true;
    }

    void* srv = vp->draw(scene, w, h);
    if (!srv) {
        ImGui::TextDisabled("3D render failed.");
        return false;
    }

    // Host the rendered image in its own child window flagged NoScrollWithMouse +
    // NoScrollbar. With both flags set, ImGui absorbs the mouse wheel over this child
    // instead of forwarding it to the scrolling panel, so hovering + wheel zooms the
    // model directly (no need to scroll the page up first). Zero window padding so the
    // image fills the child exactly (no offset / no spurious scroll region).
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
    ImGui::BeginChild("##vp3d", ImVec2(static_cast<float>(w), static_cast<float>(h)),
                      ImGuiChildFlags_None,
                      ImGuiWindowFlags_NoScrollWithMouse | ImGuiWindowFlags_NoScrollbar);
    ImGui::PopStyleVar();
    ImGui::Image(reinterpret_cast<ImTextureID>(srv), ImVec2(static_cast<float>(w), static_cast<float>(h)));
    const bool hovered = ImGui::IsItemHovered();

    // Camera input - only while the image is hovered.
    ImGuiIO& io = ImGui::GetIO();
    if (hovered) {
        if (ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
            const ImVec2 d = io.MouseDelta;
            vp->orbit(d.x, d.y);
        }
        if (ImGui::IsMouseDragging(ImGuiMouseButton_Right) ||
            ImGui::IsMouseDragging(ImGuiMouseButton_Middle)) {
            const ImVec2 d = io.MouseDelta;
            vp->pan(d.x, d.y);
        }
        if (io.MouseWheel != 0.0f) vp->zoom(io.MouseWheel);
    }
    ImGui::EndChild();

    // Legend (only elements that actually appear).
    bool seen[120] = {false};
    for (int i = 0; i < conf.size(); ++i) {
        const int z = conf.z[i];
        if (z >= 0 && z < 120) seen[z] = true;
    }
    struct LE { int z; const char* sym; };
    const LE order[] = {{6, "C"}, {1, "H"}, {7, "N"}, {8, "O"}, {16, "S"}, {15, "P"},
                        {9, "F"}, {17, "Cl"}, {35, "Br"}, {53, "I"}};
    bool first = true;
    for (const auto& e : order) {
        if (!seen[e.z]) continue;
        if (e.z == 1 && !ui.showH) continue;
        if (!first) ImGui::SameLine(0, 12);
        first = false;
        cpkSwatch(e.z, e.sym);
    }
    return true;
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

    // Embed the selected molecule to a 3D conformer, caching by id so we only
    // re-embed when the active compound changes (embedding is cheap but not free).
    static std::string lastId;
    static chem::Conformer conf;
    static ViewerUiState viewUi;
    if (m.id != lastId) {
        lastId = m.id;
        conf = chem::Conformer{};
        if (auto parsed = chem::parseSmiles(m.smiles)) conf = chem::embed3D(*parsed);
    }

    if (ImGui::BeginTable("idprop", 2, ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_Resizable)) {
        ImGui::TableSetupColumn("Identity", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("3D structure", ImGuiTableColumnFlags_WidthFixed, 360.0f);
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
        ImGui::TextDisabled("3D STRUCTURE");
        if (!molViewer3D(shell, conf, m.id, viewUi, 300.0f)) {
            // Graceful fallback when no GPU device is available: 2D schematic.
            moleculeSchematic(m, ImVec2(260, 220));
        }
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

// --------------------------------------------------------------- Molecule Input
// Free-text SMILES -> the in-house parser -> full descriptors + 3D embedding, so
// arbitrary structures (not just library compounds) can be analyzed. "Set as
// active compound" routes the parsed molecule into every other panel.
void moleculeInput(AppShell& shell) {
    static char smilesBuf[256] = "CC(N)Cc1ccc(Cl)cc1";  // 4-chloroamphetamine example
    static char nameBuf[96] = "My candidate";
    static std::string status;
    static chem::Conformer previewConf;
    static ViewerUiState viewUi;
    static Molecule preview;
    static bool valid = false;
    static std::string lastSmiles;

    ImGui::TextWrapped(
        "Enter any SMILES to analyze an arbitrary structure with the full StimLab engine "
        "(descriptors, 3D embedding, ADMET, similarity, docking). Analysis only - no synthesis "
        "guidance.");
    ImGui::Spacing();
    ImGui::SetNextItemWidth(440);
    ImGui::InputText("SMILES", smilesBuf, sizeof(smilesBuf));
    ImGui::SetNextItemWidth(260);
    ImGui::InputText("Name", nameBuf, sizeof(nameBuf));
    ImGui::SameLine();
    const bool analyze = ImGui::Button("Analyze");

    if (analyze || lastSmiles != smilesBuf) {
        lastSmiles = smilesBuf;
        auto parsed = chem::parseSmiles(smilesBuf);
        if (parsed && !parsed->empty()) {
            valid = true;
            preview = Molecule{};
            preview.id = "__custom__";
            preview.name = nameBuf[0] ? std::string(nameBuf) : "(custom)";
            preview.smiles = smilesBuf;
            preview.formula = chem::molecularFormula(*parsed);
            preview.molWeight = chem::molecularWeight(*parsed);
            preview.logP = chem::crippenLogP(*parsed);
            preview.tpsa = chem::tpsa(*parsed);
            preview.hbd = chem::hbdCount(*parsed);
            preview.hba = chem::hbaCount(*parsed);
            preview.rotatableBonds = chem::rotatableBondCount(*parsed);
            preview.drugClass = "User input";
            preview.legalStatus = "(unscheduled / unknown)";
            preview.notes = "Entered via the Molecule Input panel.";
            previewConf = chem::embed3D(*parsed);
            status = "Parsed OK - " + preview.formula;
        } else {
            valid = false;
            status = "Invalid SMILES - check brackets, ring closures and atom symbols.";
        }
    }

    ImGui::Spacing();
    ImGui::TextColored(theme::verdictColor(valid ? 1 : 2), "%s", status.c_str());
    ImGui::Separator();
    if (!valid) return;

    if (ImGui::BeginTable("inprop", 2, ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_Resizable)) {
        ImGui::TableSetupColumn("Properties", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("3D structure", ImGuiTableColumnFlags_WidthFixed, 360.0f);
        ImGui::TableNextRow();

        ImGui::TableSetColumnIndex(0);
        ImGui::TextDisabled("Formula: %s", preview.formula.c_str());
        ImGui::Spacing();
        if (ImGui::BeginTable("inp", 6, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
            for (const char* h : {"MW", "logP", "TPSA", "HBD", "HBA", "RotB"})
                ImGui::TableSetupColumn(h);
            ImGui::TableHeadersRow();
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(f2(preview.molWeight).c_str());
            ImGui::TableSetColumnIndex(1); ImGui::TextUnformatted(f2(preview.logP).c_str());
            ImGui::TableSetColumnIndex(2); ImGui::TextUnformatted(f2(preview.tpsa).c_str());
            ImGui::TableSetColumnIndex(3); ImGui::Text("%d", preview.hbd);
            ImGui::TableSetColumnIndex(4); ImGui::Text("%d", preview.hba);
            ImGui::TableSetColumnIndex(5); ImGui::Text("%d", preview.rotatableBonds);
            ImGui::EndTable();
        }
        ImGui::Spacing();
        if (ImGui::Button("Set as active compound")) {
            shell.state().customMolecule = preview;
            shell.state().hasCustom = true;
            shell.state().selectedMolecule = "__custom__";
            shell.requestHighlight("Structure",
                "Your structure is now the active compound - every panel (Structure, ADMET, "
                "Docking, ...) analyzes it.");
        }
        ImGui::SameLine();
        ImGui::TextDisabled("routes this structure into every other panel");

        ImGui::TableSetColumnIndex(1);
        ImGui::TextDisabled("3D STRUCTURE");
        molViewer3D(shell, previewConf, "custom|" + lastSmiles, viewUi, 300.0f);
        ImGui::EndTable();
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

    // ---- Real-engine provisioning (WP-F/WP-G): vina.exe + prepared receptors -----
    {
        auto& prov = shell.provisioner();
        // When a user-triggered (download) provision finishes, drop any dock computed
        // while the engine was still absent so the next frame re-docks for real. The
        // cheap startup locate-only probe doesn't invalidate (the dock path already
        // reads prepared receptors from disk), avoiding a wasteful first-entry restart.
        static bool wasProvisioning = false;
        static bool userProvision = false;
        static std::string readyKey;   // target whose receptor-readiness is cached below
        static bool selReady = false;
        if (wasProvisioning && !prov.running()) {
            if (userProvision) shell.invalidateDock();
            userProvision = false;
            readyKey.clear();          // a provision just finished -> re-probe readiness
        }
        wasProvisioning = prov.running();
        if (st.dockTarget != readyKey) {  // selection changed (or invalidated above)
            readyKey = st.dockTarget;
            selReady = shell.receptorReady(st.dockTarget);
        }
        const bool vina = prov.vinaReady();
        const int recReady = prov.receptorsReady();
        const int recTotal = prov.receptorsTotal();
        const bool realCapable = vina && recReady > 0;
        ImGui::PushStyleColor(ImGuiCol_Text, theme::verdictColor(realCapable ? 1 : 2));
        ImGui::TextUnformatted(realCapable
            ? "Real docking engine: READY (vina + prepared receptor)"
            : "Real docking engine: not provisioned - showing labeled descriptor estimate");
        ImGui::PopStyleColor();
        ImGui::TextDisabled("vina.exe: %s   receptors (last provision): %d/%d   receptor prep: %s",
                            vina ? "present" : "absent", recReady, recTotal,
                            prov.obabelReady() ? "obabel (+H)" : "built-in");
        ImGui::TextDisabled("selected target %s: receptor %s", st.dockTarget.c_str(),
                            selReady ? "prepared (real dock)"
                                     : "not prepared (shows descriptor estimate)");
        if (prov.running()) {
            ImGui::TextDisabled("Working: %s", prov.status().c_str());
        } else {
            if (ImGui::Button("Provision engine + headline receptors")) {
                shell.provisionDocking();
                userProvision = true;
            }
            if (!selReady) {
                // The selected target isn't prepared yet: offer an on-demand provision
                // of just this receptor (any of the 29 CNS presets, not only headlines).
                ImGui::SameLine();
                const std::string lbl = "Provision " + st.dockTarget;
                if (ImGui::Button(lbl.c_str()) && shell.provisionTarget(st.dockTarget))
                    userProvision = true;
            }
            ImGui::SameLine();
            ImGui::TextDisabled("%s", prov.status().c_str());
        }
        ImGui::TextDisabled("Downloads vina.exe (size-checked) + prepares receptor PDBQTs from RCSB "
                            "under %%APPDATA%%/StimLab/runtime. Headline = DAT/NET/SERT/TAAR1; any "
                            "other target prepares on demand. Best-effort; needs network.");
    }
    ImGui::Spacing();

    // One async dock per (molecule, target): a real engine run happens on a worker
    // thread (never the UI thread), and both the table/plot and the 3D pose overlay
    // read this single cached result so the two views always agree.
    bool computing = false;
    const DockJobResult d = shell.dockFor(m, st.dockTarget, computing);
    const bool hasResult = !d.poses.empty();

    if (computing && !hasResult) {
        ImGui::TextDisabled("Docking %s into the %s box on a worker thread...", m.name.c_str(),
                            st.dockTarget.c_str());
        ImGui::Spacing();
    }

    // Selected pose index, shared by the chart highlight, the table, the slider and the
    // 3D overlay so every view always agrees. Declared here so the table/chart can drive
    // it; clamped to a valid range whenever the result set changes.
    static int poseSel = 0;
    if (hasResult) { if (poseSel >= static_cast<int>(d.poses.size())) poseSel = 0; }
    else poseSel = 0;

    const std::string summary =
        d.real ? ("Best affinity " + f2(d.bestAffinity()) + " kcal/mol at " + st.dockTarget +
                  " (docked with " + d.engine + ").")
               : ("Estimated affinity " + f2(d.bestAffinity()) + " kcal/mol at " + st.dockTarget +
                  " (" + d.engine + " - structure-descriptor model, not a docked score).");

    theme::sectionHeader("DOCKING RESULT");

    // Headline metric cards. The affinity card is deliberately large + high-contrast so
    // the primary number reads at a glance; companions summarise pose count / provenance.
    statCard("BEST AFFINITY", hasResult ? f2(d.bestAffinity()) : "--",
             "kcal/mol  (more negative = stronger)", 300.0f, 120.0f, 2.7f);
    ImGui::SameLine();
    statCard("POSES", hasResult ? std::to_string(d.poses.size()) : "--",
             "ranked binding modes", 150.0f, 120.0f, 2.0f);
    ImGui::SameLine();
    statCard("ENGINE", hasResult ? (d.real ? "REAL" : "ESTIMATE") : "--",
             d.real ? "docked score" : "descriptor model", 160.0f, 120.0f, 1.6f);
    if (hasResult && d.real) {
        ImGui::SameLine();
        statCard("CONFIDENCE", d.converged ? "HIGH" : "MODERATE",
                 d.converged ? "search converged" : "budget reached", 175.0f, 120.0f, 1.6f);
    }

    ImGui::Spacing();
    static std::string runSaved;
    if (ImGui::Button("Save to Runs") && hasResult) {
        RunRecord rec;
        rec.kind = "Docking";
        rec.subject = m.name + " -> " + st.dockTarget;
        rec.status = "complete";
        rec.summary = summary;
        s.runs->record(rec);
        runSaved = "Saved to run history.";
    }
    if (!runSaved.empty()) { ImGui::SameLine(); ImGui::TextDisabled("%s", runSaved.c_str()); }
    ImGui::Spacing();

    // ---- Run info: how this number was produced, consolidated in one place. ----------
    if (hasResult && ImGui::CollapsingHeader("Run info", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ImGui::BeginTable("runinfo", 2, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                                                ImGuiTableFlags_SizingStretchProp)) {
            ImGui::TableSetupColumn("Field", ImGuiTableColumnFlags_WidthFixed, 200.0f);
            ImGui::TableSetupColumn("Value");
            auto row = [](const char* k, const std::string& v) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::TextDisabled("%s", k);
                ImGui::TableSetColumnIndex(1); ImGui::TextUnformatted(v.c_str());
            };
            row("Engine", d.engine);
            row("Result type", d.real ? "real engine dock" : "descriptor estimate (not a docked score)");
            row("Target", st.dockTarget);
            if (const ReceptorTarget* preset = docking::findPreset(st.dockTarget)) {
                row("Receptor PDB", preset->pdb.empty() ? "(box only)" : preset->pdb);
                char box[96];
                std::snprintf(box, sizeof box, "%.0f x %.0f x %.0f A  @ (%.1f, %.1f, %.1f)",
                              preset->box.sx, preset->box.sy, preset->box.sz,
                              preset->box.cx, preset->box.cy, preset->box.cz);
                row("Search box", box);
            }
            if (d.real) {
                row("Independent runs", std::to_string(d.searchRuns));
                row("Convergence", d.converged ? "converged within tolerance"
                                               : "stopped at budget / run limit");
                row("Best-affinity spread", f2(d.affinitySpread) + " kcal/mol");
                row("Flexible torsions", std::to_string(d.torsions));
            }
            if (poseSel < static_cast<int>(d.poses.size())) {
                const auto& p = d.poses[poseSel];
                row("Selected pose", "rank " + std::to_string(p.rank) + " of " +
                                         std::to_string(d.poses.size()));
                row("  affinity", f2(p.affinityKcalPerMol) + " kcal/mol");
                row("  RMSD l.b. / u.b.", f2(p.rmsdLb) + " / " + f2(p.rmsdUb) + " A");
            }
            ImGui::EndTable();
        }
    }
    ImGui::Spacing();

    // ---- Affinity-by-rank chart. The axes are pinned so the origin sits at 0: pose
    // ranks start at 1 (never negative) and affinity is drawn downward from 0. The
    // selected pose is overlaid highlighted so the chart and table stay in sync. ------
    if (hasResult) {
        theme::sectionHeader("AFFINITY BY POSE");
        if (ImPlot::BeginPlot("##dock", ImVec2(-1, 200),
                              ImPlotFlags_NoMouseText | ImPlotFlags_NoLegend)) {
            const int n = static_cast<int>(d.poses.size());
            std::vector<double> ys;
            ys.reserve(n);
            double mn = 0.0;
            for (const auto& p : d.poses) {
                ys.push_back(p.affinityKcalPerMol);
                mn = std::min(mn, p.affinityKcalPerMol);
            }
            ImPlot::SetupAxes("pose rank", "affinity (kcal/mol)",
                              ImPlotAxisFlags_Lock, ImPlotAxisFlags_Lock);
            ImPlot::SetupAxisLimits(ImAxis_X1, 0.0, n + 1.0, ImPlotCond_Always);
            ImPlot::SetupAxisLimits(ImAxis_Y1, mn * 1.12 - 0.5, 0.0, ImPlotCond_Always);
            ImPlotSpec barSpec;
            barSpec.FillColor = ImVec4(0.40f, 0.55f, 0.78f, 0.85f);  // muted base bars
            ImPlot::PlotBars("affinity", ys.data(), n, 0.6, 1.0, barSpec);
            if (poseSel < n) {  // highlight the active pose on top
                const double sx = poseSel + 1.0, sy = ys[poseSel];
                ImPlotSpec selSpec;
                selSpec.FillColor = theme::verdictColor(1);
                ImPlot::PlotBars("selected", &sx, &sy, 1, 0.6, selSpec);
            }
            ImPlot::EndPlot();
        }
    }

    // ---- Pose table. Rows are selectable and drive the shared pose selection. --------
    if (hasResult && ImGui::BeginTable("poses", 4,
                                       ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Rank");
        ImGui::TableSetupColumn("Affinity (kcal/mol)");
        ImGui::TableSetupColumn("RMSD l.b. (A)");
        ImGui::TableSetupColumn("RMSD u.b. (A)");
        ImGui::TableHeadersRow();
        for (int i = 0; i < static_cast<int>(d.poses.size()); ++i) {
            const auto& p = d.poses[i];
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            char lbl[32];
            std::snprintf(lbl, sizeof lbl, "%d##pose%d", p.rank, i);
            if (ImGui::Selectable(lbl, poseSel == i, ImGuiSelectableFlags_SpanAllColumns))
                poseSel = i;
            ImGui::TableSetColumnIndex(1); ImGui::TextUnformatted(f2(p.affinityKcalPerMol).c_str());
            ImGui::TableSetColumnIndex(2); ImGui::TextUnformatted(f2(p.rmsdLb).c_str());
            ImGui::TableSetColumnIndex(3); ImGui::TextUnformatted(f2(p.rmsdUb).c_str());
        }
        ImGui::EndTable();
    }

    // 3D pose overlay: the docked ligand conformer for the selected pose (real engine
    // poses carry their own coordinates; the labeled estimate carries the embedded
    // ligand as a stand-in so the viewport always has geometry). Engine is labeled.
    ImGui::Spacing();
    theme::sectionHeader("3D POSE / LIGAND");
    if (hasResult) {
        ImGui::SetNextItemWidth(160);
        ImGui::SliderInt("Pose##dock", &poseSel, 0, static_cast<int>(d.poses.size()) - 1);
        ImGui::SameLine();
        ImGui::TextDisabled("engine: %s%s", d.engine.c_str(),
                            d.real ? "" : "  (descriptor estimate, not a docked score)");
    }
    static std::string dockKey;
    static chem::Conformer dockConf;
    static chem::Conformer dockPocket;  // receptor binding-pocket overlay (real docks only)
    static ViewerUiState dockViewUi;
    // Data key (rebuilds the geometry) DOES include the pose, but the camera key passed
    // to the viewer below does NOT - so flipping between poses re-renders the new pose
    // without snapping the camera back to a re-fit (the poses share a coordinate frame).
    const std::string dataKey = m.id + "|" + st.dockTarget + "|" + std::to_string(poseSel) +
                                (hasResult ? "|r" : "|p") + (d.real ? "|e" : "");
    if (dataKey != dockKey) {
        dockKey = dataKey;
        dockConf = chem::Conformer{};
        dockPocket = chem::Conformer{};
        if (poseSel < static_cast<int>(d.poses.size()) && !d.poses[poseSel].ligand.empty())
            dockConf = d.poses[poseSel].ligand;
        else if (auto parsed = chem::parseSmiles(m.smiles))
            dockConf = chem::embed3D(*parsed);
        // Overlay the receptor pocket only for a REAL dock (the pose shares the
        // receptor's coordinate frame); the descriptor estimate is not receptor-aligned.
        if (d.real && !dockConf.empty()) {
            if (const ReceptorTarget* preset = docking::findPreset(st.dockTarget)) {
                const auto rp = docking::locatePreparedReceptor(preset->id);
                if (rp.ready) dockPocket = loadReceptorPocket(rp.path, dockConf, 5.5);
            }
        }
    }
    // Camera-fit key: stable across pose changes (molecule + target + render mode only).
    const std::string camKey = m.id + "|" + st.dockTarget + (d.real ? "|e" : "|p");
    molViewer3D(shell, dockConf, camKey, dockViewUi, 280.0f,
                dockPocket.empty() ? nullptr : &dockPocket);
    if (!dockPocket.empty())
        ImGui::TextDisabled("Receptor binding pocket shown (muted); toggle with the Receptor checkbox.");
    // Screenshot/automation hook (mirrors STIMLAB_PANEL/TARGET): scroll the pose viewer
    // into view for capture tooling. Harmless in normal use.
    static const bool kScroll3d = std::getenv("STIMLAB_DOCK_SCROLL3D") != nullptr;
    if (kScroll3d) ImGui::SetScrollHereY(1.0f);

    ImGui::Spacing();
    ImGui::TextDisabled("Note: binding affinity is a target-engagement signal, never a make-it signal.");
}

// ------------------------------------------------------------------- Workflows
void workflows(AppShell& shell) {
    const Molecule m = shell.currentMolecule();
    Services& s = shell.services();
    UiState& st = shell.state();
    if (!s.docking) return;

    auto targets = s.docking->targets();
    if (st.dockTarget.empty() && !targets.empty()) st.dockTarget = targets.front();

    ImGui::TextWrapped(
        "Run the prep->dock pipeline as a content-cached, cancellable DAG. Nodes run on a worker "
        "thread (the UI never blocks); re-running unchanged inputs is an instant cache hit, and "
        "provisioning a receptor/engine re-runs only the affected nodes.");
    ImGui::Spacing();

    ImGui::TextDisabled("TARGET (binding/pharmacology only)");
    ImGui::SetNextItemWidth(360);
    if (ImGui::BeginCombo("##wftarget", st.dockTarget.c_str())) {
        for (const auto& t : targets) {
            const bool sel = (t == st.dockTarget);
            if (ImGui::Selectable(t.c_str(), sel)) st.dockTarget = t;
            if (sel) ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }
    ImGui::Spacing();

    // Receptor readiness for the selected target: a prep->dock run still completes for
    // an unprepared target but yields the labeled estimate, so offer an on-demand
    // provision (the same single-target path as the Docking panel) before the run.
    {
        auto& prov = shell.provisioner();
        static std::string wfReadyKey;
        static bool wfSelReady = false;
        static bool wfWasRunning = false;
        const bool provRunning = prov.running();
        if (st.dockTarget != wfReadyKey || (wfWasRunning && !provRunning)) {
            wfReadyKey = st.dockTarget;
            wfSelReady = shell.receptorReady(st.dockTarget);
        }
        wfWasRunning = provRunning;
        ImGui::TextDisabled("receptor for %s: %s", st.dockTarget.c_str(),
                            wfSelReady ? "prepared (real dock)" : "not prepared (estimate)");
        if (provRunning) {
            ImGui::TextDisabled("provisioning: %s", prov.status().c_str());
        } else if (!wfSelReady) {
            const std::string lbl = "Provision " + st.dockTarget;
            if (ImGui::Button(lbl.c_str())) shell.provisionTarget(st.dockTarget);
        }
    }
    ImGui::Spacing();

    const auto snap = shell.workflowSnapshot();
    if (snap.running) {
        if (ImGui::Button("Cancel run")) shell.cancelWorkflow();
        ImGui::SameLine();
        ImGui::TextDisabled("running %s ...", snap.label.c_str());
    } else {
        if (ImGui::Button("Run workflow")) {
            std::string tid = st.dockTarget;
            if (const auto* p = docking::findPreset(st.dockTarget)) tid = p->id;
            shell.runWorkflow(m.smiles, tid, m.name + " -> " + st.dockTarget);
        }
        ImGui::SameLine();
        ImGui::TextDisabled("ligand: %s", m.name.c_str());
    }
    ImGui::Spacing();

    // ---- live DAG diagram ----------------------------------------------------
    auto colorOf = [](workflow::NodeStatus s) -> ImU32 {
        using NS = workflow::NodeStatus;
        switch (s) {
            case NS::Done:      return IM_COL32(76, 199, 107, 255);   // green
            case NS::Cached:    return IM_COL32(76, 179, 189, 255);   // teal
            case NS::Running:   return IM_COL32(92, 158, 245, 255);   // blue
            case NS::Failed:    return IM_COL32(230, 92, 86, 255);    // red
            case NS::Cancelled:
            case NS::Skipped:   return IM_COL32(230, 158, 71, 255);   // orange
            default:            return IM_COL32(115, 122, 140, 255);  // gray
        }
    };

    const auto& nodes = snap.nodes;
    if (!nodes.empty()) {
        // Column = topological depth (longest path from a root); rows stack per column.
        std::unordered_map<std::string, int> idx;
        for (int i = 0; i < static_cast<int>(nodes.size()); ++i) idx[nodes[i].id] = i;
        std::vector<int> depth(nodes.size(), 0);
        for (size_t pass = 0; pass < nodes.size(); ++pass)
            for (size_t i = 0; i < nodes.size(); ++i)
                for (const auto& d : nodes[i].deps) {
                    auto it = idx.find(d);
                    if (it != idx.end())
                        depth[i] = std::max(depth[i], depth[it->second] + 1);
                }
        int maxDepth = 0;
        for (int d : depth) maxDepth = std::max(maxDepth, d);
        std::vector<int> rowInCol(maxDepth + 1, 0);
        std::vector<ImVec2> centerL(nodes.size()), centerR(nodes.size()), p0(nodes.size()),
            p1(nodes.size());

        const float boxW = 168.0f, boxH = 52.0f, hgap = 64.0f, vgap = 20.0f;
        const ImVec2 origin = ImGui::GetCursorScreenPos();
        int maxRows = 0;
        for (size_t i = 0; i < nodes.size(); ++i) {
            const int col = depth[i];
            const int row = rowInCol[col]++;
            maxRows = std::max(maxRows, row + 1);
            const float x = origin.x + col * (boxW + hgap);
            const float y = origin.y + row * (boxH + vgap);
            p0[i] = ImVec2(x, y);
            p1[i] = ImVec2(x + boxW, y + boxH);
            centerL[i] = ImVec2(x, y + boxH * 0.5f);
            centerR[i] = ImVec2(x + boxW, y + boxH * 0.5f);
        }

        ImDrawList* dl = ImGui::GetWindowDrawList();
        // Edges first (under the boxes).
        for (size_t i = 0; i < nodes.size(); ++i)
            for (const auto& d : nodes[i].deps) {
                auto it = idx.find(d);
                if (it == idx.end()) continue;
                const ImVec2 a = centerR[it->second], b = centerL[i];
                dl->AddBezierCubic(a, ImVec2(a.x + hgap * 0.6f, a.y),
                                   ImVec2(b.x - hgap * 0.6f, b.y), b,
                                   IM_COL32(120, 128, 146, 220), 2.0f);
            }
        // Boxes.
        const double t = ImGui::GetTime();
        for (size_t i = 0; i < nodes.size(); ++i) {
            const ImU32 col = colorOf(nodes[i].status);
            float thick = 2.0f;
            if (nodes[i].status == workflow::NodeStatus::Running)
                thick = 2.0f + 1.6f * static_cast<float>(0.5 + 0.5 * std::sin(t * 6.0));
            dl->AddRectFilled(p0[i], p1[i], IM_COL32(26, 30, 38, 255), 6.0f);
            dl->AddRect(p0[i], p1[i], col, 6.0f, 0, thick);
            const ImVec2 idsz = ImGui::CalcTextSize(nodes[i].id.c_str());
            dl->AddText(ImVec2(p0[i].x + (boxW - idsz.x) * 0.5f, p0[i].y + 9.0f),
                        IM_COL32(232, 235, 240, 255), nodes[i].id.c_str());
            const char* slab = workflow::toString(nodes[i].status);
            const ImVec2 ssz = ImGui::CalcTextSize(slab);
            dl->AddText(ImVec2(p0[i].x + (boxW - ssz.x) * 0.5f, p0[i].y + 29.0f), col, slab);
        }
        ImGui::Dummy(ImVec2((maxDepth + 1) * (boxW + hgap), maxRows * (boxH + vgap) + 4.0f));
    }

    // ---- legend + summary + node details ------------------------------------
    ImGui::Spacing();
    if (snap.hasResult && !snap.running) {
        ImGui::Text("Last run: ran %d  -  cached %d  -  failed %d  -  cancelled %d", snap.ran,
                    snap.cached, snap.failed, snap.cancelled);
        if (snap.cached > 0 && snap.ran == 0)
            ImGui::TextDisabled("(fully resumed from cache - no node re-executed)");
    }
    ImGui::Spacing();
    if (ImGui::BeginTable("wfnodes", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Node", ImGuiTableColumnFlags_WidthFixed, 130.0f);
        ImGui::TableSetupColumn("Status", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Detail / output");
        ImGui::TableHeadersRow();
        for (const auto& n : snap.nodes) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(n.id.c_str());
            ImGui::TableSetColumnIndex(1);
            ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(colorOf(n.status)), "%s",
                               workflow::toString(n.status));
            ImGui::TableSetColumnIndex(2); ImGui::TextWrapped("%s", n.detail.c_str());
        }
        ImGui::EndTable();
    }
    ImGui::Spacing();
    ImGui::TextDisabled("Nodes are cached by hash(module, version, inputs) under "
                        "%%APPDATA%%/StimLab/cache. Analysis only - no synthesis content.");
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
    const auto rows = s.runs->recent();
    ImGui::TextDisabled("%zu run(s) - persisted to SQLite under %%APPDATA%%/StimLab/stimlab.db",
                        rows.size());
    ImGui::Spacing();
    if (ImGui::BeginTable("runs", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("ID", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Kind", ImGuiTableColumnFlags_WidthFixed, 110.0f);
        ImGui::TableSetupColumn("Subject");
        ImGui::TableSetupColumn("Status", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Summary");
        ImGui::TableHeadersRow();
        for (const auto& r : rows) {
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
        "(each with a real PDB reference + a binding-site box) is a built-in C++ table; receptor "
        "PDBQTs are prepared on demand into %APPDATA%/StimLab/runtime/receptors.");
    ImGui::Separator();
    if (s.docking) {
        for (const auto& t : s.docking->targets()) ImGui::BulletText("%s", t.c_str());
    }
}

// -------------------------------------------------------------------- Settings
void settings(AppShell& shell) {
    ImGui::TextDisabled("AI ASSISTANT");
    ImGui::TextWrapped(
        "The assistant explains panels, reads the selected compound's real structure-derived "
        "properties, and navigates/highlights the UI. It never provides synthesis, route, or "
        "manufacturability guidance - that is out of scope by design.");
    ImGui::Spacing();

    int provider = shell.agentProviderIndex();
    const char* providers[] = {"Anthropic (Claude)", "Offline (no API key)"};
    ImGui::SetNextItemWidth(280);
    if (ImGui::Combo("Provider", &provider, providers, IM_ARRAYSIZE(providers)))
        shell.setAgentProviderIndex(provider);

    if (!shell.anthropicTransport()) {
        ImGui::PushStyleColor(ImGuiCol_Text, theme::verdictColor(2));
        ImGui::TextWrapped(
            "This build has no networking (default 'windows' preset). Rebuild with the "
            "'windows-science' preset to enable the live Anthropic provider; the offline assistant "
            "still navigates and explains.");
        ImGui::PopStyleColor();
    }

    // Model (free text so any model the key supports can be used).
    static char modelBuf[96];
    static bool modelInit = false;
    if (!modelInit) {
        std::snprintf(modelBuf, sizeof(modelBuf), "%s", shell.agentModel().c_str());
        modelInit = true;
    }
    ImGui::SetNextItemWidth(280);
    ImGui::InputText("Model", modelBuf, sizeof(modelBuf));
    ImGui::SameLine();
    if (ImGui::Button("Set##model")) shell.setAgentModel(modelBuf);
    ImGui::TextDisabled("Default claude-opus-4-8; claude-haiku-4-5 is faster/cheaper for UI help.");
    ImGui::Spacing();

    // API key - encrypted at rest via DPAPI; plaintext never persisted or shown.
    ImGui::TextDisabled("API KEY (encrypted at rest via Windows DPAPI)");
    static char keyBuf[256] = {0};
    ImGui::SetNextItemWidth(360);
    ImGui::InputTextWithHint(
        "##key", shell.hasApiKey() ? "A key is stored - type to replace..." : "Paste Anthropic API key...",
        keyBuf, sizeof(keyBuf), ImGuiInputTextFlags_Password);
    ImGui::SameLine();
    if (ImGui::Button("Save key") && keyBuf[0] != '\0') {
        shell.saveApiKey(keyBuf);
        for (char& c : keyBuf) c = '\0';  // wipe plaintext from the input buffer
    }
    if (shell.hasApiKey()) {
        ImGui::SameLine();
        if (ImGui::Button("Clear key")) shell.clearApiKey();
    }

    if (shell.anthropicReady()) {
        ImGui::PushStyleColor(ImGuiCol_Text, theme::verdictColor(1));
        ImGui::Text("Live provider active: %s", shell.activeProviderLabel().c_str());
        ImGui::PopStyleColor();
    } else {
        ImGui::TextDisabled("Active provider: %s", shell.activeProviderLabel().c_str());
    }
    ImGui::Spacing();

    ImGui::TextDisabled("BEHAVIOR");
    bool autop = shell.autopilot();
    if (ImGui::Checkbox("Autopilot - run navigate/highlight tools automatically", &autop))
        shell.setAutopilot(autop);
    ImGui::Spacing();

    ImGui::TextDisabled("COMPUTE (docking engine)");
    int compute = shell.computeMode();
    bool computeChanged = false;
    computeChanged |= ImGui::RadioButton("Auto", &compute, 0); ImGui::SameLine();
    computeChanged |= ImGui::RadioButton("GPU", &compute, 1); ImGui::SameLine();
    computeChanged |= ImGui::RadioButton("CPU", &compute, 2);
    if (computeChanged) shell.setComputeMode(compute);
#ifdef STIMLAB_HAVE_CUDA
    ImGui::TextWrapped("Auto: AutoDock Vina (CPU) first, then the GPU engines. GPU: Vina-GPU "
                       "(OpenCL - the real Vina search on the GPU) when provisioned, else the "
                       "first-party CUDA rigid-grid dock on this NVIDIA GPU. CPU: Vina/smina only.");
#else
    ImGui::TextWrapped("Auto: AutoDock Vina (CPU) first, then the GPU engine. GPU: Vina-GPU "
                       "(OpenCL - the real Vina search on the GPU) when provisioned, else the labeled "
                       "descriptor estimate (the first-party CUDA dock needs the windows-cuda build). "
                       "CPU: Vina/smina only.");
#endif
    // Vina-GPU (OpenCL) is a self-contained subprocess engine available in EVERY build;
    // provisioning downloads its binaries + compiles a kernel for this GPU (off-thread).
    {
        const bool busy = shell.vinaGpuProvisioning();
        const bool vgReady = shell.vinaGpuReady();  // cache-only fs check
        ImGui::BeginDisabled(busy);
        if (ImGui::Button(vgReady ? "Re-provision Vina-GPU (OpenCL)"
                                  : "Provision Vina-GPU (OpenCL GPU engine)"))
            shell.provisionVinaGpu();
        ImGui::EndDisabled();
        ImGui::SameLine();
        if (vgReady) {
            ImGui::PushStyleColor(ImGuiCol_Text, theme::verdictColor(1));
            ImGui::TextUnformatted("ready");
            ImGui::PopStyleColor();
        } else {
            ImGui::TextDisabled(busy ? "working..." : "not provisioned");
        }
        ImGui::TextDisabled("%s", shell.vinaGpuStatus().c_str());
    }
    ImGui::Spacing();

    ImGui::TextDisabled("STORAGE");
    ImGui::TextWrapped("All state lives under %%APPDATA%%/StimLab (db, artifacts, runtime, presets, logs).");
    ImGui::Spacing();

    ImGui::TextDisabled("RUNTIME (self-provisioned components)");
    auto fmtManifest = [](const ManifestStatus& st) {
        if (st.total == 0)
            return std::string("Nothing provisioned yet - use the Docking panel's Provision button.");
        std::string s = std::to_string(st.present) + "/" + std::to_string(st.total) +
                        " components verified (size + content hash)";
        if (!st.missing.empty()) s += ", " + std::to_string(st.missing.size()) + " missing";
        if (!st.corrupt.empty()) s += ", " + std::to_string(st.corrupt.size()) + " corrupt (healed)";
        return s + ".";
    };
    static std::string mfStatus;
    static bool mfInit = false;
    if (!mfInit) {
        mfInit = true;
        mfStatus = fmtManifest(Manifest::load(AppPaths::instance().manifest()).verify());
    }
    ImGui::TextWrapped("%s", mfStatus.c_str());
    if (ImGui::Button("Verify + heal runtime")) mfStatus = fmtManifest(docking::selfHealManifest());
    ImGui::SameLine();
    ImGui::TextDisabled("Deletes any corrupt engine/receptor so it re-provisions; manifest.json is the source of truth.");
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

// ===================================================================== Sketcher
// A lightweight 2D structure editor: the user places atoms and bonds on a canvas
// and the graph is serialised to SMILES, which flows through the SAME tested
// parser -> descriptor/embed3D pipeline as every other compound (so nothing is
// faked - every predicted number is computed from the drawn structure).
namespace {

struct SkAtom { ImVec2 p; int z = 6; };
struct SkBond { int a = 0, b = 0; int order = 1; bool aromatic = false; };

struct Sketch {
    std::vector<SkAtom> atoms;
    std::vector<SkBond> bonds;
    int  tool = 0;          // 0 Atom, 1 Bond, 2 Move, 3 Erase
    int  element = 6;       // element Z placed / retyped onto atoms
    int  order = 1;         // current bond order
    bool aromatic = false;  // current bond is aromatic
    int  pendingBond = -1;  // first-picked atom while bonding
    int  dragAtom = -1;     // atom being dragged in Move mode
    std::string smiles;     // last serialised SMILES
};

struct SkElem { int z; const char* sym; };
const SkElem kSkElems[] = {{6, "C"}, {7, "N"},  {8, "O"},  {16, "S"}, {15, "P"},
                           {9, "F"}, {17, "Cl"}, {35, "Br"}, {53, "I"}};

ImU32 skAtomColor(int z) {
    const std::uint32_t c = render::cpkColor(z);  // 0xAABBGGRR
    return IM_COL32(c & 0xFF, (c >> 8) & 0xFF, (c >> 16) & 0xFF, 0xFF);
}

void removeSketchAtom(Sketch& sk, int idx) {
    sk.bonds.erase(std::remove_if(sk.bonds.begin(), sk.bonds.end(),
                                  [&](const SkBond& b) { return b.a == idx || b.b == idx; }),
                   sk.bonds.end());
    for (auto& b : sk.bonds) { if (b.a > idx) --b.a; if (b.b > idx) --b.b; }
    sk.atoms.erase(sk.atoms.begin() + idx);
    if (sk.pendingBond == idx) sk.pendingBond = -1;
    else if (sk.pendingBond > idx) --sk.pendingBond;
}

// Templates (reused by the toolbar buttons + the screenshot/automation hook).
void sketchAddBenzene(Sketch& sk, ImVec2 ctr) {
    const int base = static_cast<int>(sk.atoms.size());
    const float r = 46.0f;
    for (int i = 0; i < 6; ++i) {
        const float a = static_cast<float>(i) / 6.0f * 6.2831853f - 1.5708f;
        sk.atoms.push_back({ImVec2(ctr.x + r * std::cos(a), ctr.y + r * std::sin(a)), 6});
    }
    for (int i = 0; i < 6; ++i) sk.bonds.push_back({base + i, base + (i + 1) % 6, 1, true});
}

void sketchSeedPea(Sketch& sk, float canvasH) {
    const ImVec2 ctr(130.0f, canvasH * 0.5f);
    const int base = static_cast<int>(sk.atoms.size());
    sketchAddBenzene(sk, ctr);
    const ImVec2 r0 = sk.atoms[base].p;
    const int ic1 = static_cast<int>(sk.atoms.size()); sk.atoms.push_back({ImVec2(r0.x + 36, r0.y - 20), 6});
    const int ic2 = static_cast<int>(sk.atoms.size()); sk.atoms.push_back({ImVec2(r0.x + 70, r0.y - 2), 6});
    const int iN  = static_cast<int>(sk.atoms.size()); sk.atoms.push_back({ImVec2(r0.x + 104, r0.y - 22), 7});
    sk.bonds.push_back({base, ic1, 1, false});
    sk.bonds.push_back({ic1, ic2, 1, false});
    sk.bonds.push_back({ic2, iN, 1, false});
}

// Serialise the sketch graph to SMILES (spanning-tree DFS + ring-closure digits).
std::string sketchToSmiles(const Sketch& sk) {
    const int n = static_cast<int>(sk.atoms.size());
    if (n == 0) return "";
    const int nb = static_cast<int>(sk.bonds.size());

    std::vector<std::vector<std::pair<int, int>>> adj(n);  // (neighbor, bondIdx)
    for (int bi = 0; bi < nb; ++bi) {
        const auto& b = sk.bonds[bi];
        if (b.a < 0 || b.a >= n || b.b < 0 || b.b >= n || b.a == b.b) continue;
        adj[b.a].push_back({b.b, bi});
        adj[b.b].push_back({b.a, bi});
    }
    std::vector<char> arom(n, 0);
    for (const auto& b : sk.bonds)
        if (b.aromatic && b.a >= 0 && b.a < n && b.b >= 0 && b.b < n) { arom[b.a] = 1; arom[b.b] = 1; }

    std::vector<char> visited(n, 0), treeBond(nb, 0), ringBond(nb, 0);
    std::vector<std::vector<std::pair<int, int>>> ringTok(n);  // (digit, bondIdx)
    int nextDigit = 1;
    std::function<void(int, int)> findRings = [&](int u, int parentBond) {
        visited[u] = 1;
        for (auto [v, bi] : adj[u]) {
            if (bi == parentBond) continue;
            if (!visited[v]) { treeBond[bi] = 1; findRings(v, bi); }
            else if (!ringBond[bi] && !treeBond[bi]) {
                ringBond[bi] = 1;
                const int d = nextDigit++;
                ringTok[u].push_back({d, bi});
                ringTok[v].push_back({d, bi});
            }
        }
    };
    for (int i = 0; i < n; ++i)
        if (!visited[i]) findRings(i, -1);

    auto digitTok = [](int d) -> std::string {
        if (d < 10) return std::string(1, static_cast<char>('0' + d));
        char b[8]; std::snprintf(b, sizeof b, "%%%02d", d); return b;
    };
    auto bondSym = [&](int bi) -> std::string {
        if (sk.bonds[bi].aromatic) return "";
        if (sk.bonds[bi].order == 2) return "=";
        if (sk.bonds[bi].order == 3) return "#";
        return "";
    };
    auto atomTok = [&](int u) -> std::string {
        const char* sym = chem::symbolByZ(sk.atoms[u].z);
        std::string s = sym ? sym : "C";
        if (arom[u] && s.size() == 1)
            s[0] = static_cast<char>(std::tolower(static_cast<unsigned char>(s[0])));
        return s;
    };

    std::vector<char> emitted(n, 0);
    std::string out;
    std::function<void(int)> emit = [&](int u) {
        emitted[u] = 1;
        out += atomTok(u);
        for (auto [d, bi] : ringTok[u]) { out += bondSym(bi); out += digitTok(d); }
        std::vector<std::pair<int, int>> kids;
        for (auto [v, bi] : adj[u])
            if (treeBond[bi] && !emitted[v]) kids.push_back({v, bi});
        for (size_t k = 0; k < kids.size(); ++k) {
            const int v = kids[k].first, bi = kids[k].second;
            if (emitted[v]) continue;
            const bool last = (k + 1 == kids.size());
            if (!last) out += "(";
            out += bondSym(bi);
            emit(v);
            if (!last) out += ")";
        }
    };
    for (int i = 0; i < n; ++i)
        if (!emitted[i]) { if (!out.empty()) out += "."; emit(i); }
    return out;
}

// Draw the sketch toolbar + canvas. Returns true when the connectivity (not just
// the 2D layout) changed, so the caller re-serialises + re-analyses.
bool drawSketchCanvas(Sketch& sk, float height) {
    bool changed = false;
    auto toggleBtn = [](const char* label, bool active) -> bool {
        if (active) {
            ImGui::PushStyleColor(ImGuiCol_Button, theme::kAccentSoft);
            ImGui::PushStyleColor(ImGuiCol_Text, theme::kAccent);
        }
        const bool clicked = ImGui::Button(label);
        if (active) ImGui::PopStyleColor(2);
        return clicked;
    };

    // --- toolbar ---
    ImGui::TextDisabled("Tool"); ImGui::SameLine();
    const char* tools[] = {"Atom", "Bond", "Move", "Erase"};
    for (int t = 0; t < 4; ++t) {
        ImGui::SameLine();
        if (toggleBtn(tools[t], sk.tool == t)) { sk.tool = t; sk.pendingBond = -1; }
    }

    ImGui::TextDisabled("Atom"); ImGui::SameLine();
    for (const auto& e : kSkElems) {
        ImGui::SameLine();
        if (toggleBtn(e.sym, sk.element == e.z)) sk.element = e.z;
    }

    ImGui::TextDisabled("Bond"); ImGui::SameLine();
    struct BT { const char* l; int order; bool ar; };
    const BT bts[] = {{"Single", 1, false}, {"Double", 2, false}, {"Triple", 3, false},
                      {"Aromatic", 1, true}};
    for (const auto& bt : bts) {
        ImGui::SameLine();
        const bool active = bt.ar ? sk.aromatic : (!sk.aromatic && sk.order == bt.order);
        if (toggleBtn(bt.l, active)) { sk.order = bt.order; sk.aromatic = bt.ar; }
    }

    if (ImGui::Button("Benzene")) { sketchAddBenzene(sk, ImVec2(150.0f, height * 0.5f)); changed = true; }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Drop an aromatic 6-ring");
    ImGui::SameLine();
    if (ImGui::Button("PEA core")) { sketchSeedPea(sk, height); changed = true; }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Phenethylamine starter: ring-CH2-CH2-NH2");
    ImGui::SameLine();
    if (ImGui::Button("Clear")) { sk.atoms.clear(); sk.bonds.clear(); sk.pendingBond = -1; changed = true; }

    // --- canvas ---
    const ImVec2 p0 = ImGui::GetCursorScreenPos();
    const float w = std::max(120.0f, ImGui::GetContentRegionAvail().x);
    const ImVec2 size(w, height);
    ImGui::InvisibleButton("##skcanvas", size);
    const bool hovered = ImGui::IsItemHovered();
    const bool active = ImGui::IsItemActive();
    const ImVec2 p1(p0.x + size.x, p0.y + size.y);
    ImGuiIO& io = ImGui::GetIO();
    const ImVec2 m(io.MousePos.x - p0.x, io.MousePos.y - p0.y);

    // hit-test atoms (then bonds).
    int hitAtom = -1; float bestD = 16.0f * 16.0f;
    for (int i = 0; i < static_cast<int>(sk.atoms.size()); ++i) {
        const float dx = sk.atoms[i].p.x - m.x, dy = sk.atoms[i].p.y - m.y;
        const float d = dx * dx + dy * dy;
        if (d < bestD) { bestD = d; hitAtom = i; }
    }
    auto distSeg = [](ImVec2 q, ImVec2 a, ImVec2 b) -> float {
        const float vx = b.x - a.x, vy = b.y - a.y, wx = q.x - a.x, wy = q.y - a.y;
        const float c1 = vx * wx + vy * wy;
        if (c1 <= 0) return std::hypot(q.x - a.x, q.y - a.y);
        const float c2 = vx * vx + vy * vy;
        if (c2 <= c1) return std::hypot(q.x - b.x, q.y - b.y);
        const float t = c1 / c2;
        return std::hypot(q.x - (a.x + t * vx), q.y - (a.y + t * vy));
    };
    int hitBond = -1;
    if (hitAtom < 0) {
        float bd = 8.0f;
        for (int bi = 0; bi < static_cast<int>(sk.bonds.size()); ++bi) {
            const auto& b = sk.bonds[bi];
            const float d = distSeg(m, sk.atoms[b.a].p, sk.atoms[b.b].p);
            if (d < bd) { bd = d; hitBond = bi; }
        }
    }

    // interaction.
    if (hovered && ImGui::IsMouseClicked(0)) {
        if (sk.tool == 0) {  // Atom
            if (hitAtom >= 0) { if (sk.atoms[hitAtom].z != sk.element) { sk.atoms[hitAtom].z = sk.element; changed = true; } }
            else { sk.atoms.push_back({m, sk.element}); changed = true; }
        } else if (sk.tool == 1) {  // Bond
            if (hitAtom >= 0) {
                if (sk.pendingBond < 0) sk.pendingBond = hitAtom;
                else if (sk.pendingBond == hitAtom) sk.pendingBond = -1;
                else {
                    int existing = -1;
                    for (int bi = 0; bi < static_cast<int>(sk.bonds.size()); ++bi) {
                        const auto& b = sk.bonds[bi];
                        if ((b.a == sk.pendingBond && b.b == hitAtom) ||
                            (b.a == hitAtom && b.b == sk.pendingBond)) { existing = bi; break; }
                    }
                    if (existing >= 0) { sk.bonds[existing].order = sk.order; sk.bonds[existing].aromatic = sk.aromatic; }
                    else sk.bonds.push_back({sk.pendingBond, hitAtom, sk.order, sk.aromatic});
                    sk.pendingBond = -1; changed = true;
                }
            } else sk.pendingBond = -1;
        } else if (sk.tool == 2) {  // Move
            if (hitAtom >= 0) sk.dragAtom = hitAtom;
        } else if (sk.tool == 3) {  // Erase
            if (hitAtom >= 0) { removeSketchAtom(sk, hitAtom); changed = true; }
            else if (hitBond >= 0) { sk.bonds.erase(sk.bonds.begin() + hitBond); changed = true; }
        }
    }
    if (sk.tool == 2 && sk.dragAtom >= 0) {
        if (active && ImGui::IsMouseDown(0))
            sk.atoms[sk.dragAtom].p = ImVec2(std::clamp(m.x, 8.0f, size.x - 8.0f),
                                             std::clamp(m.y, 8.0f, size.y - 8.0f));
        else sk.dragAtom = -1;
    }

    // render.
    ImDrawList* dl = ImGui::GetWindowDrawList();
    dl->PushClipRect(p0, p1, true);
    dl->AddRectFilled(p0, p1, theme::kSurface, 9.0f);
    dl->AddRect(p0, p1, theme::kBorder, 9.0f);

    auto unit = [](ImVec2 d) -> ImVec2 {
        const float l = std::hypot(d.x, d.y);
        return (l < 1e-3f) ? ImVec2(0, 0) : ImVec2(d.x / l, d.y / l);
    };
    for (int bi = 0; bi < static_cast<int>(sk.bonds.size()); ++bi) {
        const auto& b = sk.bonds[bi];
        const ImVec2 a(p0.x + sk.atoms[b.a].p.x, p0.y + sk.atoms[b.a].p.y);
        const ImVec2 c(p0.x + sk.atoms[b.b].p.x, p0.y + sk.atoms[b.b].p.y);
        const ImU32 col = (bi == hitBond) ? theme::kAccent : IM_COL32(178, 188, 202, 255);
        const ImVec2 dir = unit(ImVec2(c.x - a.x, c.y - a.y));
        const ImVec2 pp(-dir.y, dir.x);
        const float o = 3.2f;
        if (b.aromatic) {
            dl->AddLine(a, c, col, 2.0f);
            dl->AddLine(ImVec2(a.x + pp.x * o + dir.x * 8, a.y + pp.y * o + dir.y * 8),
                        ImVec2(c.x + pp.x * o - dir.x * 8, c.y + pp.y * o - dir.y * 8), col, 1.4f);
        } else if (b.order == 2) {
            dl->AddLine(ImVec2(a.x + pp.x * o, a.y + pp.y * o), ImVec2(c.x + pp.x * o, c.y + pp.y * o), col, 2.0f);
            dl->AddLine(ImVec2(a.x - pp.x * o, a.y - pp.y * o), ImVec2(c.x - pp.x * o, c.y - pp.y * o), col, 2.0f);
        } else if (b.order == 3) {
            dl->AddLine(a, c, col, 2.0f);
            dl->AddLine(ImVec2(a.x + pp.x * o * 1.6f, a.y + pp.y * o * 1.6f), ImVec2(c.x + pp.x * o * 1.6f, c.y + pp.y * o * 1.6f), col, 1.6f);
            dl->AddLine(ImVec2(a.x - pp.x * o * 1.6f, a.y - pp.y * o * 1.6f), ImVec2(c.x - pp.x * o * 1.6f, c.y - pp.y * o * 1.6f), col, 1.6f);
        } else {
            dl->AddLine(a, c, col, 2.0f);
        }
    }
    for (int i = 0; i < static_cast<int>(sk.atoms.size()); ++i) {
        const ImVec2 c(p0.x + sk.atoms[i].p.x, p0.y + sk.atoms[i].p.y);
        const int z = sk.atoms[i].z;
        const char* sym = chem::symbolByZ(z);
        if (!sym) sym = "C";
        if (i == sk.pendingBond) dl->AddCircle(c, 14.0f, theme::kAccent, 18, 2.5f);
        if (i == hitAtom && hovered) dl->AddCircle(c, 14.0f, IM_COL32(255, 255, 255, 90), 18, 1.5f);
        if (z == 6) {
            dl->AddCircleFilled(c, 4.0f, IM_COL32(150, 160, 175, 255));
        } else {
            const ImVec2 ts = ImGui::CalcTextSize(sym);
            const float rad = std::max(10.0f, ts.x * 0.5f + 5.0f);
            dl->AddCircleFilled(c, rad, skAtomColor(z));
            dl->AddText(ImVec2(c.x - ts.x * 0.5f, c.y - ts.y * 0.5f), IM_COL32(15, 19, 26, 255), sym);
        }
    }
    if (sk.atoms.empty()) {
        const char* hint = "Click to place atoms - then use Bond to connect them, or drop a template above.";
        const ImVec2 ts = ImGui::CalcTextSize(hint);
        dl->AddText(ImVec2(p0.x + (size.x - ts.x) * 0.5f, p0.y + size.y * 0.5f - ts.y),
                    theme::kTextDim, hint);
    }
    dl->PopClipRect();

    ImGui::TextDisabled("%zu atoms, %zu bonds   (drag with Move; click two atoms in Bond mode)",
                        sk.atoms.size(), sk.bonds.size());
    return changed;
}

}  // namespace

// ------------------------------------------------------------ Analog Explorer
void analogExplorer(AppShell& shell) {
    Services& s = shell.services();
    if (!s.library || !s.stability || !s.absorption || !s.admet || !s.similarity || !s.legal) return;
    const auto lib = s.library->all();

    static int   mode = 0;                          // 0 = model by properties, 1 = draw
    static int   parentIdx = 1;                      // amphetamine
    static float logP = 1.76f, tpsa = 26.0f, mw = 135.0f;
    static int   hbd = 1, hba = 1, rot = 2;
    static bool  ester = false, catechol = false, arylKetone = false, mdoxy = false;
    static int   lastParent = -999;
    static Sketch sketch;
    static ViewerUiState aeViewer;
    static std::string   aeConfKey;
    static chem::Conformer aeConf;

    // Screenshot/automation hook (mirrors STIMLAB_PANEL/TARGET): STIMLAB_ANALOG_DRAW
    // opens straight to the sketcher seeded with a phenethylamine core. Harmless in
    // normal use; only the capture tooling sets it.
    static const bool kAutoDraw = std::getenv("STIMLAB_ANALOG_DRAW") != nullptr;
    static bool autoDrawApplied = false;
    if (kAutoDraw && !autoDrawApplied && sketch.atoms.empty()) sketchSeedPea(sketch, 300.0f);

    ImGui::TextWrapped(
        "Design a candidate analog - model it from physicochemical knobs or DRAW the structure "
        "atom-by-atom - then preview it in 3D and screen it against existing samples and predicted "
        "byproducts. (Analysis only - no synthesis guidance.)");
    ImGui::Spacing();

    if (ImGui::BeginTable("aelayout", 2,
                          ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_Resizable)) {
        ImGui::TableSetupColumn("design", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("preview", ImGuiTableColumnFlags_WidthFixed, 360.0f);

        Molecule c;
        c.id = "candidate";
        c.name = "Candidate analog";
        bool haveGraph = false;

        // ---- column 0: define the candidate (two input modes) ----------------
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        if (ImGui::BeginTabBar("aemode")) {
            if (ImGui::BeginTabItem("Model by properties")) {
                mode = 0;
                ImGui::Spacing();
                ImGui::SetNextItemWidth(240);
                const std::string pcur = (parentIdx >= 0 && parentIdx < static_cast<int>(lib.size()))
                                             ? lib[parentIdx].name : "(scratch)";
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
                ImGui::Spacing();
                if (ImGui::BeginTable("aeknobs", 2, ImGuiTableFlags_BordersInnerV)) {
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0);
                    theme::sectionHeader("PROPERTIES");
                    ImGui::SetNextItemWidth(200); ImGui::SliderFloat("logP", &logP, -3.0f, 5.0f, "%.2f");
                    ImGui::SetNextItemWidth(200); ImGui::SliderFloat("TPSA", &tpsa, 0.0f, 150.0f, "%.0f");
                    ImGui::SetNextItemWidth(200); ImGui::SliderFloat("MW", &mw, 80.0f, 400.0f, "%.0f");
                    ImGui::SetNextItemWidth(200); ImGui::SliderInt("H-bond donors", &hbd, 0, 6);
                    ImGui::SetNextItemWidth(200); ImGui::SliderInt("H-bond acceptors", &hba, 0, 10);
                    ImGui::SetNextItemWidth(200); ImGui::SliderInt("Rotatable bonds", &rot, 0, 12);
                    ImGui::TableSetColumnIndex(1);
                    theme::sectionHeader("FUNCTIONAL GROUPS");
                    ImGui::Checkbox("Ester (hydrolysis-labile)", &ester);
                    ImGui::Checkbox("Catechol (oxidation / COMT)", &catechol);
                    ImGui::Checkbox("Aryl ketone (beta-keto)", &arylKetone);
                    ImGui::Checkbox("Methylenedioxy ring", &mdoxy);
                    ImGui::EndTable();
                }
                ImGui::EndTabItem();
            }
            ImGuiTabItemFlags drawFlags = 0;
            if (kAutoDraw && !autoDrawApplied) {
                drawFlags = ImGuiTabItemFlags_SetSelected;
                autoDrawApplied = true;
            }
            if (ImGui::BeginTabItem("Draw structure", nullptr, drawFlags)) {
                mode = 1;
                ImGui::Spacing();
                drawSketchCanvas(sketch, 300.0f);
                sketch.smiles = sketchToSmiles(sketch);
                ImGui::EndTabItem();
            }
            ImGui::EndTabBar();
        }

        // ---- assemble the candidate from the active mode --------------------
        // Draw mode runs the sketch through the real SMILES -> descriptor pipeline,
        // so every predicted number is computed from the drawn structure.
        if (mode == 0) {
            c.logP = logP; c.tpsa = tpsa; c.molWeight = mw;
            c.hbd = hbd; c.hba = hba; c.rotatableBonds = rot;
            c.formula = "(modeled)";
            c.drugClass = (parentIdx >= 0 && parentIdx < static_cast<int>(lib.size()))
                              ? lib[parentIdx].drugClass : "Phenethylamine (candidate)";
            const std::string ring =
                catechol ? "c1ccc(O)c(O)c1" : (mdoxy ? "c1ccc2OCOc2c1" : "c1ccccc1");
            c.smiles = ring + "CC(N)C" + (ester ? "OC(=O)C" : "") + (arylKetone ? "C(=O)c1ccccc1" : "");
            if (auto pm = chem::parseSmiles(c.smiles)) { c.formula = chem::molecularFormula(*pm); }
            haveGraph = true;
        } else {
            c.smiles = sketch.smiles;
            c.drugClass = "User-drawn structure";
            if (auto pm = chem::parseSmiles(c.smiles); pm && !pm->empty()) {
                c.formula = chem::molecularFormula(*pm);
                c.molWeight = chem::molecularWeight(*pm);
                c.logP = chem::crippenLogP(*pm);
                c.tpsa = chem::tpsa(*pm);
                c.hbd = chem::hbdCount(*pm);
                c.hba = chem::hbaCount(*pm);
                c.rotatableBonds = chem::rotatableBondCount(*pm);
                haveGraph = true;
            }
        }

        // ---- column 1: live chemical preview --------------------------------
        ImGui::TableSetColumnIndex(1);
        theme::sectionHeader("STRUCTURE PREVIEW");
        if (haveGraph) {
            if (c.smiles != aeConfKey) {
                aeConfKey = c.smiles;
                aeConf = chem::Conformer{};
                if (auto pm = chem::parseSmiles(c.smiles)) aeConf = chem::embed3D(*pm);
            }
            if (!molViewer3D(shell, aeConf, "ae:" + c.smiles, aeViewer, 220.0f))
                moleculeSchematic(c, ImVec2(320, 200));
            ImGui::Spacing();
            ImGui::TextDisabled("Formula  %s", c.formula.c_str());
            ImGui::PushTextWrapPos(0.0f);
            ImGui::TextDisabled("SMILES  %s", c.smiles.c_str());
            ImGui::PopTextWrapPos();
            ImGui::Spacing();
            if (ImGui::Button("Set as active compound", ImVec2(-1, 0))) {
                Molecule act = c;
                act.id = "__custom__";
                act.name = (mode == 1) ? "Drawn analog" : "Modeled analog";
                act.legalStatus = "(unscheduled / unknown)";
                act.notes = "Created in the Analog Explorer.";
                if (auto pm = chem::parseSmiles(act.smiles)) {
                    act.formula = chem::molecularFormula(*pm);
                    act.molWeight = chem::molecularWeight(*pm);
                    act.logP = chem::crippenLogP(*pm);
                    act.tpsa = chem::tpsa(*pm);
                    act.hbd = chem::hbdCount(*pm);
                    act.hba = chem::hbaCount(*pm);
                    act.rotatableBonds = chem::rotatableBondCount(*pm);
                }
                shell.state().customMolecule = act;
                shell.state().hasCustom = true;
                shell.state().selectedMolecule = "__custom__";
                shell.requestHighlight("Structure",
                    "Your analog is now the active compound across every panel.");
            }
        } else {
            ImGui::Dummy(ImVec2(0, 24));
            ImGui::TextColored(theme::verdictColor(2), "Draw atoms and bonds to build a structure.");
            ImGui::TextDisabled("Tip: drop a Benzene or PEA core, then edit it with the Atom/Bond tools.");
        }
        ImGui::EndTable();

        if (!haveGraph) return;  // empty/invalid sketch - nothing to analyze yet

        // ---- analysis (full width, computed from `c`) -----------------------
        const auto st = s.stability->analyze(c);
        const auto ab = s.absorption->predict(c);
        const auto ad = s.admet->screen(c);
        const auto si = s.similarity->search(c);
        const auto lg = s.legal->score(c);

        theme::sectionHeader("PREDICTED PROFILE");
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

        if (ImGui::BeginTable("aescreen", 2, ImGuiTableFlags_BordersInnerV)) {
            ImGui::TableSetupColumn("near", ImGuiTableColumnFlags_WidthStretch);
            ImGui::TableSetupColumn("legal", ImGuiTableColumnFlags_WidthStretch);
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            theme::sectionHeader("NEAREST KNOWN SAMPLE");
            if (!si.hits.empty()) {
                ImGui::TextColored(theme::verdictColor(1), "%s", si.nearestName.c_str());
                ImGui::TextDisabled("Tanimoto %.2f  -  %s", si.nearestScore,
                                    si.hits.front().referenceClass.c_str());
                ImGui::TextDisabled("%s", si.hits.front().legalStatus.c_str());
            } else {
                ImGui::TextDisabled("No close reference.");
            }
            ImGui::TableSetColumnIndex(1);
            theme::sectionHeader("LEGAL-ANALOG SCORE");
            const Verdict lband = lg.substantialSimilarity >= 75 ? Verdict::Danger
                                  : (lg.substantialSimilarity >= 50 ? Verdict::Warn : Verdict::Good);
            ImGui::PushStyleColor(ImGuiCol_PlotHistogram, theme::verdictColor(static_cast<int>(lband)));
            ImGui::ProgressBar(static_cast<float>(lg.substantialSimilarity / 100.0), ImVec2(-1, 0),
                               (f0(lg.substantialSimilarity) + " / 100").c_str());
            ImGui::PopStyleColor();
            ImGui::TextColored(theme::verdictColor(static_cast<int>(lband)), "%s",
                               lg.classification.c_str());
            ImGui::EndTable();
        }

        ImGui::Spacing();
        theme::sectionHeader("PREDICTED BYPRODUCTS / INTERACTIONS");
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
}

}  // namespace stimlab::panels
