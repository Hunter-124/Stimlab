#include "ui/Panels.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <imgui.h>
#include <implot.h>

#include "assay/Design.h"
#include "bio/Cartoon.h"
#include "bio/CifReader.h"
#include "bio/Connectivity.h"
#include "bio/Imgt.h"
#include "bio/Liabilities.h"
#include "bio/NucSeq.h"
#include "bio/PdbReader.h"
#include "bio/Structure.h"
#include "chem/Canonical.h"
#include "chem/Formula.h"
#include "chem/Rings.h"
#include "chem/Aromaticity.h"
#include "chem/Solubility.h"
#include "chem/Descriptors.h"
#include "chem/Perceive.h"
#include "chem/Embed3D.h"
#include "chem/Smiles.h"
#include "core/AppPaths.h"
#include "core/Manifest.h"
#include "modules/IonizationModule.h"
#include "modules/MechanismModule.h"
#include "modules/Metabolites.h"
#include "modules/NucleicModule.h"
#include "modules/docking/Presets.h"
#include "modules/docking/Provisioning.h"
#include "modules/docking/ReceptorPrep.h"
#include "render/MolViewport.h"
#include "ui/AppShell.h"
#include "ui/Theme.h"
#include "workflow/Dag.h"

namespace biocad::panels {
namespace {

std::string f2(double v) { char b[40]; std::snprintf(b, sizeof b, "%.2f", v); return b; }
std::string f0(double v) { char b[40]; std::snprintf(b, sizeof b, "%.0f", v); return b; }

// A one-line rendering of a Quantity for a stat card or a table cell, where
// drawQuantity's full colour-plus-provenance row does not fit. "not computed" is
// spelled out rather than left blank: an empty cell reads as zero.
std::string quantityShort(const Quantity& q) {
    if (q.provenance == Provenance::NotComputed) return "not computed";
    std::string s = f2(q.value);
    if (q.error > 0.0) s += " +/- " + f2(q.error);
    if (!q.unit.empty()) s += " " + q.unit;
    return s;
}

void verdictText(Verdict v) {
    ImGui::TextColored(theme::verdictColor(static_cast<int>(v)), "%s", verdictLabel(v));
}

}  // namespace

void drawQuantity(const char* label, const Quantity& q) {
    ImGui::PushStyleColor(ImGuiCol_Text, ImGui::ColorConvertU32ToFloat4(theme::kTextDim));
    ImGui::TextUnformatted(label);
    ImGui::PopStyleColor();
    ImGui::SameLine();

    const ImVec4 col = theme::provenanceColor(q.provenance);
    if (q.provenance == Provenance::NotComputed) {
        ImGui::TextColored(col, "not computed%s%s", q.source.empty() ? "" : " - needs ",
                           q.source.c_str());
        return;
    }

    std::string text = f2(q.value);
    if (q.error > 0.0) text += " +/- " + f2(q.error);
    if (!q.unit.empty()) text += " " + q.unit;
    text += "   (" + std::string(provenanceLabel(q.provenance));
    if (!q.source.empty()) text += " - " + q.source;
    text += ")";
    ImGui::TextColored(col, "%s", text.c_str());
}

namespace {

void statCard(const char* title, const std::string& value, const char* sub, float w = 168.0f,
              float h = 84.0f, float valueScale = 1.6f) {
    // valueScale parameter kept for call-site compat; ignored — theme fonts govern size.
    (void)valueScale;
    const bool visible = theme::beginCard(title, ImVec2(w, h > 84.0f ? h : 104.0f));
    if (!visible) {
        theme::endCard();
        return;
    }
    theme::pushSmallStrong();
    ImGui::PushStyleColor(ImGuiCol_Text, ImGui::ColorConvertU32ToFloat4(theme::kTextDim));
    ImGui::TextUnformatted(title);
    ImGui::PopStyleColor();
    theme::popFont();
    theme::pushValue();
    ImGui::PushStyleColor(ImGuiCol_Text, ImGui::ColorConvertU32ToFloat4(theme::kTextHi));
    ImGui::TextUnformatted(value.c_str());
    ImGui::PopStyleColor();
    theme::popFont();
    if (sub) {
        ImGui::PushStyleColor(ImGuiCol_Text, ImGui::ColorConvertU32ToFloat4(theme::kTextDim));
        ImGui::TextUnformatted(sub);
        ImGui::PopStyleColor();
    }
    theme::endCard();
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

// Host one already-built MolScene in an interactive off-screen viewport. Extracted from
// molViewer3D so the protein cartoon (which builds its scene from a bio::Structure, not a
// chem::Conformer) reuses exactly the same camera, hover and wheel-capture behaviour instead of
// growing a second copy of it.
bool sceneViewport(AppShell& shell, const render::MolScene& scene, ViewerUiState& ui,
                   const std::string& key, float height) {
    render::MolViewport* vp = shell.viewer();
    if (!vp) {
        ImGui::TextDisabled("3D viewer unavailable (no DirectX device).");
        return false;
    }
    // Auto-fit only when the framed content changes, so a user's orbit survives every frame.
    if (key != ui.lastKey) {
        vp->resetCamera(scene);
        ui.lastKey = key;
    }

    const float avail = ImGui::GetContentRegionAvail().x;
    const int w = std::max(64, static_cast<int>(avail));
    const int h = std::max(64, static_cast<int>(height));

    if (scene.atoms.empty() && scene.mesh.empty()) {
        ImGui::Dummy(ImVec2(static_cast<float>(w), static_cast<float>(h)));
        ImGui::TextDisabled("(nothing to display)");
        return true;
    }

    void* srv = vp->draw(scene, w, h);
    if (!srv) {
        ImGui::TextDisabled("3D render failed.");
        return false;
    }

    // Host the rendered image in its own child window flagged NoScrollWithMouse +
    // NoScrollbar so ImGui absorbs the wheel here (zoom) instead of scrolling the panel.
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
    ImGui::BeginChild(("##vp3d_" + key).c_str(),
                      ImVec2(static_cast<float>(w), static_cast<float>(h)),
                      ImGuiChildFlags_None,
                      ImGuiWindowFlags_NoScrollWithMouse | ImGuiWindowFlags_NoScrollbar);
    ImGui::PopStyleVar();
    ImGui::Image(reinterpret_cast<ImTextureID>(srv),
                 ImVec2(static_cast<float>(w), static_cast<float>(h)));
    const bool hovered = ImGui::IsItemHovered();
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
    return true;
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

    // Auto-fit when the molecule (or a toggle that changes framing) changes.
    const std::string key = cacheKey + (ui.showH ? "|H" : "") + (ui.spacefill ? "|S" : "") +
                            (overlayReceptor ? "|R" : "");
    if (!sceneViewport(shell, scene, ui, key, height)) return false;
    if (conf.empty() || scene.atoms.empty()) return true;

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
    line("# BioCAD report - " + m.name);
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
        line("## Stability: " + f0(r.overallScore) + "/100 (rank ordering, no units)");
        line("- Shelf life: " + quantityShort(r.shelfLife) +
             (r.shelfLife.provenance == Provenance::NotComputed
                  ? " - needs " + r.shelfLife.source
                  : std::string(" (") + provenanceLabel(r.shelfLife.provenance) + ")"));
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

    // ---- Hero card ----------------------------------------------------------
    const bool heroVisible = theme::beginCard("hero", ImVec2(-1.0f, 132.0f));
    if (heroVisible) {
        // Compound name
        theme::pushH2();
        ImGui::PushStyleColor(ImGuiCol_Text, ImGui::ColorConvertU32ToFloat4(theme::kTextHi));
        ImGui::TextUnformatted(m.name.c_str());
        ImGui::PopStyleColor();
        theme::popFont();
        ImGui::SameLine(0, 10);
        theme::badge(m.drugClass.c_str());
        ImGui::SameLine(0, 6);
        theme::badge(m.legalStatus.c_str(), theme::kTextDim, theme::kSurfaceHi);
        // Formula + truncated SMILES
        std::string smilesTrunc = m.smiles;
        if (smilesTrunc.size() > 60) smilesTrunc = smilesTrunc.substr(0, 60) + "\xe2\x80\xa6";  // UTF-8 ellipsis
        ImGui::PushStyleColor(ImGuiCol_Text, ImGui::ColorConvertU32ToFloat4(theme::kTextDim));
        ImGui::Text("%s   %s", m.formula.c_str(), smilesTrunc.c_str());
        ImGui::PopStyleColor();
        // Inline physchem chips
        ImGui::PushStyleColor(ImGuiCol_Text, ImGui::ColorConvertU32ToFloat4(theme::kTextDim));
        ImGui::Text("MW %.1f   logP %.2f   TPSA %.0f   HBD %d   HBA %d",
                    m.molWeight, m.logP, m.tpsa, m.hbd, m.hba);
        ImGui::PopStyleColor();
    }
    theme::endCard();
    ImGui::Spacing();

    // ---- Responsive metric grid ---------------------------------------------
    if (s.stability && s.absorption && s.admet) {
        const auto stab = s.stability->analyze(m);
        const auto abs  = s.absorption->predict(m);
        const auto adm  = s.admet->screen(m);

        const float spacing = 8.0f;
        const float avail   = ImGui::GetContentRegionAvail().x;
        const int cols      = std::max(1, static_cast<int>(avail / (215.0f + spacing)));
        const float cw      = (avail - spacing * (cols - 1)) / static_cast<float>(cols);

        // Helper to lay out a tile, advancing row as needed.
        int tileIdx = 0;
        auto nextTile = [&]() {
            const int col = tileIdx % cols;
            if (col != 0) ImGui::SameLine(0, spacing);
            ++tileIdx;
        };

        // Tile 1 — Stability
        nextTile();
        theme::metricCard("STABILITY", (f0(stab.overallScore) + "/100").c_str(),
                          quantityShort(stab.shelfLife).c_str(), theme::kTextHi, cw);

        // Tile 2 — Oral F
        nextTile();
        theme::metricCard("ORAL F", (f0(abs.bioavailabilityPct) + "%").c_str(),
                          abs.cnsPenetrant ? "BBB-permeant" : "low BBB partition",
                          theme::kTextHi, cw);

        // Tile 3 — HIA
        nextTile();
        theme::metricCard("HIA", (f0(abs.hiaPct) + "%").c_str(), "intestinal abs.",
                          theme::kTextHi, cw);

        // Tile 4 — ADMET
        {
            nextTile();
            const std::string epCount = std::to_string(adm.endpoints.size()) + " endpoints";
            const ImVec4 admetCol = theme::verdictColor(static_cast<int>(adm.overall));
            const ImU32 admetColU32 = ImGui::ColorConvertFloat4ToU32(admetCol);
            theme::metricCard("ADMET", verdictLabel(adm.overall), epCount.c_str(),
                              admetColU32, cw);
        }

        // Tile 5 — Library
        nextTile();
        theme::metricCard("LIBRARY",
                          std::to_string(s.library ? s.library->count() : 0).c_str(),
                          "compounds", theme::kTextHi, cw);

        // Tile 6 — Runs
        nextTile();
        theme::metricCard("RUNS",
                          std::to_string(s.runs ? s.runs->recent().size() : 0).c_str(),
                          "this session", theme::kTextHi, cw);

        // ---- Summaries ------------------------------------------------------
        ImGui::Spacing();
        theme::sectionHeader("Snapshot");

        // Two side-by-side summary cards, each half width.
        const float halfW = (avail - spacing) * 0.5f;
        const bool stabilitySummaryVisible =
            theme::beginCard("##stabsum", ImVec2(halfW, 0.0f), true);
        if (stabilitySummaryVisible) {
            ImGui::PushStyleColor(ImGuiCol_Text, ImGui::ColorConvertU32ToFloat4(theme::kText));
            ImGui::TextWrapped("%s", stab.summary.c_str());
            ImGui::PopStyleColor();
        }
        theme::endCard();
        ImGui::SameLine(0, spacing);
        const bool absorptionSummaryVisible =
            theme::beginCard("##abssum", ImVec2(halfW, 0.0f), true);
        if (absorptionSummaryVisible) {
            ImGui::PushStyleColor(ImGuiCol_Text, ImGui::ColorConvertU32ToFloat4(theme::kText));
            ImGui::TextWrapped("%s", abs.summary.c_str());
            ImGui::PopStyleColor();
        }
        theme::endCard();
    } else {
        // Services not yet ready — show library/runs tiles only.
        const float spacing = 8.0f;
        const float avail   = ImGui::GetContentRegionAvail().x;
        const int cols      = std::max(1, static_cast<int>(avail / (215.0f + spacing)));
        const float cw      = (avail - spacing * (cols - 1)) / static_cast<float>(cols);
        theme::metricCard("LIBRARY",
                          std::to_string(s.library ? s.library->count() : 0).c_str(),
                          "compounds", theme::kTextHi, cw);
        ImGui::SameLine(0, spacing);
        theme::metricCard("RUNS",
                          std::to_string(s.runs ? s.runs->recent().size() : 0).c_str(),
                          "this session", theme::kTextHi, cw);
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
        if (auto parsed = chem::parsePerceived(m.smiles)) conf = chem::embed3D(*parsed);
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
        theme::sectionHeader("Physicochemical");
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
        theme::sectionHeader("3D Structure");
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
        "Enter any SMILES to analyze an arbitrary structure with the full BioCAD engine "
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
        auto parsed = chem::parsePerceived(smilesBuf);
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
        theme::sectionHeader("3D Structure");
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
    // Not a stat card: a shelf life without its provenance and its missing
    // prerequisite is the number this phase deleted.
    drawQuantity("Shelf life", r.shelfLife);
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
    theme::sectionHeader("Likely Degradants");
    for (const auto& d : r.degradants)
        ImGui::BulletText("%s  -  %s (%s)", d.name.c_str(), d.pathway.c_str(), d.note.c_str());
}

// ------------------------------------------------------------------ Absorption
void absorption(AppShell& shell) {
    const Molecule m = shell.currentMolecule();
    if (!shell.services().absorption) return;
    const auto r = shell.services().absorption->predict(m);

    statCard("ORAL F", f0(r.bioavailabilityPct) + "%", "assumed-CLint, rank order", 190.0f);
    ImGui::SameLine();
    statCard("HIA", f0(r.hiaPct) + "%", "intestinal absorption", 170.0f);
    ImGui::SameLine();
    statCard("logBB", f2(r.logBB), r.cnsPenetrant ? "BBB-permeant" : "low BBB partition", 150.0f);
    ImGui::SameLine();
    statCard("P-gp", r.pgpSubstrate ? "substrate" : "no", "efflux", 150.0f);
    ImGui::Spacing();

    // Hepatic availability rests on an ASSUMED fu.CLint, so it is a rank-ordering
    // score with no unit - saying "%" here would be the exact dishonesty the
    // provenance rule exists to prevent.
    drawQuantity("Hepatic availability (rank order)",
                 makeQuantity(r.bioavailabilityPct / 100.0, "", 0.0, Provenance::Heuristic,
                              "well-stirred model, Q_H = 90 L/h, assumed fu.CLint"));
    drawQuantity("Absorbed fraction (rank order)",
                 makeQuantity(r.hiaPct / 100.0, "", 0.0, Provenance::Heuristic,
                              "Veber/Egan descriptor model"));
    // Aqueous solubility comes back NotComputed without a melting point, and that is
    // exactly what the reader must see instead of a fabricated logS.
    drawQuantity("Aqueous solubility (log10 mol/L)", r.logS);
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

// ------------------------------------------------------------ Structural Alerts
// The banner is fixed and unconditional. Everything in this panel is a LIABILITY
// FLAG: a substructure the literature has associated with reactive-metabolite
// formation. Nothing here is a toxicity verdict, and the "no alerts matched" case
// is stated as the non-claim it is rather than as a clean bill of health.
void alerts(AppShell& shell) {
    const Molecule m = shell.currentMolecule();
    if (!shell.services().alerts) return;
    const AlertReport r = shell.services().alerts->screen(m);

    ImGui::TextColored(theme::verdictColor(static_cast<int>(Verdict::Warn)),
                       "LIABILITY FLAGS, NOT A TOXICITY VERDICT");
    ImGui::TextWrapped(
        "A flag means the matched substructure has been ASSOCIATED with reactive-metabolite "
        "formation in the medicinal-chemistry literature. It does not say this compound is "
        "toxic, and it is not a prediction that bioactivation occurs: that depends on the "
        "enzymes present, the competing clearance routes, the dose and the detoxication "
        "capacity, none of which a substructure knows. Widely used marketed drugs match "
        "several of these alerts.");
    ImGui::Separator();

    ImGui::TextWrapped("%s", r.summary.c_str());

    if (r.flags.empty()) {
        ImGui::Spacing();
        ImGui::TextWrapped(
            "No alert matched. Read that literally: it means none of the motifs in this short, "
            "in-house pack is present. It is NOT a safety claim, and it says nothing about "
            "motifs the pack does not list.");
        return;
    }

    ImGui::Spacing();
    if (ImGui::BeginTable("alerts", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Flag", ImGuiTableColumnFlags_WidthFixed, 250.0f);
        ImGui::TableSetupColumn("Level", ImGuiTableColumnFlags_WidthFixed, 60.0f);
        ImGui::TableSetupColumn("Atoms", ImGuiTableColumnFlags_WidthFixed, 50.0f);
        ImGui::TableSetupColumn("Metabolic route it is associated with, and the source");
        ImGui::TableHeadersRow();
        for (const auto& f : r.flags) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0); ImGui::TextWrapped("%s", f.label.c_str());
            ImGui::TableSetColumnIndex(1); verdictText(f.severity);
            ImGui::TableSetColumnIndex(2); ImGui::Text("%d", f.atomCount);
            ImGui::TableSetColumnIndex(3);
            ImGui::TextWrapped("%s", f.mechanism.c_str());
            // The citation travels with the flag: an alert whose source is hidden is
            // indistinguishable from an invented rule.
            ImGui::TextDisabled("%s", f.citation.c_str());
        }
        ImGui::EndTable();
    }
}

// ----------------------------------------------------------------------- PK / PD
namespace {

// An assumed default is a CONSTRUCTED artefact, not a prediction: F, ka and fu have
// no credible structure-only predictor, so they are Model-tier and the source string
// says out loud that nothing computed them.
Quantity pkAssumed(double v, const char* unit, const char* why) {
    return makeQuantity(v, unit, 0.0, Provenance::Model,
                        std::string("assumed default - ") + why);
}

PkModelSpec defaultPkSpec() {
    PkModelSpec s;
    s.model              = PkModel::OralOneCompartment;
    s.bioavailability    = pkAssumed(0.80, "", "F is not predictable from structure");
    s.absorptionRate     = pkAssumed(1.20, "1/h", "ka is not predictable from structure");
    s.clearance          = pkAssumed(5.00, "L/h", "enter a measured CL to leave this tier");
    s.volume             = pkAssumed(50.0, "L", "enter a measured V to leave this tier");
    s.volumePeripheral   = pkAssumed(30.0, "L", "two-compartment only");
    s.intercompartmental = pkAssumed(4.00, "L/h", "two-compartment only");
    s.unboundFraction    = pkAssumed(0.50, "", "fu is not predictable from structure");
    s.vmax               = pkAssumed(0.00, "mg/h", "0 disables Michaelis-Menten elimination");
    s.km                 = pkAssumed(1.00, "mg/L", "used only when Vmax > 0");
    return s;
}

// One editable parameter row. Editing a value makes it user-entered, so the row
// colour changes from "model" to "measured" the moment a real number is supplied.
void pkParamRow(const char* id, const char* label, Quantity& q) {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::TextColored(theme::provenanceColor(q.provenance), "%s", label);
    ImGui::TableSetColumnIndex(1);
    ImGui::SetNextItemWidth(-1.0f);
    if (ImGui::InputDouble(id, &q.value, 0.0, 0.0, "%.4f")) {
        q.provenance = Provenance::Measured;
        q.source     = "user-entered";
    }
    ImGui::TableSetColumnIndex(2);
    ImGui::TextUnformatted(q.unit.empty() ? "-" : q.unit.c_str());
    ImGui::TableSetColumnIndex(3);
    ImGui::TextColored(theme::provenanceColor(q.provenance), "%s - %s",
                       provenanceLabel(q.provenance), q.source.c_str());
}

}  // namespace

void pkpd(AppShell& shell) {
    Services& s = shell.services();
    if (!s.pharmacodynamics) return;

    // Panel state is intentionally static: the regimen and the parameter set are the
    // user's scenario, not a property of the selected molecule.
    static PkModelSpec spec  = defaultPkSpec();
    static std::vector<DoseEvent> doses = {{0.0, 100.0, 0.0}, {12.0, 100.0, 0.0}};
    static Quantity    kd    = pkAssumed(0.05, "mg/L", "supply a measured Kd for a real occupancy");

    // Fixed banner - always the first thing drawn, never behind a scroll or a tab.
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic),
                       "EXPOSURE SCENARIO ONLY. This panel integrates a concentration-time and a "
                       "target-occupancy curve under the assumptions listed under the plot. It is "
                       "not, and will never be, a dose recommendation or a regimen suggestion.");
    ImGui::PopTextWrapPos();
    ImGui::Spacing();

    const char* kModelNames[] = {"IV bolus", "IV infusion", "Oral 1-compartment",
                                 "Oral 2-compartment"};
    int modelIdx = static_cast<int>(spec.model);
    ImGui::SetNextItemWidth(220.0f);
    if (ImGui::Combo("Structural model", &modelIdx, kModelNames, 4))
        spec.model = static_cast<PkModel>(modelIdx);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(110.0f);
    ImGui::InputDouble("Horizon (h)", &spec.horizonH, 1.0, 6.0, "%.1f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(110.0f);
    ImGui::InputDouble("Step (h)", &spec.stepH, 0.005, 0.05, "%.3f");
    ImGui::Spacing();

    DoseRegimen regimen;
    regimen.doses = doses;
    const PkProfile      prof = s.pharmacodynamics->simulate(spec, regimen);
    const OccupancyCurve occ  = s.pharmacodynamics->occupancy(prof, kd);

    statCard("Cmax", f2(prof.cmax.value), "mg/L (simulated)", 170.0f);
    ImGui::SameLine();
    statCard("Tmax", f2(prof.tmax.value), "h", 150.0f);
    ImGui::SameLine();
    statCard("AUC", f2(prof.auc.value), "mg*h/L over horizon", 190.0f);
    ImGui::SameLine();
    statCard("t 1/2", f2(prof.halfLife.value), "h", 150.0f);
    ImGui::SameLine();
    statCard("PEAK OCC.", f2(occ.peakOccupancy.value * 100.0) + "%", "fraction of target bound",
             200.0f);
    ImGui::Spacing();

    if (prof.flipFlop) {
        ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic),
                           "Flip-flop kinetics: ka < ke, so the terminal slope reflects "
                           "absorption, not elimination.");
    }

    const int n = static_cast<int>(prof.timeH.size());
    if (n > 1 && ImPlot::BeginPlot("##pkpd-conc", ImVec2(-1, 220))) {
        ImPlot::SetupAxes("time (h)", "concentration (mg/L)");
        ImPlot::PlotLine("total", prof.timeH.data(), prof.concentrationMgPerL.data(), n);
        if (prof.unboundMgPerL.size() == prof.timeH.size())
            ImPlot::PlotLine("unbound", prof.timeH.data(), prof.unboundMgPerL.data(), n);
        ImPlot::EndPlot();
    }

    const int on = static_cast<int>(occ.timeH.size());
    if (on > 1 && ImPlot::BeginPlot("##pkpd-occ", ImVec2(-1, 200))) {
        ImPlot::SetupAxes("time (h)", "fractional occupancy");
        ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, 1.0, ImPlotCond_Always);
        ImPlot::PlotLine("occupancy", occ.timeH.data(), occ.occupancy.data(), on);
        ImPlot::EndPlot();
    }

    // The assumptions belong UNDER the plot, verbatim: a curve whose assumptions are
    // not visible is a curve that misleads.
    theme::sectionHeader("ASSUMPTIONS THIS CURVE RESTS ON");
    if (prof.assumptions.empty()) {
        ImGui::TextDisabled("(the module reported no assumptions)");
    } else {
        for (const auto& a : prof.assumptions) ImGui::BulletText("%s", a.c_str());
    }
    if (!prof.note.empty()) ImGui::TextWrapped("%s", prof.note.c_str());
    if (!occ.note.empty()) ImGui::TextWrapped("%s", occ.note.c_str());
    ImGui::Spacing();

    theme::sectionHeader("DERIVED EXPOSURE QUANTITIES");
    drawQuantity("Cmax", prof.cmax);
    drawQuantity("Tmax", prof.tmax);
    drawQuantity("AUC", prof.auc);
    drawQuantity("Half-life", prof.halfLife);
    drawQuantity("Accumulation (Rac)", prof.accumulation);
    drawQuantity("Peak occupancy", occ.peakOccupancy);
    drawQuantity("Time above 50% occupancy", occ.timeAbove50Pct);
    ImGui::Spacing();

    theme::sectionHeader("DOSE EVENTS");
    if (ImGui::BeginTable("pkpd-doses", 4,
                          ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Time (h)", ImGuiTableColumnFlags_WidthFixed, 150.0f);
        ImGui::TableSetupColumn("Amount (mg)", ImGuiTableColumnFlags_WidthFixed, 150.0f);
        ImGui::TableSetupColumn("Duration (h)", ImGuiTableColumnFlags_WidthFixed, 150.0f);
        ImGui::TableSetupColumn("");
        ImGui::TableHeadersRow();
        int removeIdx = -1;
        for (int i = 0; i < static_cast<int>(doses.size()); ++i) {
            ImGui::PushID(i);
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::SetNextItemWidth(-1.0f);
            ImGui::InputDouble("##t", &doses[i].timeH, 0.0, 0.0, "%.2f");
            ImGui::TableSetColumnIndex(1);
            ImGui::SetNextItemWidth(-1.0f);
            ImGui::InputDouble("##amt", &doses[i].amountMg, 0.0, 0.0, "%.2f");
            ImGui::TableSetColumnIndex(2);
            ImGui::SetNextItemWidth(-1.0f);
            ImGui::InputDouble("##dur", &doses[i].durationH, 0.0, 0.0, "%.2f");
            ImGui::TableSetColumnIndex(3);
            if (ImGui::SmallButton("Remove")) removeIdx = i;
            ImGui::PopID();
        }
        ImGui::EndTable();
        if (removeIdx >= 0) doses.erase(doses.begin() + removeIdx);
    }
    if (ImGui::Button("Add dose event")) {
        const double last = doses.empty() ? 0.0 : doses.back().timeH;
        doses.push_back({last + 12.0, doses.empty() ? 100.0 : doses.back().amountMg, 0.0});
    }
    ImGui::Spacing();

    theme::sectionHeader("PARAMETERS (ROW COLOUR = PROVENANCE)");
    if (ImGui::BeginTable("pkpd-params", 4,
                          ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Parameter", ImGuiTableColumnFlags_WidthFixed, 220.0f);
        ImGui::TableSetupColumn("Value", ImGuiTableColumnFlags_WidthFixed, 140.0f);
        ImGui::TableSetupColumn("Unit", ImGuiTableColumnFlags_WidthFixed, 70.0f);
        ImGui::TableSetupColumn("Provenance");
        ImGui::TableHeadersRow();
        pkParamRow("##F",  "Bioavailability F",      spec.bioavailability);
        pkParamRow("##ka", "Absorption rate ka",     spec.absorptionRate);
        pkParamRow("##CL", "Clearance CL",           spec.clearance);
        pkParamRow("##V",  "Central volume V",       spec.volume);
        pkParamRow("##V2", "Peripheral volume V2",   spec.volumePeripheral);
        pkParamRow("##Q",  "Intercompartmental Q",   spec.intercompartmental);
        pkParamRow("##fu", "Unbound fraction fu",    spec.unboundFraction);
        pkParamRow("##Vm", "Vmax (0 = linear)",      spec.vmax);
        pkParamRow("##Km", "Km",                     spec.km);
        pkParamRow("##Kd", "Target Kd (occupancy)",  kd);
        ImGui::EndTable();
    }
}

// ------------------------------------------------- Population PK (Phase 13)
namespace {

// The population panel does NOT re-simulate every frame: 300 subjects over a 48 h
// horizon is tens of millions of RK4 steps, and running that at 60 Hz would make the
// UI unusable and the seed meaningless to the user (they could never read a band
// before it was replaced). The result is cached and recomputed only when the user
// presses Simulate.
struct PopPkState {
    PkModelSpec      spec = defaultPkSpec();
    std::vector<DoseEvent> doses = {{0.0, 100.0, 0.0}};
    VariabilitySpec  variability;
    PopulationProfile profile;
    Quantity         kd = pkAssumed(0.05, "mg/L", "supply a measured Kd for a real occupancy");
    bool             haveResult = false;

    PopPkState() {
        variability.betweenSubject = true;
        variability.parameters = {"CL", "V"};
        variability.omega = {0.09, 0.0, 0.0, 0.04};
        variability.parameterCovariance = {0.01, 0.0, 0.0, 0.01};
        variability.proportionalResidualCv = 0.15;
        variability.additiveResidualSd = 0.0;
        variability.seed = 20260101;
        variability.subjects = 200;
        variability.sampler = "latin-hypercube";
        spec.model = PkModel::IvBolus;
        spec.horizonH = 48.0;
        spec.stepH = 0.05;
    }
};

// Omega is edited as variances and a correlation rather than as raw covariances,
// because a user typing four numbers into a covariance matrix produces a
// non-positive-definite matrix within about three keystrokes, and the engine
// (correctly) refuses to repair one.
void omegaEditor(VariabilitySpec& v) {
    if (v.parameters.size() != 2 || v.omega.size() != 4) return;
    double sd1 = std::sqrt(std::max(v.omega[0], 0.0));
    double sd2 = std::sqrt(std::max(v.omega[3], 0.0));
    double corr = (sd1 > 0 && sd2 > 0) ? v.omega[1] / (sd1 * sd2) : 0.0;
    bool changed = false;
    ImGui::SetNextItemWidth(140.0f);
    changed |= ImGui::InputDouble(("omega SD " + v.parameters[0]).c_str(), &sd1, 0.01, 0.1,
                                  "%.3f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(140.0f);
    changed |= ImGui::InputDouble(("omega SD " + v.parameters[1]).c_str(), &sd2, 0.01, 0.1,
                                  "%.3f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(160.0f);
    const double corrLo = -0.95, corrHi = 0.95;
    changed |= ImGui::SliderScalar("eta correlation", ImGuiDataType_Double, &corr, &corrLo,
                                   &corrHi, "%.2f");
    if (changed) {
        sd1 = std::max(sd1, 0.0);
        sd2 = std::max(sd2, 0.0);
        corr = std::clamp(corr, -0.95, 0.95);
        v.omega = {sd1 * sd1, corr * sd1 * sd2, corr * sd1 * sd2, sd2 * sd2};
    }
    ImGui::TextDisabled("Omega is the covariance of the LOG-scale random effects: an SD of "
                        "0.30 is roughly a 30%% coefficient of variation on the parameter.");
}

}  // namespace

void popPk(AppShell& shell) {
    Services& s = shell.services();
    if (!s.populationPk) return;
    static PopPkState st;

    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic),
                       "EXPOSURE SCENARIO ONLY. The band below is the variability that was "
                       "ENTERED - Omega, the fit covariance and the residual error - pushed "
                       "through the stated PK model. It is not a prediction about an "
                       "individual, and this panel never emits a dose or a dose adjustment.");
    ImGui::PopTextWrapPos();
    ImGui::Spacing();

    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputInt("Subjects", &st.variability.subjects);
    st.variability.subjects = std::clamp(st.variability.subjects, 1, 5000);
    ImGui::SameLine();
    int seed = static_cast<int>(st.variability.seed);
    ImGui::SetNextItemWidth(160.0f);
    if (ImGui::InputInt("Seed", &seed))
        st.variability.seed = static_cast<std::uint64_t>(std::max(seed, 0));
    ImGui::SameLine();
    ImGui::SetNextItemWidth(110.0f);
    ImGui::InputDouble("Horizon (h)", &st.spec.horizonH, 1.0, 6.0, "%.1f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(110.0f);
    ImGui::InputDouble("Step (h)", &st.spec.stepH, 0.01, 0.05, "%.3f");

    ImGui::Checkbox("Between-subject variability", &st.variability.betweenSubject);
    ImGui::SameLine();
    ImGui::Checkbox("Parameter uncertainty", &st.variability.parameterUncertainty);
    ImGui::SameLine();
    ImGui::Checkbox("Residual error", &st.variability.residualError);
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputDouble("Proportional residual CV", &st.variability.proportionalResidualCv, 0.01,
                       0.05, "%.3f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputDouble("Additive residual SD (mg/L)", &st.variability.additiveResidualSd, 0.01,
                       0.05, "%.3f");
    omegaEditor(st.variability);
    ImGui::Spacing();

    if (ImGui::Button("Simulate population")) {
        DoseRegimen regimen;
        regimen.doses = st.doses;
        st.profile = s.populationPk->simulate(st.spec, regimen, st.variability);
        st.haveResult = true;
    }
    ImGui::SameLine();
    ImGui::TextDisabled("the simulation runs on demand, so the seed you see is the seed that "
                        "produced the band you see");
    if (!st.haveResult) {
        ImGui::Spacing();
        ImGui::TextDisabled("(no population simulated yet)");
        return;
    }

    const PopulationProfile& p = st.profile;
    const int n = static_cast<int>(p.times.size());
    if (n > 1 && ImPlot::BeginPlot("##poppk-band", ImVec2(-1, 260))) {
        ImPlot::SetupAxes("time (h)", "concentration (mg/L)");
        std::vector<double> lo(n), mid(n), hi(n);
        for (int i = 0; i < n; ++i) {
            lo[i] = p.bands[static_cast<std::size_t>(i)].p5;
            mid[i] = p.bands[static_cast<std::size_t>(i)].p50;
            hi[i] = p.bands[static_cast<std::size_t>(i)].p95;
        }
        // The faint individual trajectories go UNDER the band: they are what the band
        // is a summary of, and drawing them on top would imply they are the result.
        for (const auto& tr : p.sampleTrajectories) {
            if (static_cast<int>(tr.size()) != n) continue;
            // Unlabelled ("##") so fifty subjects do not produce fifty legend rows;
            // no per-item style call, because ImPlot's styling entry points differ
            // between the versions this tree may be built against and a legend the
            // reader can use matters more than the exact alpha.
            ImPlot::PlotLine("##subject", p.times.data(), tr.data(), n);
        }
        ImPlot::PlotShaded("5th-95th percentile", p.times.data(), lo.data(), hi.data(), n);
        ImPlot::PlotLine("median", p.times.data(), mid.data(), n);
        ImPlot::EndPlot();
    }
    ImGui::TextDisabled("%d of %d simulated subjects are drawn as faint lines; the percentiles "
                        "use every subject.",
                        static_cast<int>(p.sampleTrajectories.size()), p.spec.subjects);

    // Occupancy of the MEDIAN profile, not of a median occupancy: theta is nonlinear
    // in concentration, so the two differ and only the first is a curve of anything.
    if (s.pharmacodynamics && n > 1) {
        PkProfile medianProfile;
        medianProfile.timeH = p.times;
        medianProfile.concentrationMgPerL.resize(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i)
            medianProfile.concentrationMgPerL[static_cast<std::size_t>(i)] =
                p.bands[static_cast<std::size_t>(i)].p50;
        medianProfile.unboundMgPerL = medianProfile.concentrationMgPerL;
        const OccupancyCurve occ = s.pharmacodynamics->occupancy(medianProfile, st.kd);
        if (occ.timeH.size() > 1 && ImPlot::BeginPlot("##poppk-occ", ImVec2(-1, 180))) {
            ImPlot::SetupAxes("time (h)", "fractional occupancy of the median profile");
            ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, 1.0, ImPlotCond_Always);
            ImPlot::PlotLine("median occupancy", occ.timeH.data(), occ.occupancy.data(),
                             static_cast<int>(occ.timeH.size()));
            ImPlot::EndPlot();
        }
        if (!occ.note.empty()) ImGui::TextWrapped("%s", occ.note.c_str());
    }

    theme::sectionHeader("POPULATION SUMMARY");
    drawQuantity("Median AUC", p.medianAuc);
    drawQuantity("Median Cmax", p.medianCmax);
    drawQuantity("AUC CV across subjects", p.aucCv);
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::Model), "%s",
                       p.provenanceStatement.c_str());
    ImGui::PopTextWrapPos();

    // NCA of the median profile, so the same panel shows the model-free reading of
    // the curve it just simulated.
    theme::sectionHeader("NONCOMPARTMENTAL ANALYSIS OF THE MEDIAN PROFILE");
    ConcentrationSeries series;
    series.subjectId = "median";
    series.timeH = p.times;
    series.concentration.resize(p.bands.size());
    for (std::size_t i = 0; i < p.bands.size(); ++i) series.concentration[i] = p.bands[i].p50;
    double total = 0;
    for (const auto& d : st.doses) total += d.amountMg;
    series.dose = total;
    series.intravenous = st.spec.model == PkModel::IvBolus || st.spec.model == PkModel::IvInfusion;
    const NcaResult nca = s.populationPk->nca(series);
    if (ImGui::BeginTable("poppk-nca", 2, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Parameter", ImGuiTableColumnFlags_WidthFixed, 260.0f);
        ImGui::TableSetupColumn("Value");
        ImGui::TableHeadersRow();
        const std::pair<const char*, const Quantity*> rows[] = {
            {"Cmax", &nca.cmax},          {"Tmax", &nca.tmax},
            {"AUClast", &nca.aucLast},    {"AUCinf", &nca.aucInfinity},
            {"% extrapolated", &nca.percentExtrapolated},
            {"lambda_z", &nca.lambdaZ},   {"Half-life", &nca.halfLife},
            {"CL (or CL/F)", &nca.clearance}, {"Vz (or Vz/F)", &nca.volumeZ},
            {"Vss (IV only)", &nca.volumeSteadyState},
            {"AUMC", &nca.aumc},          {"MRT", &nca.meanResidenceTime},
            {"AUCtau", &nca.aucTau},      {"Cavg", &nca.cAverage},
            {"Swing", &nca.swing}};
        for (const auto& [label, q] : rows) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextUnformatted(label);
            ImGui::TableSetColumnIndex(1);
            drawQuantity("", *q);
        }
        ImGui::EndTable();
    }
    ImGui::Text("lambda_z used %d points, adjusted R-squared %.4f", nca.lambdaZPointCount,
                nca.adjustedRSquared);
    if (nca.extrapolationUnreliable)
        ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic),
                           "Extrapolation above 20%%: AUCinf and everything derived from it "
                           "are unreliable for this sampling window.");
    for (const auto& w : nca.warnings) ImGui::BulletText("%s", w.c_str());

    theme::sectionHeader("ASSUMPTIONS THIS BAND RESTS ON");
    for (const auto& a : p.assumptions) ImGui::BulletText("%s", a.c_str());
    for (const auto& w : p.warnings)
        ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic), "%s", w.c_str());

    theme::sectionHeader("PARAMETERS (ROW COLOUR = PROVENANCE)");
    if (ImGui::BeginTable("poppk-params", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Parameter", ImGuiTableColumnFlags_WidthFixed, 220.0f);
        ImGui::TableSetupColumn("Value", ImGuiTableColumnFlags_WidthFixed, 140.0f);
        ImGui::TableSetupColumn("Unit", ImGuiTableColumnFlags_WidthFixed, 70.0f);
        ImGui::TableSetupColumn("Provenance");
        ImGui::TableHeadersRow();
        pkParamRow("##pF", "Bioavailability F", st.spec.bioavailability);
        pkParamRow("##pka", "Absorption rate ka", st.spec.absorptionRate);
        pkParamRow("##pCL", "Clearance CL", st.spec.clearance);
        pkParamRow("##pV", "Central volume V", st.spec.volume);
        pkParamRow("##pfu", "Unbound fraction fu", st.spec.unboundFraction);
        pkParamRow("##pKd", "Target Kd (occupancy)", st.kd);
        ImGui::EndTable();
    }

    theme::sectionHeader("DOSE EVENTS");
    if (ImGui::BeginTable("poppk-doses", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Time (h)", ImGuiTableColumnFlags_WidthFixed, 150.0f);
        ImGui::TableSetupColumn("Amount (mg)", ImGuiTableColumnFlags_WidthFixed, 150.0f);
        ImGui::TableSetupColumn("Duration (h)");
        ImGui::TableHeadersRow();
        for (int i = 0; i < static_cast<int>(st.doses.size()); ++i) {
            ImGui::PushID(i);
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::SetNextItemWidth(-1.0f);
            ImGui::InputDouble("##t", &st.doses[i].timeH, 0.0, 0.0, "%.2f");
            ImGui::TableSetColumnIndex(1);
            ImGui::SetNextItemWidth(-1.0f);
            ImGui::InputDouble("##a", &st.doses[i].amountMg, 0.0, 0.0, "%.2f");
            ImGui::TableSetColumnIndex(2);
            ImGui::SetNextItemWidth(-1.0f);
            ImGui::InputDouble("##d", &st.doses[i].durationH, 0.0, 0.0, "%.2f");
            ImGui::PopID();
        }
        ImGui::EndTable();
    }
    if (ImGui::Button("Add dose event"))
        st.doses.push_back({st.doses.empty() ? 0.0 : st.doses.back().timeH + 12.0, 100.0, 0.0});
}

// -------------------------------------------- Interaction scenarios (Phase 13)
void interactionScenarios(AppShell& shell) {
    Services& s = shell.services();
    if (!s.populationPk) return;

    static PerpetratorSpec perp = [] {
        PerpetratorSpec p;
        p.label = "perpetrator";
        p.enzyme = "CYP3A4";
        p.ki = 0.5;
        p.kinact = -1;
        p.kI = -1;
        p.indEmax = -1;
        p.indEc50 = -1;
        p.indD = 1.0;
        p.unboundHepaticInletUM = 1.0;
        p.unboundSystemicUM = 0.5;
        p.enterocyteUM = 10.0;
        p.source = "enter the in vitro citation";
        return p;
    }();
    static VictimSpec victim = [] {
        VictimSpec v;
        v.label = "victim";
        // fm and Fg start ABSENT on purpose: the panel must show what a missing fm
        // does before it shows an AUCR.
        v.fractionMetabolizedByEnzyme = -1;
        v.intestinalAvailability = -1;
        v.fractionExcretedUnchanged = -1;
        v.source = "enter the clinical citation";
        return v;
    }();
    static double horizonH = 168.0;
    static double renalRatio = 1.0;
    static double hepaticRatio = 1.0;

    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic),
                       "AUC RATIOS ONLY. Every number here is an exposure ratio computed from "
                       "in vitro parameters you supplied. There is deliberately no risk score "
                       "and no red/green verdict: a 2-fold AUC ratio is not a hazard class, "
                       "and this panel never suggests a dose or a dose change.");
    ImGui::PopTextWrapPos();
    ImGui::Spacing();

    char enzyme[32];
    std::snprintf(enzyme, sizeof enzyme, "%s", perp.enzyme.c_str());
    ImGui::SetNextItemWidth(160.0f);
    if (ImGui::InputText("Enzyme", enzyme, sizeof enzyme)) perp.enzyme = enzyme;
    ImGui::SameLine();
    ImGui::TextDisabled("kdeg is looked up as <enzyme>_hepatic in assets/packs/physiology.json");

    theme::sectionHeader("PERPETRATOR IN VITRO PARAMETERS (-1 = ABSENT)");
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputDouble("Ki (uM)", &perp.ki, 0.0, 0.0, "%.4f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputDouble("kinact (1/h)", &perp.kinact, 0.0, 0.0, "%.4f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputDouble("KI (uM)", &perp.kI, 0.0, 0.0, "%.4f");
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputDouble("Induction Emax", &perp.indEmax, 0.0, 0.0, "%.4f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputDouble("Induction EC50 (uM)", &perp.indEc50, 0.0, 0.0, "%.4f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputDouble("Induction d", &perp.indD, 0.0, 0.0, "%.2f");
    ImGui::SetNextItemWidth(160.0f);
    ImGui::InputDouble("[I]h,u (uM)", &perp.unboundHepaticInletUM, 0.0, 0.0, "%.4f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(160.0f);
    ImGui::InputDouble("[I]sys,u (uM)", &perp.unboundSystemicUM, 0.0, 0.0, "%.4f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(160.0f);
    ImGui::InputDouble("[I]g (uM)", &perp.enterocyteUM, 0.0, 0.0, "%.4f");

    theme::sectionHeader("VICTIM DISPOSITION (-1 = ABSENT, WHICH IS NOT THE SAME AS 1)");
    ImGui::SetNextItemWidth(160.0f);
    ImGui::InputDouble("fm", &victim.fractionMetabolizedByEnzyme, 0.0, 0.0, "%.3f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(160.0f);
    ImGui::InputDouble("Fg", &victim.intestinalAvailability, 0.0, 0.0, "%.3f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(160.0f);
    ImGui::InputDouble("fe", &victim.fractionExcretedUnchanged, 0.0, 0.0, "%.3f");

    const InteractionReport rep = s.populationPk->interaction(perp, victim);

    theme::sectionHeader("FDA BASIC-MODEL SCREENING R-VALUES");
    drawQuantity("R1 (hepatic, reversible)", rep.r1);
    drawQuantity("R1,gut", rep.r1Gut);
    drawQuantity("R2 (time-dependent inactivation)", rep.r2);
    drawQuantity("R (induction)", rep.rInduction);

    theme::sectionHeader("MECHANISTIC STATIC AUC RATIO");
    drawQuantity("AUCR (hepatic x gut)", rep.aucRatio);
    drawQuantity("AUCR (hepatic only)", rep.aucRatioHepaticOnly);
    drawQuantity("Theoretical ceiling 1/(1-fm)", rep.theoreticalCeiling);
    ImGui::Text("Gut term included: %s", rep.gutIncluded ? "yes" : "no (Fg absent)");
    ImGui::TextWrapped("Dominant mechanism: %s", rep.dominantMechanism.c_str());
    for (const auto& w : rep.warnings)
        ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic), "%s", w.c_str());
    for (const auto& a : rep.assumptions) ImGui::BulletText("%s", a.c_str());

    theme::sectionHeader("DYNAMIC ENZYME ACTIVITY");
    ImGui::SetNextItemWidth(160.0f);
    ImGui::InputDouble("Horizon (h)", &horizonH, 12.0, 24.0, "%.0f");
    const EnzymeTimeCourse tc = s.populationPk->enzymeTimeCourse(perp, horizonH);
    const int tn = static_cast<int>(tc.timeH.size());
    if (tn > 1 && ImPlot::BeginPlot("##ddi-enzyme", ImVec2(-1, 220))) {
        ImPlot::SetupAxes("time (h)", "relative to uninhibited baseline");
        ImPlot::PlotLine("enzyme activity E", tc.timeH.data(), tc.relativeActivity.data(), tn);
        ImPlot::PlotLine("CLint(t)/CLint0", tc.timeH.data(), tc.clearanceRatio.data(), tn);
        ImPlot::EndPlot();
    }
    ImGui::Text("Dynamic steady-state activity %.9f versus the static model's %.9f "
                "(difference %.3g)",
                tc.steadyStateActivity, tc.staticModelActivity, tc.agreement);
    ImGui::Text("kdeg = %.5f /h from %s", tc.kdegUsed, tc.kdegSource.c_str());
    for (const auto& a : tc.assumptions) ImGui::BulletText("%s", a.c_str());

    theme::sectionHeader("ORGAN-IMPAIRMENT EXPOSURE SCENARIOS");
    ImGui::SetNextItemWidth(220.0f);
    const double ratioLo = 0.05, renalHi = 1.0, hepaticHi = 2.0;
    ImGui::SliderScalar("Renal function ratio", ImGuiDataType_Double, &renalRatio, &ratioLo,
                        &renalHi, "%.2f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(220.0f);
    ImGui::SliderScalar("Hepatic CLint ratio", ImGuiDataType_Double, &hepaticRatio, &ratioLo,
                        &hepaticHi, "%.2f");
    const ImpairmentScenario imp = s.populationPk->impairment(victim, renalRatio, hepaticRatio);
    drawQuantity("AUC impaired / AUC normal", imp.exposureRatio);
    for (const auto& a : imp.assumptions) ImGui::BulletText("%s", a.c_str());
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic), "%s",
                       imp.boundaryStatement.c_str());
    ImGui::PopTextWrapPos();
}

// ------------------------------------------------------- Sequence Compare
namespace {

// The app loads Segoe UI, which is proportional, so a plain TextUnformatted would
// stagger the three alignment rows against each other and make the midline point at
// the wrong column. Each glyph is therefore placed at a fixed advance, which is what
// makes the block readable as an alignment rather than as three sentences.
float monoCell() { return ImGui::CalcTextSize("M").x; }

void monoRow(const std::string& text, std::size_t from, std::size_t count, ImU32 col,
             float cell) {
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const ImVec2 origin = ImGui::GetCursorScreenPos();
    char buf[2] = {0, 0};
    for (std::size_t i = 0; i < count && from + i < text.size(); ++i) {
        buf[0] = text[from + i];
        dl->AddText(ImVec2(origin.x + static_cast<float>(i) * cell, origin.y), col, buf, buf + 1);
    }
    ImGui::Dummy(ImVec2(static_cast<float>(count) * cell, ImGui::GetTextLineHeight()));
}

void alignmentBlock(const SequenceAlignment& a) {
    const float cell = monoCell();
    if (cell <= 0.0f || a.aligned1.empty()) {
        ImGui::TextDisabled("(no alignment)");
        return;
    }
    const float avail = ImGui::GetContentRegionAvail().x - 60.0f;
    std::size_t perRow = static_cast<std::size_t>(avail / cell);
    if (perRow < 20) perRow = 20;

    for (std::size_t off = 0; off < a.aligned1.size(); off += perRow) {
        const std::size_t n = std::min(perRow, a.aligned1.size() - off);
        monoRow(a.aligned1, off, n, theme::kTextHi, cell);
        monoRow(a.midline, off, n, theme::kAccent2, cell);
        monoRow(a.aligned2, off, n, theme::kTextHi, cell);
        ImGui::PushStyleColor(ImGuiCol_Text, ImGui::ColorConvertU32ToFloat4(theme::kTextFaint));
        ImGui::Text("columns %zu-%zu", off + 1, off + n);
        ImGui::PopStyleColor();
        ImGui::Spacing();
    }
}

}  // namespace

void sequenceCompare(AppShell& shell) {
    Services& s = shell.services();
    if (!s.sequence) return;

    // The two sequences are the user's scenario, not a property of the selected
    // molecule, so they live in panel state like the PK regimen does.
    static char seqA[8192] =
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR";
    static char seqB[8192] =
        "MKTAYIAKQRQISFVKSHFSRQEEERLGLIEVQAAILSRVGDGTQDNLSGCEKAVQVKVKALPDAQFEVVHSLAKWKQ";
    static bool local = true;

    theme::sectionHeader("SEQUENCES (ONE-LETTER, GAPS AND WHITESPACE IGNORED)");
    ImGui::TextUnformatted("Sequence A");
    ImGui::InputTextMultiline("##seqA", seqA, sizeof(seqA), ImVec2(-1, 70));
    ImGui::TextUnformatted("Sequence B");
    ImGui::InputTextMultiline("##seqB", seqB, sizeof(seqB), ImVec2(-1, 70));

    int mode = local ? 1 : 0;
    ImGui::RadioButton("Global (Needleman-Wunsch)", &mode, 0);
    ImGui::SameLine();
    ImGui::RadioButton("Local (Smith-Waterman)", &mode, 1);
    local = (mode == 1);
    ImGui::Spacing();

    // Strip anything that is not a residue letter: a pasted FASTA header or a
    // wrapped line would otherwise be aligned as if it were sequence.
    auto clean = [](const char* raw) {
        std::string out;
        for (const char* p = raw; *p; ++p)
            if (std::isalpha(static_cast<unsigned char>(*p)))
                out += static_cast<char>(std::toupper(static_cast<unsigned char>(*p)));
        return out;
    };
    const std::string a = clean(seqA);
    const std::string b = clean(seqB);
    if (a.empty() || b.empty()) {
        ImGui::TextDisabled("Enter two sequences to align.");
        return;
    }

    const SequenceAlignment r = local ? s.sequence->alignLocal(a, b)
                                      : s.sequence->alignGlobal(a, b);

    statCard("IDENTITY", f2(r.identityPct.value) + "%", "identical / aligned columns", 210.0f);
    ImGui::SameLine();
    statCard("SIMILARITY", f2(r.similarityPct.value) + "%", "positive-scoring columns", 220.0f);
    ImGui::SameLine();
    statCard("SCORE", f0(r.score.value), "half-bits (BLOSUM62)", 190.0f);
    ImGui::SameLine();
    statCard("GAP OPENS", std::to_string(r.gapOpens), "gap runs, both sequences", 190.0f);
    ImGui::Spacing();

    theme::sectionHeader("ALIGNMENT");
    alignmentBlock(r);

    theme::sectionHeader("NUMBERS AND THEIR PROVENANCE");
    drawQuantity("Score", r.score);
    drawQuantity("Identity", r.identityPct);
    drawQuantity("Similarity", r.similarityPct);
    // The E-value ROW IS ABSENT for a global alignment, not blank: Karlin-Altschul
    // statistics have no global analogue, and an empty row invites the reader to
    // assume the number merely failed to compute this time.
    if (!local) drawQuantity("E-value", r.eValue);
    ImGui::TextColored(theme::provenanceColor(Provenance::Measured),
                       "Aligned columns: %d", r.alignedLength);
    if (!r.note.empty()) ImGui::TextWrapped("%s", r.note.c_str());
}

// ------------------------------------------------------- Known Metabolites
// Facts only. Every row here came from a citable source, which is the entire
// reason the table is coloured Measured; nothing on this surface is enumerated,
// and the note at the bottom says why in numbers.
void metabolites(AppShell& shell) {
    Services& s = shell.services();
    if (!s.metabolismFacts) return;
    const Molecule m = shell.currentMolecule();
    const MetabolismReport r = s.metabolismFacts->known(m);

    theme::sectionHeader("CURATED, CITED BIOTRANSFORMATIONS");
    ImGui::TextWrapped("%s", r.summary.c_str());
    ImGui::Spacing();

    if (!r.known.empty() &&
        ImGui::BeginTable("metabolite-facts", 5,
                          ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                              ImGuiTableFlags_Resizable)) {
        ImGui::TableSetupColumn("Metabolite", ImGuiTableColumnFlags_WidthFixed, 210.0f);
        ImGui::TableSetupColumn("Enzyme", ImGuiTableColumnFlags_WidthFixed, 190.0f);
        ImGui::TableSetupColumn("Reaction", ImGuiTableColumnFlags_WidthFixed, 230.0f);
        ImGui::TableSetupColumn("Significance");
        ImGui::TableSetupColumn("Citation", ImGuiTableColumnFlags_WidthFixed, 260.0f);
        ImGui::TableHeadersRow();
        const ImVec4 measured = theme::provenanceColor(Provenance::Measured);
        for (const auto& f : r.known) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextColored(measured, "%s", f.metaboliteName.c_str());
            // A structure is shown only where the pack authored one: an omitted
            // SMILES means the structure was not certain, and inventing one to fill
            // the cell would be exactly the fabrication this panel exists to avoid.
            if (!f.metaboliteSmiles.empty()) {
                ImGui::TextDisabled("%s", f.metaboliteSmiles.c_str());
            } else {
                ImGui::TextDisabled("(structure not authored)");
            }
            ImGui::TableSetColumnIndex(1);
            ImGui::TextColored(measured, "%s", f.enzyme.c_str());
            if (f.polymorphic) {
                ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic),
                                   "polymorphic enzyme");
                ImGui::TextWrapped("Genotype changes this route's flux, so the exposure differs "
                                   "between phenotypes rather than being one number.");
            }
            ImGui::TableSetColumnIndex(2);
            ImGui::TextWrapped("%s", f.reaction.c_str());
            ImGui::TableSetColumnIndex(3);
            ImGui::TextWrapped("%s", f.significance.c_str());
            ImGui::TableSetColumnIndex(4);
            ImGui::TextWrapped("%s", f.citation.c_str());
        }
        ImGui::EndTable();
        ImGui::Spacing();
        ImGui::TextColored(theme::provenanceColor(Provenance::Measured),
                           "Every row above is %s - a characterised transformation with the "
                           "reference attached.", provenanceLabel(Provenance::Measured));
    }

    // ALWAYS rendered, with or without facts: absence of a curated entry must never
    // read as absence of metabolism.
    ImGui::Spacing();
    theme::sectionHeader("COVERAGE");
    ImGui::TextWrapped("%s", r.coverageNote.c_str());

    ImGui::Spacing();
    theme::sectionHeader("WHY NOTHING HERE IS ENUMERATED");
    ImGui::TextWrapped("%s", metaboliteNoEnumerationNote());
}

// --------------------------------------------------- Ionization & Solubility
// Every curve on this surface rests on a pKa, a melting point, a Ksp or a rate
// constant, and every one of those is an INPUT. So the panel's first job is to
// say which inputs it had and which it did not: a missing prerequisite is named
// in the body, in the place the curve would have been, never tucked into a
// tooltip a reader can miss. Formula, exact mass and the isotope envelope are
// deliberately drawn first because they are the one part of this panel that is
// always available - they are arithmetic on measured isotope masses, not a model.
void ionization(AppShell& shell) {
    Services& s = shell.services();
    if (!s.ionization) return;
    const Molecule          m = shell.currentMolecule();
    const IonizationReport  r = s.ionization->analyze(m);

    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::Measured), "%s", ionizationInputNote());
    ImGui::PopTextWrapPos();
    ImGui::Spacing();

    // A pack that failed to parse is louder than a compound that is merely absent,
    // because the two call for opposite responses and look identical downstream.
    if (const auto* real = dynamic_cast<const RealIonization*>(s.ionization)) {
        const auto& pack = real->pack();
        if (!pack.errors.empty()) {
            theme::sectionHeader("IONIZATION PACK FAILED TO LOAD");
            for (const auto& e : pack.errors) ImGui::TextColored(theme::verdictColor(3), "%s", e.c_str());
            ImGui::Spacing();
        }
    }

    // ---------------------------------------------- formula and exact mass
    theme::sectionHeader("FORMULA AND EXACT MASS");
    statCard("FORMULA", r.mass.formula, "Hill order", 190.0f);
    ImGui::SameLine();
    statCard("MONOISOTOPIC", f2(r.mass.monoisotopic.value), "Da (most abundant isotopes)", 210.0f);
    ImGui::SameLine();
    statCard("AVERAGE", f2(r.mass.average.value), "Da (standard atomic weights)", 200.0f);
    ImGui::SameLine();
    statCard("RDBE", f2(r.mass.unsaturation), "rings + double bonds", 160.0f);
    ImGui::Spacing();
    drawQuantity("Monoisotopic mass", r.mass.monoisotopic);
    drawQuantity("Average mass", r.mass.average);
    // The m/z row is present even for a neutral, where it reads "not computed -
    // needs charge is zero": a blank row would invite the reader to assume the
    // number merely failed this time.
    drawQuantity("m/z", r.mass.mz);
    ImGui::TextColored(theme::provenanceColor(Provenance::Measured), "Electrons: %d",
                       r.mass.electrons);
    for (const auto& w : r.mass.warnings)
        ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic), "%s", w.c_str());
    ImGui::Spacing();

    // ------------------------------------------------- isotope envelope
    theme::sectionHeader("THEORETICAL ISOTOPE ENVELOPE");
    if (r.envelope.peaks.empty()) {
        ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed), "%s",
                           r.envelope.source.c_str());
    } else {
        std::vector<double> mz, rel;
        mz.reserve(r.envelope.peaks.size());
        rel.reserve(r.envelope.peaks.size());
        for (const auto& p : r.envelope.peaks) {
            mz.push_back(p.mass);
            rel.push_back(p.intensity);
        }
        if (ImPlot::BeginPlot("##ion-env", ImVec2(-1, 190),
                              ImPlotFlags_NoMouseText | ImPlotFlags_NoLegend)) {
            ImPlot::SetupAxes("m/z (Da)", "relative intensity");
            ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, 1.05, ImPlotCond_Always);
            ImPlot::PlotStems("isotopologue", mz.data(), rel.data(),
                              static_cast<int>(mz.size()));
            ImPlot::EndPlot();
        }
        ImGui::TextWrapped("Peaks below %.1e of the base peak were pruned. Isotope masses and "
                           "abundances: %s", r.envelope.prunedBelow, r.envelope.source.c_str());
    }
    ImGui::Spacing();

    // ------------------------------------------- microspecies vs pH
    theme::sectionHeader("MICROSPECIES DISTRIBUTION VS pH");
    const std::size_t nL = r.speciation.labels.size();
    const std::size_t nP = r.speciation.points.size();
    if (nP < 2 || nL == 0) {
        // The missing input, in the body, where the plot would have been.
        ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed),
                           "No microspecies plot: this needs %s.",
                           r.speciation.isoelectricPoint.source.empty()
                               ? "a cited pKa"
                               : r.speciation.isoelectricPoint.source.c_str());
    } else {
        std::vector<double> xs(nP);
        // Stacked areas need cumulative bounds, so cum[k] is the running sum below
        // species k and cum[k+1] the sum including it.
        std::vector<std::vector<double>> cum(nL + 1, std::vector<double>(nP, 0.0));
        for (std::size_t p = 0; p < nP; ++p) {
            xs[p] = r.speciation.points[p].pH;
            for (std::size_t k = 0; k < nL; ++k) {
                const auto& fr = r.speciation.points[p].microspeciesFractions;
                cum[k + 1][p] = cum[k][p] + (k < fr.size() ? fr[k] : 0.0);
            }
        }
        if (ImPlot::BeginPlot("##ion-spec", ImVec2(-1, 240))) {
            ImPlot::SetupAxes("pH", "fraction of compound");
            ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, 1.0, ImPlotCond_Always);
            for (std::size_t k = 0; k < nL; ++k)
                ImPlot::PlotShaded(r.speciation.labels[k].c_str(), xs.data(), cum[k].data(),
                                   cum[k + 1].data(), static_cast<int>(nP));
            ImPlot::EndPlot();
        }

        std::vector<double> charge(nP), logD(nP);
        for (std::size_t p = 0; p < nP; ++p) {
            charge[p] = r.speciation.points[p].netCharge;
            logD[p]   = r.speciation.points[p].logD;
        }
        if (ImPlot::BeginPlot("##ion-charge", ImVec2(-1, 190))) {
            ImPlot::SetupAxes("pH", "net charge (e) / log D");
            ImPlot::PlotLine("net charge", xs.data(), charge.data(), static_cast<int>(nP));
            // logD shares the axis on purpose: the pH at which the charge leaves
            // zero is exactly the pH at which logD departs from logP, and splitting
            // the two plots hides that they are the same event.
            ImPlot::PlotLine("log D", xs.data(), logD.data(), static_cast<int>(nP));
            ImPlot::EndPlot();
        }
    }
    drawQuantity("Input log P", r.speciation.logP);
    drawQuantity("log D at pH 7.4", r.speciation.logDAtPh74);
    drawQuantity("Isoelectric point", r.speciation.isoelectricPoint);
    for (const auto& w : r.speciation.warnings)
        ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic), "%s", w.c_str());
    ImGui::Spacing();

    // ------------------------------------------------------ pH-solubility
    theme::sectionHeader("pH-SOLUBILITY PROFILE");
    if (r.solubility.curve.size() < 2) {
        ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed),
                           "No solubility curve: this needs %s.",
                           r.solubility.intrinsic.source.empty()
                               ? "a melting point and a log P, or a measured intrinsic solubility"
                               : r.solubility.intrinsic.source.c_str());
    } else {
        std::vector<double> sx, sy, saltX, saltY;
        for (const auto& p : r.solubility.curve) {
            sx.push_back(p.pH);
            sy.push_back(p.logS);
            if (p.saltLimited) {
                saltX.push_back(p.pH);
                saltY.push_back(p.logS);
            }
        }
        // logS is plotted directly rather than putting S on a log axis: the value is
        // already log10 mol/L, and a decade grid on a linear axis of logs is the
        // same picture with no axis-scale API in the way.
        if (ImPlot::BeginPlot("##ion-sol", ImVec2(-1, 230))) {
            ImPlot::SetupAxes("pH", "log10 S (mol/L, total dissolved)");
            ImPlot::PlotLine("pH-dependent", sx.data(), sy.data(), static_cast<int>(sx.size()));
            if (!saltX.empty())
                ImPlot::PlotLine("salt-limited plateau", saltX.data(), saltY.data(),
                                 static_cast<int>(saltX.size()));
            ImPlot::EndPlot();
        }
        if (saltX.empty()) {
            ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed),
                               "No pHmax kink is drawn. The salt plateau needs a Ksp and a "
                               "counter-ion concentration, and a kink at an unknown pH would be "
                               "the most misleading feature on the plot.");
        }
    }
    drawQuantity("Intrinsic solubility S0", r.solubility.intrinsic);
    drawQuantity("pHmax (salt kink)", r.solubility.pHmax);
    drawQuantity("Solubility at pH 7.4", r.solubility.solubilityAtPh74);
    drawQuantity("Dose number Do", r.solubility.doseNumber);
    drawQuantity("Dissolution number Dn", r.solubility.dissolutionNumber);
    drawQuantity("Absorption number An", r.solubility.absorptionNumber);
    for (const auto& w : r.solubility.warnings) ImGui::BulletText("%s", w.c_str());
    ImGui::Spacing();

    // ----------------------------------------------------- buffer capacity
    theme::sectionHeader("BUFFER CAPACITY (VAN SLYKE)");
    if (r.buffer.curve.size() < 2) {
        ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed),
                           "No buffer-capacity curve: this needs %s.",
                           r.buffer.betaAtPh74.source.empty() ? "a cited pKa"
                                                              : r.buffer.betaAtPh74.source.c_str());
    } else {
        std::vector<double> bx, by;
        for (const auto& p : r.buffer.curve) {
            bx.push_back(p.pH);
            by.push_back(p.beta);
        }
        if (ImPlot::BeginPlot("##ion-buffer", ImVec2(-1, 200))) {
            ImPlot::SetupAxes("pH", "beta (mol/L per pH)");
            ImPlot::PlotLine("buffer value", bx.data(), by.data(), static_cast<int>(bx.size()));
            ImPlot::EndPlot();
        }
    }
    drawQuantity("beta at pH 7.4", r.buffer.betaAtPh74);
    drawQuantity("Maximum beta", r.buffer.maxCapacity);
    drawQuantity("pH of maximum beta", r.buffer.maxCapacityPh);
    ImGui::Spacing();

    // ------------------------------------------- dissolution / precipitation
    // Dissolution belongs to a FORMULATION, not to a molecule, so these four
    // numbers are the user's and analyze() refuses to invent them. Once they are
    // entered the time course is real physics on real inputs.
    theme::sectionHeader("DISSOLUTION AND pH-SHIFT PRECIPITATION");
    static double doseMg = 0.0, radiusUm = 0.0, densityGPerCm3 = 0.0, diffusivity = 0.0;
    static double filmUm = 30.0, kppt = 0.0;
    static bool   wantPrecipitation = false;
    ImGui::SetNextItemWidth(130.0f);
    ImGui::InputDouble("Dose (mg)", &doseMg, 0.0, 0.0, "%.1f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(130.0f);
    ImGui::InputDouble("Radius (um)", &radiusUm, 0.0, 0.0, "%.2f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(150.0f);
    ImGui::InputDouble("Density (g/cm3)", &densityGPerCm3, 0.0, 0.0, "%.3f");
    ImGui::SetNextItemWidth(160.0f);
    ImGui::InputDouble("Diffusivity (cm2/s)", &diffusivity, 0.0, 0.0, "%.2e");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(150.0f);
    ImGui::InputDouble("Film h (um)", &filmUm, 0.0, 0.0, "%.1f");
    ImGui::SameLine();
    ImGui::Checkbox("pH-shift precipitation", &wantPrecipitation);
    if (wantPrecipitation) {
        ImGui::SetNextItemWidth(160.0f);
        ImGui::InputDouble("kppt (1/s)", &kppt, 0.0, 0.0, "%.4e");
    }

    const bool haveSolubility = r.solubility.solubilityAtPh74.provenance != Provenance::NotComputed;
    const bool haveMass       = r.mass.average.provenance != Provenance::NotComputed;
    std::vector<std::string> missing;
    if (!(doseMg > 0.0)) missing.emplace_back("a dose in mg");
    if (!(radiusUm > 0.0)) missing.emplace_back("an initial particle radius");
    if (!(densityGPerCm3 > 0.0)) missing.emplace_back("a solid density");
    if (!(diffusivity > 0.0)) missing.emplace_back("a diffusivity");
    if (!(filmUm > 0.0)) missing.emplace_back("a diffusion-layer thickness");
    if (!haveSolubility) missing.emplace_back("a solubility to dissolve toward (needs a pKa)");
    if (!haveMass) missing.emplace_back("a molecular weight");
    if (wantPrecipitation && !(kppt > 0.0))
        missing.emplace_back("a precipitation rate constant kppt, which is never predicted");

    if (!missing.empty()) {
        ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed),
                           "No dissolution time course. Still needed:");
        for (const auto& mi : missing) ImGui::BulletText("%s", mi.c_str());
        drawQuantity("Time to 85% dissolved", r.dissolution.timeTo85Pct);
    } else {
        chem::DissolutionInput in;
        in.doseMg               = doseMg;
        in.molWeight            = r.mass.average.value;
        in.initialRadiusUm      = radiusUm;
        in.densityGPerCm3       = densityGPerCm3;
        in.diffusivityCm2PerS   = diffusivity;
        in.diffusionLayerUm     = filmUm;
        in.solubilityMolar      = r.solubility.solubilityAtPh74.value;
        in.precipitation        = wantPrecipitation;
        in.kpptPerS             = kppt;
        in.hasKppt              = wantPrecipitation && kppt > 0.0;
        in.precipSolubilityMolar = r.solubility.solubilityAtPh74.value;
        const DissolutionReport d = chem::dissolutionTimeCourse(in);

        std::vector<double> t, dissolved, solid, precipitated;
        for (const auto& p : d.points) {
            t.push_back(p.timeS / 60.0);   // minutes read better than seconds
            dissolved.push_back(p.dissolvedMolar);
            solid.push_back(p.solidMg);
            precipitated.push_back(p.precipitatedMg);
        }
        if (t.size() > 1 && ImPlot::BeginPlot("##ion-diss", ImVec2(-1, 220))) {
            ImPlot::SetupAxes("time (min)", "dissolved (mol/L)");
            ImPlot::PlotLine("dissolved", t.data(), dissolved.data(), static_cast<int>(t.size()));
            ImPlot::EndPlot();
        }
        if (t.size() > 1 && ImPlot::BeginPlot("##ion-diss-mass", ImVec2(-1, 200))) {
            ImPlot::SetupAxes("time (min)", "mass (mg)");
            ImPlot::PlotLine("undissolved solid", t.data(), solid.data(),
                             static_cast<int>(t.size()));
            ImPlot::PlotLine("precipitated", t.data(), precipitated.data(),
                             static_cast<int>(t.size()));
            ImPlot::EndPlot();
        }
        drawQuantity("Time to 85% dissolved", d.timeTo85Pct);
        // The mass-balance residual is shown, not asserted in private: a
        // dissolution model that loses mass is wrong in a way a single curve hides.
        ImGui::TextColored(theme::provenanceColor(Provenance::Measured),
                           "Worst solid+dissolved+precipitated mass imbalance: %.3e mg",
                           d.maxMassImbalance);
        for (const auto& a : d.assumptions) ImGui::BulletText("%s", a.c_str());
        for (const auto& w : d.warnings)
            ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic), "%s", w.c_str());
    }
    ImGui::Spacing();

    // ------------------------------------------------------- assumptions
    theme::sectionHeader("WHAT THESE CURVES REST ON");
    for (const auto& a : r.speciation.assumptions) ImGui::BulletText("%s", a.c_str());
    for (const auto& a : r.solubility.assumptions) ImGui::BulletText("%s", a.c_str());
    for (const auto& a : r.buffer.assumptions) ImGui::BulletText("%s", a.c_str());
    for (const auto& a : r.dissolution.assumptions) ImGui::BulletText("%s", a.c_str());
}

// ------------------------------------------------------ Assay Workbench
namespace {

// The default contents of the import box. It is a real long CSV, not a description
// of one: the fastest way to learn an import format is to see a working example you
// can edit in place.
const char* kAssaySampleCsv =
    "plate_id,well,role,sample_id,series_id,concentration,conc_unit,replicate,readout,readout_unit\n"
    "P1,A1,negative,ctrl,,0,M,0,10,RFU\n"
    "P1,B1,negative,ctrl,,0,M,1,11,RFU\n"
    "P1,C1,negative,ctrl,,0,M,2,9,RFU\n"
    "P1,D1,negative,ctrl,,0,M,3,10.5,RFU\n"
    "P1,A2,positive,ctrl,,0,M,0,100,RFU\n"
    "P1,B2,positive,ctrl,,0,M,1,102,RFU\n"
    "P1,C2,positive,ctrl,,0,M,2,98,RFU\n"
    "P1,D2,positive,ctrl,,0,M,3,101,RFU\n"
    "P1,A3,sample,cmpd1,s1,1e-5,M,0,11.0,RFU\n"
    "P1,A4,sample,cmpd1,s1,3.16e-6,M,0,12.7,RFU\n"
    "P1,A5,sample,cmpd1,s1,1e-6,M,0,18.4,RFU\n"
    "P1,A6,sample,cmpd1,s1,3.16e-7,M,0,35.5,RFU\n"
    "P1,A7,sample,cmpd1,s1,1e-7,M,0,64.5,RFU\n"
    "P1,A8,sample,cmpd1,s1,3.16e-8,M,0,81.6,RFU\n"
    "P1,A9,sample,cmpd1,s1,1e-8,M,0,87.3,RFU\n"
    "P1,A10,sample,cmpd1,s1,1e-9,M,0,89.6,RFU\n";

// A hollow marker is drawn with a transparent FILL and a visible outline. An
// excluded well is never removed from a plot: the reader has to be able to see the
// point that a rule threw away, and which rule threw it.
ImPlotSpec hollowMarkerSpec(ImPlotMarker marker, float size) {
    ImPlotSpec s;
    s.Marker = marker;
    s.MarkerSize = size;
    s.MarkerFillColor = ImVec4(0.0f, 0.0f, 0.0f, 0.0f);
    s.MarkerLineColor = ImVec4(1.0f, 0.75f, 0.25f, 1.0f);
    s.LineWeight = 2.0f;
    return s;
}

ImPlotSpec markerSpec(ImPlotMarker marker, float size) {
    ImPlotSpec s;
    s.Marker = marker;
    s.MarkerSize = size;
    return s;
}

std::string sci(double v) {
    char b[48];
    std::snprintf(b, sizeof b, "%.4g", v);
    return b;
}

const char* const kAssayModelNames[] = {"4PL", "5PL", "Michaelis-Menten", "Hill",
                                        "Substrate inhibition", "Morrison tight binding",
                                        "Langmuir 1:1 (SPR/BLI)", "Mass transport (SPR)",
                                        "Boltzmann melt (DSF)", "Two-state melt (DSF)",
                                        "Wiseman isotherm (ITC)"};

// The plate as a heat map. Excluded wells are drawn as an overlay of hollow squares
// on top of the map rather than being blanked out, for the same reason as above.
void assayHeatMap(const Plate& p) {
    if (p.rows <= 0 || p.columns <= 0) return;
    std::vector<double> grid(static_cast<std::size_t>(p.rows) * p.columns,
                             std::numeric_limits<double>::quiet_NaN());
    double lo = 0.0, hi = 0.0;
    bool first = true;
    std::vector<double> exX, exY;
    for (const auto& w : p.wells) {
        if (w.row < 0 || w.row >= p.rows || w.column < 0 || w.column >= p.columns) continue;
        grid[static_cast<std::size_t>(w.row) * p.columns + w.column] = w.readout;
        if (first || w.readout < lo) lo = w.readout;
        if (first || w.readout > hi) hi = w.readout;
        first = false;
        if (w.excluded) {
            exX.push_back((w.column + 0.5) / p.columns);
            exY.push_back(1.0 - (w.row + 0.5) / p.rows);
        }
    }
    if (first) return;
    if (hi <= lo) hi = lo + 1.0;
    if (ImPlot::BeginPlot("##assay-heat", ImVec2(-1, 240), ImPlotFlags_NoLegend)) {
        ImPlot::SetupAxes("column", "row",
                          ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoTickLabels,
                          ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoTickLabels);
        ImPlot::PlotHeatmap("readout", grid.data(), p.rows, p.columns, lo, hi, nullptr,
                            ImPlotPoint(0, 0), ImPlotPoint(1, 1));
        if (!exX.empty()) {
            ImPlot::PlotScatter("excluded", exX.data(), exY.data(),
                                static_cast<int>(exX.size()),
                                hollowMarkerSpec(ImPlotMarker_Square, 6.0f));
        }
        ImPlot::EndPlot();
    }
    ImGui::TextDisabled("Heat map spans %s to %s %s; hollow squares are wells an opt-in "
                        "outlier rule excluded (they are never deleted).",
                        sci(lo).c_str(), sci(hi).c_str(),
                        p.readoutUnit.empty() ? "readout units" : p.readoutUnit.c_str());
}

// The log10 interval the FIT reported for its half-maximal concentration, profile
// first, Wald second, nothing third. Nothing is a legitimate answer: an interval
// invented by the panel would be indistinguishable from one the fitter earned.
bool fitLog10Interval(const FitResult& f, double& lo, double& hi) {
    for (const auto& p : f.parameters) {
        if (p.name != "log10EC50" && p.name != "log10C") continue;
        if (p.profileComputed && p.profileUpper > p.profileLower) {
            lo = p.profileLower;
            hi = p.profileUpper;
            return true;
        }
        if (p.value.error > 0.0) {
            lo = p.value.value - 1.959963984540054 * p.value.error;
            hi = p.value.value + 1.959963984540054 * p.value.error;
            return true;
        }
    }
    return false;
}

void assayCurve(const std::vector<Well>& series, const FitResult& fit) {
    std::vector<double> ix, iy, ex, ey;
    for (const auto& w : series) {
        if (!(w.concentration > 0.0)) continue;
        const double x = std::log10(w.concentration);
        if (w.excluded) {
            ex.push_back(x);
            ey.push_back(w.readout);
        } else {
            ix.push_back(x);
            iy.push_back(w.readout);
        }
    }
    std::vector<double> fx, fy;
    for (std::size_t i = 0; i < fit.fittedX.size() && i < fit.fittedY.size(); ++i) {
        if (!(fit.fittedX[i] > 0.0)) continue;
        fx.push_back(std::log10(fit.fittedX[i]));
        fy.push_back(fit.fittedY[i]);
    }
    if (ImPlot::BeginPlot("##assay-curve", ImVec2(-1, 250))) {
        ImPlot::SetupAxes("log10 concentration (mol/L)", "readout");
        double lo = 0.0, hi = 0.0;
        if (fitLog10Interval(fit, lo, hi)) {
            // The band is the reported CONFIDENCE INTERVAL ON THE MIDPOINT, drawn
            // where it lives - on the concentration axis. It is not a pointwise
            // prediction band, and labelling it as one would overstate the fit.
            const double bx[2] = {lo, hi};
            const double top[2] = {1e30, 1e30};
            const double bot[2] = {-1e30, -1e30};
            ImPlot::PlotShaded("95% CI on the midpoint", bx, bot, top, 2);
            const double mid = 0.5 * (lo + hi);
            ImPlot::PlotInfLines("reported midpoint", &mid, 1);
        }
        if (!fx.empty())
            ImPlot::PlotLine("fit", fx.data(), fy.data(), static_cast<int>(fx.size()));
        ImPlot::PlotScatter("measured", ix.data(), iy.data(), static_cast<int>(ix.size()),
                            markerSpec(ImPlotMarker_Circle, 5.0f));
        if (!ex.empty()) {
            ImPlot::PlotScatter("excluded", ex.data(), ey.data(), static_cast<int>(ex.size()),
                                hollowMarkerSpec(ImPlotMarker_Circle, 6.0f));
        }
        ImPlot::EndPlot();
    }

    if (!fit.residuals.empty() && !ix.empty()) {
        std::vector<double> rx;
        for (std::size_t i = 0; i < fit.residuals.size() && i < ix.size(); ++i)
            rx.push_back(ix[i]);
        if (ImPlot::BeginPlot("##assay-resid", ImVec2(-1, 140), ImPlotFlags_NoLegend)) {
            ImPlot::SetupAxes("log10 concentration (mol/L)", "residual");
            const double zero = 0.0;
            ImPlotSpec zeroSpec;
            zeroSpec.Flags = ImPlotInfLinesFlags_Horizontal;
            ImPlot::PlotInfLines("##zero", &zero, 1, zeroSpec);
            ImPlot::PlotScatter("residual", rx.data(), fit.residuals.data(),
                                static_cast<int>(std::min(rx.size(), fit.residuals.size())),
                                markerSpec(ImPlotMarker_Diamond, 4.0f));
            ImPlot::EndPlot();
        }
        ImGui::TextDisabled("Residuals belong next to the curve: a high R-squared with "
                            "structured residuals is a wrong model that fits well.");
    }
}

// Sensorgram: response vs time, one trace per analyte concentration.
void assaySensorgram(const std::vector<Well>& wells) {
    std::map<double, std::pair<std::vector<double>, std::vector<double>>> byConc;
    for (const auto& w : wells) byConc[w.concentration].first.push_back(w.timeS);
    for (const auto& w : wells) byConc[w.concentration].second.push_back(w.readout);
    if (byConc.empty()) return;
    if (ImPlot::BeginPlot("##assay-spr", ImVec2(-1, 220))) {
        ImPlot::SetupAxes("time (s)", "response (RU)");
        for (auto& [conc, xy] : byConc) {
            const std::string label = sci(conc) + " M";
            ImPlot::PlotLine(label.c_str(), xy.first.data(), xy.second.data(),
                             static_cast<int>(std::min(xy.first.size(), xy.second.size())));
        }
        ImPlot::EndPlot();
    }
}

void assayMelt(const std::vector<Well>& wells) {
    std::vector<double> t, y;
    for (const auto& w : wells) {
        t.push_back(w.temperatureC);
        y.push_back(w.readout);
    }
    if (t.size() < 2) return;
    if (ImPlot::BeginPlot("##assay-dsf", ImVec2(-1, 200))) {
        ImPlot::SetupAxes("temperature (C)", "signal");
        ImPlot::PlotLine("melt", t.data(), y.data(), static_cast<int>(t.size()));
        ImPlot::EndPlot();
    }
}

// The [S] x [I] velocity matrix as a heat map. The modality verdict is a global fit
// over this whole surface, not over any one row of it, which is why the matrix is
// what gets drawn.
void assayInhibitionMatrix(const std::vector<Well>& wells) {
    std::vector<double> subs, inhibs;
    for (const auto& w : wells) {
        subs.push_back(w.concentration);
        inhibs.push_back(std::atof(w.seriesId.c_str()));
    }
    std::sort(subs.begin(), subs.end());
    subs.erase(std::unique(subs.begin(), subs.end()), subs.end());
    std::sort(inhibs.begin(), inhibs.end());
    inhibs.erase(std::unique(inhibs.begin(), inhibs.end()), inhibs.end());
    if (subs.size() < 2 || inhibs.empty()) return;
    std::vector<double> grid(subs.size() * inhibs.size(), 0.0);
    double lo = 0.0, hi = 0.0;
    bool first = true;
    for (const auto& w : wells) {
        const auto si = std::lower_bound(subs.begin(), subs.end(), w.concentration) - subs.begin();
        const double iv = std::atof(w.seriesId.c_str());
        const auto ii = std::lower_bound(inhibs.begin(), inhibs.end(), iv) - inhibs.begin();
        if (si >= static_cast<long>(subs.size()) || ii >= static_cast<long>(inhibs.size()))
            continue;
        grid[static_cast<std::size_t>(ii) * subs.size() + si] = w.readout;
        if (first || w.readout < lo) lo = w.readout;
        if (first || w.readout > hi) hi = w.readout;
        first = false;
    }
    if (hi <= lo) hi = lo + 1.0;
    if (ImPlot::BeginPlot("##assay-matrix", ImVec2(-1, 220), ImPlotFlags_NoLegend)) {
        ImPlot::SetupAxes("[S] (ascending)", "[I] (ascending)",
                          ImPlotAxisFlags_NoGridLines, ImPlotAxisFlags_NoGridLines);
        ImPlot::PlotHeatmap("velocity", grid.data(), static_cast<int>(inhibs.size()),
                            static_cast<int>(subs.size()), lo, hi, nullptr, ImPlotPoint(0, 0),
                            ImPlotPoint(1, 1));
        ImPlot::EndPlot();
    }
    ImGui::TextDisabled("A Lineweaver-Burk plot is available as a DIAGNOSTIC only; no fitter "
                        "in BioCAD attaches to a transformed axis, because the transform "
                        "distorts the error structure it would be fitting.");
}

// ITC in the standard layout: raw differential power against time on top, the
// integrated per-injection heat against molar ratio underneath. Anything else is
// not an ITC figure a calorimetrist can check.
void assayItc(const std::vector<Well>& wells) {
    std::vector<double> t, p;
    for (const auto& w : wells) {
        t.push_back(w.timeS);
        p.push_back(w.readout);
    }
    if (t.size() < 2) return;
    if (ImPlot::BeginPlot("##assay-itc-raw", ImVec2(-1, 180), ImPlotFlags_NoLegend)) {
        ImPlot::SetupAxes("time (s)", "differential power (ucal/s)");
        ImPlot::PlotLine("power", t.data(), p.data(), static_cast<int>(t.size()));
        ImPlot::EndPlot();
    }
    // Integrated heat per injection: the trapezoid of power over each injection's
    // own time span, which is the only integration the raw trace supports.
    std::vector<double> ratio, heat;
    double acc = 0.0;
    for (std::size_t i = 1; i < t.size(); ++i) {
        acc += 0.5 * (p[i] + p[i - 1]) * (t[i] - t[i - 1]);
        ratio.push_back(static_cast<double>(i));
        heat.push_back(acc);
    }
    if (ImPlot::BeginPlot("##assay-itc-iso", ImVec2(-1, 180), ImPlotFlags_NoLegend)) {
        ImPlot::SetupAxes("injection number (molar ratio needs the cell composition)",
                          "cumulative heat (ucal)");
        ImPlot::PlotScatter("integrated", ratio.data(), heat.data(),
                            static_cast<int>(ratio.size()),
                            markerSpec(ImPlotMarker_Square, 4.0f));
        ImPlot::EndPlot();
    }
    ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed),
                       "The x axis is the injection index, not the molar ratio: the ratio needs "
                       "the cell volume and the macromolecule and titrant concentrations, which "
                       "are experiment metadata and not well fields.");
}

}  // namespace

void assayWorkbench(AppShell& shell) {
    Services& s = shell.services();
    if (!s.assay) return;

    static std::vector<char>              csv(1 << 16);
    static bool                           seeded = false;
    static std::optional<AssayDataset>    ds;
    static std::string                    importError;
    static int                            plateIdx = 0;
    static int                            modelIdx = 0;
    static bool                           robust = false;
    static std::string                    seriesId;
    static std::optional<FitResult>        fit;
    static std::optional<ModelComparison>  modality;
    if (!seeded) {
        std::snprintf(csv.data(), csv.size(), "%s", kAssaySampleCsv);
        seeded = true;
    }

    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::Measured),
                       "MEASURED DATA IN, MODEL PARAMETERS OUT. Well readouts are measured and "
                       "stay measured; every fitted parameter is a model value with its error "
                       "bar. Excluded wells are hollowed, never deleted, and the rule that "
                       "excluded them is named.");
    ImGui::PopTextWrapPos();
    ImGui::Spacing();

    theme::sectionHeader("IMPORT (LONG CSV/TSV OR A 96/384/1536 GRID EXPORT)");
    ImGui::InputTextMultiline("##assay-csv", csv.data(), csv.size(), ImVec2(-1, 150));
    if (ImGui::Button("Import")) {
        ds = s.assay->import(std::string(csv.data()));
        importError = ds ? std::string() : "that text is not a plate table BioCAD can read";
        plateIdx = 0;
        fit.reset();
        modality.reset();
        seriesId.clear();
    }
    if (!importError.empty())
        ImGui::TextColored(theme::verdictColor(3), "%s", importError.c_str());
    if (!ds) {
        ImGui::TextDisabled("Paste a plate export and press Import. Recognised columns: "
                            "plate_id, well, role, sample_id, series_id, concentration, "
                            "conc_unit, replicate, readout, readout_unit, time_s, "
                            "temperature_c, excluded, exclusion_rule. Unknown columns survive "
                            "as metadata.");
        return;
    }
    ImGui::Text("Layout: %s   plates: %d", ds->detectedLayout.c_str(),
                static_cast<int>(ds->plates.size()));
    for (const auto& w : ds->warnings)
        ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic), "%s", w.c_str());
    if (ds->plates.empty()) return;

    plateIdx = std::clamp(plateIdx, 0, static_cast<int>(ds->plates.size()) - 1);
    if (ds->plates.size() > 1) {
        std::vector<const char*> names;
        for (const auto& p : ds->plates) names.push_back(p.id.c_str());
        ImGui::SetNextItemWidth(220.0f);
        ImGui::Combo("Plate", &plateIdx, names.data(), static_cast<int>(names.size()));
    }
    const Plate& plate = ds->plates[static_cast<std::size_t>(plateIdx)];

    theme::sectionHeader("PLATE HEAT MAP");
    assayHeatMap(plate);
    ImGui::Spacing();

    // ------------------------------------------------------------------ QC
    const QcReport qcr = s.assay->qc(plate);
    theme::sectionHeader("PLATE QC");
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(qcr.zPrime.provenance == Provenance::NotComputed
                           ? theme::provenanceColor(Provenance::NotComputed)
                           : (qcr.zPrime.value <= 0.0 ? theme::verdictColor(3)
                                                      : (qcr.zPrime.value < 0.5
                                                             ? theme::verdictColor(2)
                                                             : theme::verdictColor(0))),
                       "%s", qcr.interpretation.c_str());
    ImGui::PopTextWrapPos();
    drawQuantity("Z-prime", qcr.zPrime);
    drawQuantity("Robust Z-prime (median/MAD)", qcr.robustZPrime);
    drawQuantity("SSMD", qcr.ssmd);
    drawQuantity("Signal / background", qcr.signalToBackground);
    drawQuantity("Signal / noise", qcr.signalToNoise);
    drawQuantity("Positive control mean", qcr.positiveMean);
    drawQuantity("Positive control SD", qcr.positiveSd);
    drawQuantity("Negative control mean", qcr.negativeMean);
    drawQuantity("Negative control SD", qcr.negativeSd);
    drawQuantity("Positive control %CV", qcr.cvPositivePct);
    drawQuantity("Negative control %CV", qcr.cvNegativePct);
    drawQuantity("Edge effect p (Mann-Whitney)", qcr.edgeEffectP);
    drawQuantity("Row effect p (Kruskal-Wallis)", qcr.rowEffectP);
    drawQuantity("Column effect p (Kruskal-Wallis)", qcr.columnEffectP);
    for (const auto& w : qcr.warnings)
        ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic), "%s", w.c_str());
    ImGui::TextDisabled("Edge, row and column effects are REPORTED, never auto-corrected: "
                        "median-polishing a gradient away hides the pipetting problem that "
                        "caused it.");
    ImGui::Spacing();

    // -------------------------------------------------------------- fitting
    theme::sectionHeader("FIT ONE SERIES");
    std::vector<std::string> seriesIds;
    for (const auto& w : plate.wells)
        if (w.role == WellRole::Sample && !w.seriesId.empty() &&
            std::find(seriesIds.begin(), seriesIds.end(), w.seriesId) == seriesIds.end())
            seriesIds.push_back(w.seriesId);
    if (seriesIds.empty()) {
        ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed),
                           "No sample series on this plate: a fit needs wells with role "
                           "'sample' and a series_id.");
        return;
    }
    if (seriesId.empty()) seriesId = seriesIds.front();
    std::vector<const char*> sids;
    for (const auto& id : seriesIds) sids.push_back(id.c_str());
    int sidIdx = 0;
    for (std::size_t i = 0; i < seriesIds.size(); ++i)
        if (seriesIds[i] == seriesId) sidIdx = static_cast<int>(i);
    ImGui::SetNextItemWidth(180.0f);
    if (ImGui::Combo("Series", &sidIdx, sids.data(), static_cast<int>(sids.size())))
        seriesId = seriesIds[static_cast<std::size_t>(sidIdx)];
    ImGui::SameLine();
    ImGui::SetNextItemWidth(220.0f);
    ImGui::Combo("Model", &modelIdx, kAssayModelNames, IM_ARRAYSIZE(kAssayModelNames));
    ImGui::SameLine();
    ImGui::Checkbox("Tukey-biweight IRLS", &robust);

    std::vector<Well> series;
    for (const auto& w : plate.wells)
        if (w.role == WellRole::Sample && w.seriesId == seriesId) series.push_back(w);

    ImGui::SameLine();
    if (ImGui::Button("Fit"))
        fit = s.assay->fit(series, static_cast<AssayModel>(modelIdx), robust);
    ImGui::SameLine();
    if (ImGui::Button("Inhibition modality (whole plate)")) {
        std::vector<Well> matrix;
        for (const auto& w : plate.wells)
            if (w.role == WellRole::Sample) matrix.push_back(w);
        modality = s.assay->inhibitionModality(matrix);
    }

    const auto model = static_cast<AssayModel>(modelIdx);
    if (model == AssayModel::LangmuirKinetics || model == AssayModel::MassTransportKinetics)
        assaySensorgram(series);
    else if (model == AssayModel::BoltzmannMelt || model == AssayModel::TwoStateThermodynamic)
        assayMelt(series);
    else if (model == AssayModel::WisemanIsotherm)
        assayItc(series);

    if (fit) {
        if (model == AssayModel::FourParameterLogistic ||
            model == AssayModel::FiveParameterLogistic || model == AssayModel::Hill ||
            model == AssayModel::MichaelisMenten ||
            model == AssayModel::SubstrateInhibition ||
            model == AssayModel::MorrisonTightBinding)
            assayCurve(series, *fit);

        if (!fit->converged)
            ImGui::TextColored(theme::verdictColor(3), "Not fitted: %s", fit->note.c_str());

        // An EC50 the ladder never bracketed is GREY with the reason. It is not a
        // result, and colouring it like one is how an extrapolated potency ends up
        // in a slide deck.
        if (fit->extrapolated) {
            ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed),
                               "EC50 %s mol/L lies OUTSIDE the tested concentration range: "
                               "the curve was extrapolated past the highest or lowest well, so "
                               "this is a bound, not a measurement. Widen the ladder.",
                               sci(fit->derivedEc50.value).c_str());
        } else {
            drawQuantity("EC50 / midpoint", fit->derivedEc50);
        }
        drawQuantity("KD (kinetics)", fit->derivedKd);
        ImGui::Text("R-squared %.5f   AICc %.3f   rank %d/%d   condition %.3g   n = %d",
                    fit->rSquared, fit->aicc, static_cast<int>(fit->rank),
                    static_cast<int>(fit->parameters.size()), fit->conditionNumber,
                    static_cast<int>(fit->observations));
        if (fit->robust)
            ImGui::TextDisabled("Fitted by Tukey-biweight IRLS: down-weighted points are still "
                                "plotted.");
        if (ImGui::BeginTable("assay-params", 4,
                              ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
            ImGui::TableSetupColumn("Parameter", ImGuiTableColumnFlags_WidthFixed, 220.0f);
            ImGui::TableSetupColumn("Value", ImGuiTableColumnFlags_WidthFixed, 160.0f);
            ImGui::TableSetupColumn("Std. error", ImGuiTableColumnFlags_WidthFixed, 140.0f);
            ImGui::TableSetupColumn("Interval");
            ImGui::TableHeadersRow();
            for (const auto& p : fit->parameters) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::TextColored(theme::provenanceColor(p.value.provenance), "%s",
                                   p.name.c_str());
                ImGui::TableSetColumnIndex(1);
                ImGui::TextUnformatted(sci(p.value.value).c_str());
                ImGui::TableSetColumnIndex(2);
                ImGui::TextUnformatted(p.value.error > 0.0 ? sci(p.value.error).c_str() : "-");
                ImGui::TableSetColumnIndex(3);
                if (p.profileComputed)
                    ImGui::Text("[%s, %s] profile", sci(p.profileLower).c_str(),
                                sci(p.profileUpper).c_str());
                else
                    ImGui::TextDisabled("no profile interval requested");
            }
            ImGui::EndTable();
        }
        for (const auto& a : fit->assumptions) ImGui::BulletText("%s", a.c_str());
        for (const auto& w : fit->warnings)
            ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic), "%s", w.c_str());
    }
    ImGui::Spacing();

    if (modality) {
        theme::sectionHeader("INHIBITION MODALITY ([S] x [I] GLOBAL FIT)");
        std::vector<Well> matrix;
        for (const auto& w : plate.wells)
            if (w.role == WellRole::Sample) matrix.push_back(w);
        assayInhibitionMatrix(matrix);
        ImGui::TextColored(modality->decisive ? theme::provenanceColor(Provenance::Model)
                                             : theme::provenanceColor(Provenance::NotComputed),
                           "%s", modality->conclusion.c_str());
        ImGui::Text("delta AICc %.3f%s", modality->deltaAicc,
                    modality->decisive ? "" : " - under 2, so the modality is Unknown");
        for (const auto& c : modality->candidates)
            ImGui::BulletText("%s: AICc %.3f, R-squared %.4f", kAssayModelNames[0],
                              c.aicc, c.rSquared);
    }
}

// ------------------------------------------------------ Assay Design
void assayDesign(AppShell& shell) {
    Services& s = shell.services();
    if (!s.assay) return;

    static AssayDesignSpec spec = [] {
        AssayDesignSpec d;
        d.truthModel = AssayModel::FourParameterLogistic;
        // A, B, C, D in the letter order of the 5PL: signal at zero, slope, midpoint,
        // plateau. G is 1 for a 4PL.
        d.truthParameters = {100.0, 1.0, 1.0e-7, 0.0};
        for (int i = 0; i < 10; ++i) d.concentrations.push_back(1.0e-5 / std::pow(3.1623, i));
        d.replicates = 3;
        d.rows = 8;
        d.columns = 12;
        d.additiveNoiseSd = 2.0;
        d.proportionalNoiseCv = 0.03;
        d.pipettingCv = 0.02;
        d.plateGradientPct = 4.0;
        d.dmsoTolerancePct = 3.0;
        d.seed = 20260813;
        d.replicateRuns = 200;
        return d;
    }();
    static std::optional<AssayDesignReport> rep;

    // The dilution calculator's own state: it is arithmetic, not a simulation, so it
    // is live rather than behind the Simulate button.
    static double molarMass = 250.0;
    static double stockM = 0.01;
    static double stockVolumeMl = 1.0;
    static double topM = 1.0e-5;
    static double fold = 3.1623;
    static int    steps = 10;
    static double wellUl = 100.0;

    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic),
                       "EXPERIMENTAL DESIGN, NOT A DOSE AND NOT A PROCEDURE. This panel "
                       "forward-simulates plates from a truth model YOU state, pushes each one "
                       "through the same import -> QC -> fit path real data takes, and reports "
                       "what the design would recover. A concentration ladder for a plate is "
                       "not a regimen for a person.");
    ImGui::PopTextWrapPos();
    ImGui::Spacing();

    theme::sectionHeader("TRUTH MODEL AND ERROR STRUCTURE (YOUR STATED BELIEF)");
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputDouble("Signal at zero (A)", &spec.truthParameters[0], 0.0, 0.0, "%.3f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputDouble("Slope (B)", &spec.truthParameters[1], 0.0, 0.0, "%.3f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(160.0f);
    ImGui::InputDouble("Midpoint C (mol/L)", &spec.truthParameters[2], 0.0, 0.0, "%.3e");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputDouble("Plateau (D)", &spec.truthParameters[3], 0.0, 0.0, "%.3f");

    ImGui::SetNextItemWidth(130.0f);
    ImGui::InputDouble("Additive SD", &spec.additiveNoiseSd, 0.0, 0.0, "%.3f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(130.0f);
    ImGui::InputDouble("Proportional CV", &spec.proportionalNoiseCv, 0.0, 0.0, "%.4f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(130.0f);
    ImGui::InputDouble("Pipetting CV", &spec.pipettingCv, 0.0, 0.0, "%.4f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(130.0f);
    ImGui::InputDouble("Gradient %", &spec.plateGradientPct, 0.0, 0.0, "%.2f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(130.0f);
    ImGui::InputDouble("DMSO loss %", &spec.dmsoTolerancePct, 0.0, 0.0, "%.2f");

    ImGui::SetNextItemWidth(110.0f);
    ImGui::InputInt("Replicates", &spec.replicates);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(110.0f);
    ImGui::InputInt("Runs", &spec.replicateRuns);
    ImGui::SameLine();
    int seedI = static_cast<int>(spec.seed);
    ImGui::SetNextItemWidth(140.0f);
    if (ImGui::InputInt("Seed", &seedI)) spec.seed = static_cast<std::uint64_t>(std::max(0, seedI));
    ImGui::SameLine();
    if (ImGui::Button("Simulate")) rep = s.assay->simulate(spec);
    ImGui::TextDisabled("The same seed reproduces the same report byte for byte: the RNG is "
                        "PCG64-DXSM with a Box-Muller normal, both implemented in BioCAD, "
                        "because std::normal_distribution is not specified bit-exactly.");
    ImGui::Spacing();

    theme::sectionHeader("LADDER");
    ImGui::Text("%d concentrations, %s down to %s mol/L",
                static_cast<int>(spec.concentrations.size()),
                spec.concentrations.empty() ? "-" : sci(spec.concentrations.front()).c_str(),
                spec.concentrations.empty() ? "-" : sci(spec.concentrations.back()).c_str());
    ImGui::Spacing();

    if (rep) {
        theme::sectionHeader("WHAT THIS DESIGN WOULD RECOVER");
        // Coverage first, deliberately: it is the number that decides whether the CI
        // this design reports means anything.
        statCard("CI COVERAGE",
                 rep->empiricalCoveragePct.provenance == Provenance::NotComputed
                     ? std::string("n/a")
                     : f2(rep->empiricalCoveragePct.value) + "%",
                 "nominal 95%", 200.0f);
        ImGui::SameLine();
        statCard("CONVERGED", f2(rep->convergenceRatePct.value) + "%", "of simulated runs",
                 180.0f);
        ImGui::SameLine();
        statCard("MEDIAN Z'", f2(rep->medianZPrime.value), ">= 0.5 excellent, <= 0 unusable",
                 220.0f);
        ImGui::SameLine();
        statCard("MEDIAN EC50", sci(rep->medianEc50.value), "mol/L recovered", 200.0f);
        ImGui::SameLine();
        statCard("CI WIDTH", f2(rep->ec50CiWidthLog10.value), "log10 mol/L", 170.0f);
        ImGui::Spacing();

        drawQuantity("Empirical CI coverage", rep->empiricalCoveragePct);
        drawQuantity("Convergence rate", rep->convergenceRatePct);
        drawQuantity("Median Z-prime", rep->medianZPrime);
        drawQuantity("Median recovered EC50", rep->medianEc50);
        drawQuantity("Median CI width", rep->ec50CiWidthLog10);
        // The reference measurement, so a reader can tell an ill-conditioned design from a
        // broken fitter: the same machinery was measured at 94.4% coverage against a
        // nominal 95% over 1000 runs (truth A = 100, B = 1, C = 1e-7 mol/L, D = 0, a
        // ten-point half-log ladder from 1e-5, three replicates, additive SD 2.0).
        ImGui::TextDisabled("Reference: this machinery measured 94.4%% empirical coverage at a "
                            "nominal 95%% over 1000 runs on a well-conditioned ten-point "
                            "half-log ladder. Coverage far below that is a property of the "
                            "design you entered, not of the fitter.");

        if (rep->recoveredEc50.size() > 4) {
            std::vector<double> logs;
            for (double v : rep->recoveredEc50)
                if (v > 0.0) logs.push_back(std::log10(v));
            if (ImPlot::BeginPlot("##design-hist", ImVec2(-1, 200), ImPlotFlags_NoLegend)) {
                ImPlot::SetupAxes("log10 recovered EC50 (mol/L)", "runs");
                ImPlot::PlotHistogram("recovered", logs.data(), static_cast<int>(logs.size()));
                const double truth = std::log10(spec.truthParameters[2]);
                ImPlot::PlotInfLines("truth", &truth, 1);
                ImPlot::EndPlot();
            }
            ImGui::TextDisabled("The vertical line is the truth the plates were generated "
                                "from. A histogram centred off that line is bias, not noise.");
        }

        theme::sectionHeader("D-OPTIMAL LADDER (FEDOROV COORDINATE EXCHANGE)");
        ImGui::Text("D-efficiency gain over the entered ladder: %.4f", rep->dOptimalityGain);
        std::string opt;
        for (double c : rep->optimalConcentrations) opt += sci(c) + "  ";
        ImGui::TextWrapped("%s", opt.c_str());
        ImGui::TextDisabled("Candidates are achievable ladder points only (top / sqrt(fold)^j): "
                            "optimising over a continuum would return concentrations nobody can "
                            "pipette.");
        ImGui::Spacing();

        theme::sectionHeader("WHAT THIS SIMULATION RESTS ON");
        for (const auto& a : rep->assumptions) ImGui::BulletText("%s", a.c_str());
        for (const auto& w : rep->warnings)
            ImGui::TextColored(theme::verdictColor(2), "%s", w.c_str());
        ImGui::Spacing();
    }

    theme::sectionHeader("DILUTION AND MASS CALCULATOR");
    ImGui::SetNextItemWidth(150.0f);
    ImGui::InputDouble("Molar mass (g/mol)", &molarMass, 0.0, 0.0, "%.3f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(150.0f);
    ImGui::InputDouble("Stock (mol/L)", &stockM, 0.0, 0.0, "%.4g");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(150.0f);
    ImGui::InputDouble("Stock volume (mL)", &stockVolumeMl, 0.0, 0.0, "%.3f");
    ImGui::Text("Weigh out %.6g mg of solid.",
                1000.0 * assay::massForStock(molarMass, stockM, stockVolumeMl / 1000.0));

    ImGui::SetNextItemWidth(150.0f);
    ImGui::InputDouble("Top well (mol/L)", &topM, 0.0, 0.0, "%.4g");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputDouble("Fold per step", &fold, 0.0, 0.0, "%.4f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(110.0f);
    ImGui::InputInt("Steps", &steps);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputDouble("Well volume (uL)", &wellUl, 0.0, 0.0, "%.2f");

    const assay::DilutionPlan plan =
        assay::serialDilution(stockM, topM, fold, steps, wellUl);
    if (ImGui::BeginTable("design-dil", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Step", ImGuiTableColumnFlags_WidthFixed, 60.0f);
        ImGui::TableSetupColumn("Concentration (mol/L)", ImGuiTableColumnFlags_WidthFixed, 200.0f);
        ImGui::TableSetupColumn("Transfer (uL)", ImGuiTableColumnFlags_WidthFixed, 140.0f);
        ImGui::TableSetupColumn("Diluent (uL)", ImGuiTableColumnFlags_WidthFixed, 140.0f);
        ImGui::TableSetupColumn("Fold");
        ImGui::TableHeadersRow();
        for (const auto& st : plan.steps) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%d", st.index + 1);
            ImGui::TableSetColumnIndex(1);
            ImGui::TextUnformatted(sci(st.concentration).c_str());
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%.3f", st.transferVolumeUl);
            ImGui::TableSetColumnIndex(3);
            ImGui::Text("%.3f", st.diluentVolumeUl);
            ImGui::TableSetColumnIndex(4);
            ImGui::Text("%.4f", st.fold);
        }
        ImGui::EndTable();
    }
    ImGui::Text("Pipetting error at the last step compounds to %.3f x the CV of one transfer.",
                plan.compoundedCvAtLastStep);
    for (const auto& w : plan.warnings)
        ImGui::TextColored(theme::verdictColor(2), "%s", w.c_str());
}

// ---- Protein 3D scene assembly (shared by the Structure3D panel's two modes) ----
namespace {

// Atomic number of a PDB element symbol, through the one element table the tree already has
// (chem::parseFormula). Cached because a 100k-atom structure would otherwise parse "C" 100k
// times, and duplicated as a second symbol table is how two tables disagree.
int elementZ(const std::string& symbol) {
    static std::unordered_map<std::string, int> cache;
    if (symbol.empty()) return 6;
    const auto it = cache.find(symbol);
    if (it != cache.end()) return it->second;
    int z = 6;   // an unrecognised symbol renders as carbon grey rather than vanishing
    std::string tidy;
    for (char c : symbol)
        if (c != ' ') tidy += c;
    if (!tidy.empty()) {
        tidy[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(tidy[0])));
        for (std::size_t i = 1; i < tidy.size(); ++i)
            tidy[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(tidy[i])));
        if (const auto f = chem::parseFormula(tidy); f && f->ok && f->terms.size() == 1)
            z = f->terms.front().z;
    }
    cache[symbol] = z;
    return z;
}

// Spheres and sticks for a whole model, using the residue-template bond graph. The bonds are the
// point: without Connectivity this mode could only draw disconnected spheres, because the file
// has no bond list.
render::MolScene proteinSticksScene(const bio::Model& model, const bio::ConnectivityResult& conn) {
    render::MolScene scene;
    // Flat index per (chain, residue, atom) so a StructureBond can address the sphere list.
    std::map<std::tuple<int, int, int>, std::size_t> index;
    for (std::size_t ci = 0; ci < model.chains.size(); ++ci) {
        const bio::Chain& c = model.chains[ci];
        for (std::size_t ri = 0; ri < c.residues.size(); ++ri) {
            const bio::Residue& r = c.residues[ri];
            for (std::size_t ai = 0; ai < r.atoms.size(); ++ai) {
                const bio::Atom& a = r.atoms[ai];
                render::AtomInst inst;
                inst.x = static_cast<float>(a.x);
                inst.y = static_cast<float>(a.y);
                inst.z = static_cast<float>(a.z);
                inst.r = 0.30f;
                inst.rgba = render::cpkColor(elementZ(a.element));
                index[{static_cast<int>(ci), static_cast<int>(ri), static_cast<int>(ai)}] =
                    scene.atoms.size();
                scene.atoms.push_back(inst);
            }
        }
    }
    for (const bio::StructureBond& b : conn.bonds) {
        const auto ia = index.find({b.a.chain, b.a.residue, b.a.atom});
        const auto ib = index.find({b.b.chain, b.b.residue, b.b.atom});
        if (ia == index.end() || ib == index.end()) continue;
        const render::AtomInst& A = scene.atoms[ia->second];
        const render::AtomInst& B = scene.atoms[ib->second];
        render::BondInst stick;
        stick.ax = A.x; stick.ay = A.y; stick.az = A.z;
        stick.bx = B.x; stick.by = B.y; stick.bz = B.z;
        stick.rgbaA = A.rgba;
        stick.rgbaB = B.rgba;
        // A disulfide is drawn thicker so the one bond the geometry (not a template) produced is
        // visible as such.
        stick.radius = b.kind == bio::BondKind::Disulfide ? 0.16f : 0.11f;
        scene.bonds.push_back(stick);
    }
    render::computeSceneBounds(scene);
    return scene;
}

// The cartoon: bio::CartoonMesh straight into the generic indexed-mesh path. The palette index
// is resolved here because the renderer's vertex carries a packed colour, not an index.
render::MolScene proteinCartoonScene(const bio::CartoonMesh& mesh) {
    render::MolScene scene;
    scene.mesh.vertices.reserve(mesh.vertices.size());
    for (const bio::CartoonVertex& v : mesh.vertices) {
        render::MeshVertexRgba out;
        out.px = v.px; out.py = v.py; out.pz = v.pz;
        out.nx = v.nx; out.ny = v.ny; out.nz = v.nz;
        out.rgba = v.colorIndex < mesh.palette.size() ? mesh.palette[v.colorIndex] : 0xFFFFFFFFu;
        scene.mesh.vertices.push_back(out);
    }
    scene.mesh.indices = mesh.indices;
    render::computeSceneBounds(scene);
    return scene;
}

}  // namespace

// ------------------------------------------------------ Protein Structure
void proteinStructure(AppShell& shell) {
    Services& s = shell.services();
    if (!s.structure) return;

    static char pathBuf[1024] = "";
    static std::optional<bio::Structure> loaded;
    static std::string loadError;
    // The Structure module's load() returns coordinates only, and the cartoon needs the HELIX /
    // SHEET (or _struct_conf / _struct_sheet_range) records. They are read here with the same
    // readers, from the same file, rather than by inventing a third secondary-structure source.
    static bio::Annotations annotations;
    static bool cartoonMode = true;
    static ViewerUiState viewerUi;
    static std::string sceneKey;                 // structure + mode currently built
    static render::MolScene scene;
    static bio::CartoonMesh cartoon;
    static bio::ConnectivityResult conn;

    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::Measured),
                       "This panel reads a LOCAL .pdb / .cif file (including one already "
                       "downloaded into the cache). It does not fetch anything by itself.");
    ImGui::PopTextWrapPos();
    ImGui::Spacing();

    ImGui::SetNextItemWidth(-140.0f);
    ImGui::InputTextWithHint("##structpath", "Path to a .pdb / .cif file...", pathBuf,
                             sizeof(pathBuf));
    ImGui::SameLine();
    if (ImGui::Button("Load", ImVec2(120, 0))) {
        loaded = s.structure->load(std::filesystem::path(pathBuf));
        loadError = loaded ? std::string()
                           : "Could not read '" + std::string(pathBuf) +
                                 "'. Expected an existing .pdb, .ent, .cif or .mmcif file.";
        annotations = bio::Annotations{};
        sceneKey.clear();
        if (loaded) {
            const std::filesystem::path p(pathBuf);
            std::string ext = p.extension().string();
            for (char& c : ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            if (ext == ".cif" || ext == ".mmcif")
                (void)bio::readCifFile(p, &annotations);
            else
                (void)bio::readPdbFile(p, &annotations);
        }
    }
    if (!loadError.empty()) {
        ImGui::TextColored(theme::verdictColor(3), "%s", loadError.c_str());
    }
    if (!loaded) {
        ImGui::TextDisabled("No structure loaded.");
        return;
    }

    const bio::Structure& st = *loaded;
    const bio::Model* m = st.model(1);
    const std::size_t chains = m ? m->chains.size() : 0;
    std::size_t residues = 0;
    if (m) for (const auto& c : m->chains) residues += c.residues.size();

    statCard("MODELS", std::to_string(st.models.size()), "MODEL records", 150.0f);
    ImGui::SameLine();
    statCard("CHAINS", std::to_string(chains), "in model 1", 150.0f);
    ImGui::SameLine();
    statCard("RESIDUES", std::to_string(residues), "in model 1", 170.0f);
    ImGui::SameLine();
    statCard("ATOMS", std::to_string(st.atomCount()), "all models", 170.0f);
    ImGui::Spacing();

    // Residue numbering is ambiguous and the two schemes disagree in most mmCIF
    // entries, so the panel states which one it is showing every single time.
    ImGui::TextColored(theme::provenanceColor(Provenance::Measured),
                       "Residue numbering shown below: AUTHOR (auth_seq_id) - the numbering "
                       "papers cite. mmCIF label_seq_id is NOT shown here.");
    ImGui::Text("Entry id: %s", st.id.empty() ? "(none in file)" : st.id.c_str());
    ImGui::Text("Source: %s", st.source.c_str());
    ImGui::Spacing();

    if (!st.warnings.empty()) {
        theme::sectionHeader("PARSE WARNINGS (RECOVERED, NOT FATAL)");
        for (const auto& w : st.warnings)
            ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(theme::kWarn), "%s", w.c_str());
        ImGui::Spacing();
    }

    theme::sectionHeader("3D VIEW");
    if (m) {
        // The two modes need different geometry, so each is built once and cached under a key
        // that includes the mode: rebuilding a 77k-triangle ribbon every frame would be the one
        // avoidable cost in this panel.
        const std::string wanted = st.source + (cartoonMode ? "|cartoon" : "|sticks");
        if (wanted != sceneKey) {
            conn = bio::connect(*m);
            if (cartoonMode) {
                cartoon = bio::buildCartoon(*m, annotations);
                scene = proteinCartoonScene(cartoon);
            } else {
                cartoon = bio::CartoonMesh{};
                scene = proteinSticksScene(*m, conn);
            }
            sceneKey = wanted;
            viewerUi.lastKey.clear();   // re-frame: the two modes have different extents
        }

        if (ImGui::RadioButton("Cartoon (ribbon)", cartoonMode)) cartoonMode = true;
        ImGui::SameLine();
        if (ImGui::RadioButton("Spheres and sticks", !cartoonMode)) cartoonMode = false;
        ImGui::SameLine();
        if (ImGui::SmallButton("Reset view")) viewerUi.lastKey.clear();
        ImGui::SameLine();
        ImGui::TextDisabled("(drag: orbit  -  wheel: zoom  -  right-drag: pan)");

        sceneViewport(shell, scene, viewerUi, sceneKey, 420.0f);

        if (cartoonMode) {
            ImGui::Text("Ribbon: %zu residues, %zu vertices, %zu triangles.", cartoon.residues,
                        cartoon.vertices.size(), cartoon.triangleCount());
            // Secondary structure is the records', not this panel's: an entry with no HELIX or
            // SHEET record renders as all-coil, and saying so is the difference between "this
            // protein has no helices" and "this file did not list any".
            if (annotations.helices.empty() && annotations.strands.empty()) {
                ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(theme::kWarn),
                                   "No HELIX / SHEET records in this file: every residue is drawn "
                                   "as coil. That is missing annotation, not absent structure.");
            } else {
                ImGui::TextDisabled("Secondary structure from the file's own HELIX / SHEET "
                                    "records (%zu helices, %zu strands).",
                                    annotations.helices.size(), annotations.strands.size());
            }
            for (const auto& w : cartoon.warnings)
                ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(theme::kWarn), "%s", w.c_str());
        } else {
            ImGui::Text("Bonds: %zu from residue templates, %zu polymer links, %zu disulfides.",
                        conn.diagnostics.templateBonds, conn.diagnostics.linkBonds,
                        conn.diagnostics.disulfides);
        }
        // Connectivity diagnostics are shown in BOTH modes: what could not be bonded is a fact
        // about the file, not about the current representation.
        ImGui::Text("Unbonded atoms: %zu. Template atoms absent from the file: %zu.",
                    conn.diagnostics.unbondedAtoms, conn.diagnostics.missingTemplateAtoms);
        for (const auto& w : conn.diagnostics.warnings)
            ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(theme::kWarn), "%s", w.c_str());
        ImGui::Spacing();
    }

    theme::sectionHeader("SOLVENT-ACCESSIBLE SURFACE AREA");
    // drawQuantity prints the full parameter string (probe radius, point count,
    // radii set, hydrogen policy) that the module put in `source`: a SASA without
    // them cannot be reproduced or compared against another tool.
    drawQuantity("Total SASA", s.structure->sasa(st));
    ImGui::Spacing();

    theme::sectionHeader("CHAINS (AUTHOR NUMBERING)");
    if (m && ImGui::BeginTable("bio-chains", 4,
                               ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Chain", ImGuiTableColumnFlags_WidthFixed, 70.0f);
        ImGui::TableSetupColumn("Residues", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("First..last auth #", ImGuiTableColumnFlags_WidthFixed, 160.0f);
        ImGui::TableSetupColumn("Polymer sequence (one-letter)");
        ImGui::TableHeadersRow();
        for (const auto& c : m->chains) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextUnformatted(c.id.c_str());
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%zu", c.residues.size());
            ImGui::TableSetColumnIndex(2);
            if (c.residues.empty()) {
                ImGui::TextDisabled("-");
            } else {
                ImGui::Text("%d%c..%d%c", c.residues.front().authSeqId,
                            c.residues.front().insertionCode,
                            c.residues.back().authSeqId, c.residues.back().insertionCode);
            }
            ImGui::TableSetColumnIndex(3);
            const std::vector<char> seq = bio::sequenceOf(c);
            if (seq.empty()) {
                ImGui::TextDisabled("(no polymer residues - ligand, water or ion chain)");
            } else {
                ImGui::TextWrapped("%s", std::string(seq.begin(), seq.end()).c_str());
            }
        }
        ImGui::EndTable();
    }
}

// ------------------------------------------------------- Variant Analysis
// Three surfaces, and the order is deliberate: the homolog set FIRST, because the
// conservation numbers are meaningless without it; the substitution scores with
// their published thresholds printed beside them; and the rebuilt side chain last,
// badged model, with no energy anywhere on the panel. Below the homolog minimum the
// conservation track is not drawn at all - a greyed-out track invites the reader to
// squint at it, whereas a sentence saying "15 required, 4 supplied" does not.
void variants(AppShell& shell) {
    Services& s = shell.services();
    if (!s.variants) return;

    static char queryBuf[8192] =
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR";
    static char homologBuf[65536] = "";
    static char structPath[1024] = "";
    static int  position = 1;
    static char mutant = 'A';
    static std::optional<bio::Structure> structure;
    static std::string structureError;

    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic),
                       "Conservation scores describe the homolog set YOU supply. They are not a "
                       "pathogenicity call, not a clinical interpretation, and not advice.");
    ImGui::PopTextWrapPos();
    ImGui::Spacing();

    theme::sectionHeader("QUERY AND HOMOLOGS (ONE SEQUENCE PER LINE)");
    ImGui::TextUnformatted("Query");
    ImGui::InputTextMultiline("##varquery", queryBuf, sizeof(queryBuf), ImVec2(-1, 60));
    ImGui::TextUnformatted("Homologs");
    ImGui::InputTextMultiline("##varhomologs", homologBuf, sizeof(homologBuf), ImVec2(-1, 110));

    const auto clean = [](const std::string& raw) {
        std::string out;
        for (char c : raw)
            if (std::isalpha(static_cast<unsigned char>(c)))
                out += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        return out;
    };
    const std::string query = clean(queryBuf);
    std::vector<std::string> homologs;
    {
        std::stringstream in{std::string(homologBuf)};
        std::string line;
        while (std::getline(in, line)) {
            if (!line.empty() && line[0] == '>') continue;   // a FASTA header is not sequence
            const std::string h = clean(line);
            if (!h.empty()) homologs.push_back(h);
        }
    }
    if (query.empty()) {
        ImGui::TextDisabled("Enter a query sequence.");
        return;
    }

    const ConservationProfile profile = s.variants->conservation(query, homologs);

    const std::string minimumSub =
        "supplied; " + std::to_string(profile.minimumHomologsRequired) + " required";
    statCard("HOMOLOGS", std::to_string(profile.homologs.sequenceCount), minimumSub.c_str(),
             230.0f);
    ImGui::SameLine();
    statCard("EFFECTIVE", std::to_string(profile.homologs.effectiveSequenceCount),
             "after identity clustering", 210.0f);
    ImGui::SameLine();
    statCard("MEDIAN IDENTITY", f2(profile.homologs.medianIdentityPct) + "%", "to the query",
             220.0f);
    ImGui::SameLine();
    statCard("RANGE", f2(profile.homologs.minIdentityPct) + " - " +
                          f2(profile.homologs.maxIdentityPct) + "%",
             "min - max identity", 220.0f);
    ImGui::Spacing();

    for (const auto& w : profile.warnings) ImGui::TextWrapped("%s", w.c_str());
    ImGui::Spacing();

    theme::sectionHeader("CONSERVATION TRACK");
    if (!profile.usable) {
        // No track. See the comment at the top of this function.
        ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed),
                           "Not computed: %d homolog(s) supplied, %d required.",
                           profile.homologs.sequenceCount, profile.minimumHomologsRequired);
        ImGui::TextWrapped(
            "A column distribution estimated from a handful of ad-hoc sequences is dominated by "
            "the pseudocount rather than by the protein, so no conservation track is drawn and "
            "no SIFT or PROVEAN score is produced.");
    } else {
        // Entropy per column, in bits, on a fixed 0..log2(20) axis so two proteins
        // are comparable and a flat-looking track cannot come from autoscaling.
        std::vector<double> x, y;
        x.reserve(profile.columns.size());
        y.reserve(profile.columns.size());
        for (const auto& c : profile.columns) {
            x.push_back(static_cast<double>(c.position));
            y.push_back(c.shannonEntropy);
        }
        if (ImPlot::BeginPlot("##conservation", ImVec2(-1, 190))) {
            ImPlot::SetupAxes("Query position", "Shannon entropy (bits)");
            ImPlot::SetupAxesLimits(1.0, static_cast<double>(profile.columns.size()), 0.0,
                                    4.3219280948873626, ImPlotCond_Always);
            ImPlot::PlotBars("entropy", x.data(), y.data(), static_cast<int>(x.size()), 0.8);
            ImPlot::EndPlot();
        }
        ImGui::TextColored(theme::provenanceColor(Provenance::Measured),
                           "0 bits = fully conserved in this set; log2(20) = 4.3219 bits = "
                           "uniform. Background for the log-odds: %s",
                           profile.backgroundFrequencySource.c_str());
    }
    ImGui::Spacing();

    theme::sectionHeader("SUBSTITUTION");
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputInt("Position (1-based)", &position);
    position = std::max(1, position);
    char mut[2] = {mutant, '\0'};
    ImGui::SetNextItemWidth(140.0f);
    if (ImGui::InputText("Mutant residue", mut, sizeof(mut)))
        mutant = static_cast<char>(std::toupper(static_cast<unsigned char>(mut[0])));

    const VariantScore score = s.variants->score(profile, position, mutant);
    ImGui::Text("Wild type at %d: %c", score.position, score.wildType);
    drawQuantity("BLOSUM62 delta", score.blosum62Delta);
    drawQuantity("SIFT-style index", score.siftScore);
    ImGui::SameLine();
    ImGui::TextDisabled("(deleterious below %.2f)", score.siftDeleteriousBelow);
    drawQuantity("PROVEAN-style delta", score.proveanScore);
    ImGui::SameLine();
    ImGui::TextDisabled("(deleterious below %.3f)", score.proveanDeleteriousBelow);
    drawQuantity("Column entropy", score.columnEntropy);
    ImGui::TextColored(theme::provenanceColor(Provenance::Measured),
                       "Computed from %d homolog(s), median identity %.1f%%.",
                       score.homologs.sequenceCount, score.homologs.medianIdentityPct);
    ImGui::TextWrapped("%s", score.interpretation.c_str());
    ImGui::Spacing();

    theme::sectionHeader("POINT-MUTATION REBUILD");
    ImGui::SetNextItemWidth(-140.0f);
    ImGui::InputTextWithHint("##varstruct", "Path to a .pdb / .cif file...", structPath,
                             sizeof(structPath));
    ImGui::SameLine();
    if (ImGui::Button("Load structure", ImVec2(120, 0))) {
        structure = s.structure ? s.structure->load(std::filesystem::path(structPath))
                                : std::nullopt;
        structureError = structure ? std::string()
                                   : "Could not read '" + std::string(structPath) + "'.";
    }
    if (!structureError.empty())
        ImGui::TextColored(theme::verdictColor(3), "%s", structureError.c_str());
    if (!structure) {
        ImGui::TextDisabled("Load a structure to rebuild a side chain on it.");
        return;
    }

    static char chainId[8] = "A";
    static int  residueNumber = 1;
    ImGui::SetNextItemWidth(80.0f);
    ImGui::InputText("Chain", chainId, sizeof(chainId));
    ImGui::SameLine();
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputInt("Residue (author numbering)", &residueNumber);

    const RotamerRebuild r = s.variants->rebuild(*structure, chainId, residueNumber, mutant);
    // The badge is the point: a built side chain is a constructed artefact.
    ImGui::TextColored(theme::provenanceColor(r.provenance), "[%s]  %c%d%c",
                       provenanceLabel(r.provenance), r.wildType, r.position, r.mutant);
    if (!r.rotamerLibrarySource.empty())
        ImGui::TextWrapped("Rotamer library: %s", r.rotamerLibrarySource.c_str());
    if (!r.chiAngles.empty()) {
        std::string chis;
        for (std::size_t i = 0; i < r.chiAngles.size(); ++i)
            chis += (i ? ", " : "") + std::string("chi") + std::to_string(i + 1) + " " +
                    f2(r.chiAngles[i]) + " deg";
        ImGui::TextWrapped("%s", chis.c_str());
        ImGui::Text("Library probability %.3f   heavy-atom clashes %d", r.rotamerProbability,
                    r.clashCount);
    }
    if (!r.repackedNeighbours.empty()) {
        std::string list;
        for (const auto& n : r.repackedNeighbours) list += (list.empty() ? "" : ", ") + n;
        ImGui::TextWrapped("Repacked neighbours: %s", list.c_str());
    }
    for (const auto& a : r.assumptions) ImGui::BulletText("%s", a.c_str());
    for (const auto& w : r.warnings)
        ImGui::TextColored(theme::verdictColor(3), "%s", w.c_str());
    ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed),
                       "Stability change (ddG): not computed. No model weights ship in this "
                       "build - see docs/variants.md.");
}

// ------------------------------------------------- Reaction network / systems
// The three Phase 14 panels share one editable NetworkSpec: the Reaction Network
// panel builds or imports it, and the Metabolic Flux and Pathway Enrichment panels
// read the SAME object. Keeping one copy is the whole point of a single
// stoichiometric matrix - a flux computed over a network the user cannot see in the
// Networks panel would be a second, invisible model.
namespace {

struct NetworkWorkspace {
    NetworkSpec  spec;
    std::string  loadError;
    std::string  loadedFrom = "built-in example";
    int          exampleIndex = 0;
    // Deterministic run settings.
    int          methodIndex = 1;   // 0 = rk4, 1 = rosenbrock
    double       horizon = 10.0;
    double       relTol = 1e-8;
    double       absTol = 1e-12;
    // Stochastic run settings.
    int          replicates = 200;
    int          seed = 1;
    bool         tauLeap = false;
    // Chemical kinetics data the USER enters. Empty by default: there is no default
    // degradation rate, and pre-filling one would be the fabrication this phase
    // removed from the Stability panel.
    std::vector<std::pair<double, double>> temperatureRates;   // K, k
    std::vector<std::pair<double, double>> phRates;            // pH, k_obs
    double       storageTemperatureC = 25.0;
    double       fractionLost = 0.10;
};

const char* const kNetworkExamples[] = {"A -> B (first order)", "A <=> B (reversible)",
                                        "A -> B -> C (consecutive)",
                                        "Robertson (stiff, 3 species)",
                                        "enzyme step (Michaelis-Menten)"};

SpeciesSpec makeSpecies(const char* id, double c0, bool boundary = false,
                        const char* formula = "") {
    SpeciesSpec s;
    s.id = id;
    s.name = id;
    s.initialConcentration = c0;
    s.boundary = boundary;
    s.formula = formula;
    s.compartment = "cell";
    return s;
}

ReactionSpec makeReaction(const char* id, std::vector<std::pair<std::string, double>> reactants,
                          std::vector<std::pair<std::string, double>> products, RateLaw law,
                          std::vector<double> parameters) {
    ReactionSpec r;
    r.id = id;
    r.reactants = std::move(reactants);
    r.products = std::move(products);
    r.law = law;
    r.parameters = std::move(parameters);
    r.reversible = law == RateLaw::ReversibleMassAction ||
                   law == RateLaw::ReversibleMichaelisMenten;
    return r;
}

NetworkSpec buildExample(int index) {
    NetworkSpec n;
    switch (index) {
        case 1:
            n.id = "A <=> B";
            n.species = {makeSpecies("A", 1.0), makeSpecies("B", 0.0)};
            n.reactions = {makeReaction("R1", {{"A", 1}}, {{"B", 1}},
                                        RateLaw::ReversibleMassAction, {0.3, 0.7})};
            break;
        case 2:
            n.id = "A -> B -> C";
            n.species = {makeSpecies("A", 1.0), makeSpecies("B", 0.0), makeSpecies("C", 0.0)};
            n.reactions = {makeReaction("R1", {{"A", 1}}, {{"B", 1}}, RateLaw::MassAction, {1.0}),
                           makeReaction("R2", {{"B", 1}}, {{"C", 1}}, RateLaw::MassAction, {0.3})};
            break;
        case 3:
            // The classic stiff test problem. y3 catalyses R2 and y2 catalyses R3,
            // which is what makes y1 + y2 + y3 exactly conserved.
            n.id = "Robertson";
            n.species = {makeSpecies("y1", 1.0), makeSpecies("y2", 0.0), makeSpecies("y3", 0.0)};
            n.reactions = {
                makeReaction("R1", {{"y1", 1}}, {{"y2", 1}}, RateLaw::MassAction, {0.04}),
                makeReaction("R2", {{"y2", 1}, {"y3", 1}}, {{"y1", 1}, {"y3", 1}},
                             RateLaw::MassAction, {1e4}),
                makeReaction("R3", {{"y2", 2}}, {{"y3", 1}, {"y2", 1}}, RateLaw::MassAction,
                             {3e7})};
            break;
        case 4:
            n.id = "enzyme step";
            n.species = {makeSpecies("S", 1.0), makeSpecies("P", 0.0)};
            n.reactions = {makeReaction("E1", {{"S", 1}}, {{"P", 1}}, RateLaw::MichaelisMenten,
                                        {0.5, 0.2})};
            break;
        default:
            n.id = "A -> B";
            n.species = {makeSpecies("A", 1.0), makeSpecies("B", 0.0)};
            n.reactions = {makeReaction("R1", {{"A", 1}}, {{"B", 1}}, RateLaw::MassAction, {0.7})};
            break;
    }
    return n;
}

NetworkWorkspace& networkWorkspace() {
    static NetworkWorkspace ws = [] {
        NetworkWorkspace w;
        w.spec = buildExample(0);
        return w;
    }();
    return ws;
}

const char* rateLawName(RateLaw law) {
    switch (law) {
        case RateLaw::MassAction: return "mass action";
        case RateLaw::ReversibleMassAction: return "reversible mass action";
        case RateLaw::MichaelisMenten: return "Michaelis-Menten";
        case RateLaw::ReversibleMichaelisMenten: return "reversible Michaelis-Menten";
        case RateLaw::Hill: return "Hill";
    }
    return "?";
}

std::string sideText(const std::vector<std::pair<std::string, double>>& side) {
    std::string s;
    for (const auto& [id, stoich] : side) {
        if (!s.empty()) s += " + ";
        if (std::abs(stoich - 1.0) > 1e-9) s += f2(stoich) + " ";
        s += id;
    }
    return s.empty() ? "(none)" : s;
}

// The structural analysis, refreshed whenever the network changes. It is cheap and
// it is what the conserved-quantity audit and the Wegscheider check depend on.
void refreshAnalysis(AppShell& shell, NetworkWorkspace& ws) {
    if (shell.services().simulation) ws.spec = shell.services().simulation->analyze(ws.spec);
}

// Species-to-species edges induced by shared reactions. Every edge carries its
// evidence (the reaction id) and Provenance::Model, because the network is a
// constructed artefact - a written-down mechanism - and not a measurement.
std::vector<NetworkEdge> networkAsEdges(const NetworkSpec& spec) {
    std::vector<NetworkEdge> edges;
    for (const ReactionSpec& r : spec.reactions) {
        for (const auto& [a, sa] : r.reactants) {
            for (const auto& [b, sb] : r.products) {
                if (a == b) continue;
                NetworkEdge e;
                e.source = a;
                e.target = b;
                e.weight = 1.0;
                e.evidence = "reaction '" + r.id + "' (" + rateLawName(r.law) + ")";
                e.provenance = Provenance::Model;
                edges.push_back(std::move(e));
            }
        }
    }
    return edges;
}

}  // namespace

void networks(AppShell& shell) {
    Services& s = shell.services();
    if (!s.simulation) return;
    NetworkWorkspace& ws = networkWorkspace();

    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::Model),
                       "A reaction network is a mechanism you wrote down. Its time course is a "
                       "property of that mechanism and of the rate constants you entered - never "
                       "a measurement of a cell, and never derived from a docking score.");
    ImGui::PopTextWrapPos();
    ImGui::Spacing();

    // ------------------------------------------------------------ the network
    theme::sectionHeader("NETWORK");
    ImGui::SetNextItemWidth(280.0f);
    if (ImGui::Combo("Example", &ws.exampleIndex, kNetworkExamples, 5)) {
        ws.spec = buildExample(ws.exampleIndex);
        ws.loadedFrom = "built-in example";
        ws.loadError.clear();
        refreshAnalysis(shell, ws);
    }
    ImGui::SameLine();
    static char sbmlPath[512] = "";
    ImGui::SetNextItemWidth(320.0f);
    ImGui::InputText("SBML file", sbmlPath, sizeof(sbmlPath));
    ImGui::SameLine();
    if (ImGui::Button("Import SBML")) {
        std::ifstream in(sbmlPath, std::ios::binary);
        if (!in) {
            ws.loadError = std::string("cannot open '") + sbmlPath + "'";
        } else {
            std::ostringstream buffer;
            buffer << in.rdbuf();
            std::string error;
            const auto imported = s.simulation->importSbml(buffer.str(), &error);
            if (imported) {
                ws.spec = *imported;
                ws.loadedFrom = sbmlPath;
                ws.loadError.clear();
                refreshAnalysis(shell, ws);
            } else {
                // The refusal, in full, naming the construct. It is the most useful
                // thing this panel can say about a document it will not open.
                ws.loadError = error;
            }
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Export SBML L3V2")) {
        const std::string xml = s.simulation->exportSbml(ws.spec);
        std::ofstream out(std::string(sbmlPath).empty() ? "network.sbml" : sbmlPath,
                          std::ios::binary);
        out << xml;
    }
    if (!ws.loadError.empty()) {
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextColored(theme::verdictColor(3), "SBML refused: %s", ws.loadError.c_str());
        ImGui::PopTextWrapPos();
    }
    ImGui::TextDisabled("Loaded from: %s   |   %zu species, %zu reactions",
                        ws.loadedFrom.c_str(), ws.spec.species.size(), ws.spec.reactions.size());

    if (ImGui::BeginTable("net-rx", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Reaction", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Stoichiometry", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("Rate law", ImGuiTableColumnFlags_WidthFixed, 200.0f);
        ImGui::TableSetupColumn("Parameters", ImGuiTableColumnFlags_WidthFixed, 260.0f);
        ImGui::TableHeadersRow();
        for (std::size_t i = 0; i < ws.spec.reactions.size(); ++i) {
            ReactionSpec& r = ws.spec.reactions[i];
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextUnformatted(r.id.c_str());
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%s %s %s", sideText(r.reactants).c_str(), r.reversible ? "<=>" : "->",
                        sideText(r.products).c_str());
            ImGui::TableSetColumnIndex(2);
            ImGui::TextUnformatted(rateLawName(r.law));
            ImGui::TableSetColumnIndex(3);
            for (std::size_t p = 0; p < r.parameters.size(); ++p) {
                if (p) ImGui::SameLine();
                ImGui::PushID(static_cast<int>(i * 16 + p));
                ImGui::SetNextItemWidth(80.0f);
                if (ImGui::InputDouble("##p", &r.parameters[p], 0.0, 0.0, "%.4g"))
                    refreshAnalysis(shell, ws);
                ImGui::PopID();
            }
        }
        ImGui::EndTable();
    }
    ImGui::Spacing();

    // ------------------------------------------------- structure of the matrix
    theme::sectionHeader("STRUCTURE OF THE STOICHIOMETRIC MATRIX");
    if (ws.spec.conservationLabels.empty()) {
        ImGui::TextDisabled("No conservation law: every non-clamped species can be created or "
                            "destroyed by this network.");
    } else {
        for (const std::string& label : ws.spec.conservationLabels)
            ImGui::BulletText("conserved moiety: %s", label.c_str());
    }
    ImGui::TextDisabled("%zu thermodynamic (Wegscheider) cycle(s) in the right null space",
                        ws.spec.thermodynamicCycles.size());
    for (const std::string& v : ws.spec.wegscheiderViolations) {
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextColored(theme::verdictColor(3), "WEGSCHEIDER VIOLATION: %s", v.c_str());
        ImGui::PopTextWrapPos();
    }
    for (const std::string& w : ws.spec.warnings)
        ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic), "%s", w.c_str());
    ImGui::Spacing();

    // ------------------------------------------------------- deterministic run
    theme::sectionHeader("TIME COURSE");
    static const char* kMethods[] = {"rk4 (fixed step, explicit)",
                                     "rosenbrock ROS3P (L-stable, adaptive)"};
    ImGui::SetNextItemWidth(300.0f);
    ImGui::Combo("Solver", &ws.methodIndex, kMethods, 2);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(110.0f);
    ImGui::InputDouble("Horizon", &ws.horizon, 1.0, 10.0, "%.3f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(110.0f);
    ImGui::InputDouble("rel tol", &ws.relTol, 0.0, 0.0, "%.1e");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(110.0f);
    ImGui::InputDouble("abs tol", &ws.absTol, 0.0, 0.0, "%.1e");

    const TimeCourse tc = s.simulation->integrate(ws.spec, ws.horizon, ws.relTol, ws.absTol,
                                                  ws.methodIndex == 0 ? "rk4" : "rosenbrock");
    if (tc.times.empty()) {
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextColored(theme::verdictColor(3), "The integration was refused: %s",
                           tc.solver.note.c_str());
        for (const std::string& w : tc.warnings)
            ImGui::TextColored(theme::verdictColor(3), "%s", w.c_str());
        ImGui::PopTextWrapPos();
    } else {
        if (ImPlot::BeginPlot("##net-course", ImVec2(-1, 260))) {
            ImPlot::SetupAxes("time", "concentration");
            for (std::size_t i = 0; i < tc.trajectories.size(); ++i)
                ImPlot::PlotLine(tc.speciesIds[i].c_str(), tc.times.data(),
                                 tc.trajectories[i].data(), static_cast<int>(tc.times.size()));
            ImPlot::EndPlot();
        }
        // The solver report is not diagnostics: a curve whose rejected-step count and
        // clip count are hidden cannot be judged.
        theme::sectionHeader("WHAT THE SOLVER ACTUALLY DID");
        if (ImGui::BeginTable("net-solver", 2, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
            ImGui::TableSetupColumn("Quantity", ImGuiTableColumnFlags_WidthFixed, 240.0f);
            ImGui::TableSetupColumn("Value");
            ImGui::TableHeadersRow();
            auto row = [](const char* name, const std::string& value) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(name);
                ImGui::TableSetColumnIndex(1); ImGui::TextUnformatted(value.c_str());
            };
            row("method", tc.solver.method);
            row("relative / absolute tolerance",
                sci(tc.solver.relativeTolerance) + " / " + sci(tc.solver.absoluteTolerance));
            row("accepted steps", std::to_string(tc.solver.acceptedSteps));
            row("rejected steps", std::to_string(tc.solver.rejectedSteps));
            row("Jacobian evaluations / LU factorizations",
                std::to_string(tc.solver.jacobianEvaluations));
            row("non-negativity clips", std::to_string(tc.solver.nonNegativityClips));
            row("CPU seconds", f2(tc.solver.cpuSeconds));
            row("worst conservation drift", sci(tc.worstConservationDrift));
            row("note", tc.solver.note);
            ImGui::EndTable();
        }
        for (const std::string& w : tc.warnings) {
            ImGui::PushTextWrapPos(0.0f);
            ImGui::TextColored(theme::verdictColor(3), "%s", w.c_str());
            ImGui::PopTextWrapPos();
        }
        if (!tc.conservedQuantities.empty() &&
            ImPlot::BeginPlot("##net-conserved", ImVec2(-1, 150))) {
            ImPlot::SetupAxes("time", "conserved quantity");
            for (std::size_t i = 0; i < tc.conservedQuantities.size(); ++i)
                ImPlot::PlotLine(i < ws.spec.conservationLabels.size()
                                     ? ws.spec.conservationLabels[i].c_str()
                                     : "conserved",
                                 tc.times.data(), tc.conservedQuantities[i].data(),
                                 static_cast<int>(tc.times.size()));
            ImPlot::EndPlot();
        }
    }
    ImGui::Spacing();

    // ---------------------------------------------------------- stochastic run
    theme::sectionHeader("STOCHASTIC ENSEMBLE (COPY NUMBERS, NOT CONCENTRATIONS)");
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputInt("replicates", &ws.replicates);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputInt("seed", &ws.seed);
    ImGui::SameLine();
    ImGui::Checkbox("tau-leaping (approximate)", &ws.tauLeap);
    static bool runStochastic = false;
    ImGui::SameLine();
    if (ImGui::Button(runStochastic ? "Hide ensemble" : "Run ensemble"))
        runStochastic = !runStochastic;
    if (runStochastic) {
        const StochasticEnsemble e = s.simulation->stochastic(
            ws.spec, ws.horizon, ws.replicates, static_cast<std::uint64_t>(ws.seed), ws.tauLeap);
        if (e.mean.empty()) {
            ImGui::PushTextWrapPos(0.0f);
            ImGui::TextColored(theme::verdictColor(3), "Refused: %s", e.solver.note.c_str());
            ImGui::PopTextWrapPos();
        } else {
            // Mean plus and minus one standard deviation, because for the models where
            // the SSA matters the SPREAD is the result and a mean line alone hides it.
            if (ImPlot::BeginPlot("##net-ssa", ImVec2(-1, 240))) {
                ImPlot::SetupAxes("time", "copy number");
                for (std::size_t i = 0; i < e.mean.size(); ++i) {
                    std::vector<double> lo(e.times.size()), hi(e.times.size());
                    for (std::size_t k = 0; k < e.times.size(); ++k) {
                        const double sd = std::sqrt(std::max(0.0, e.variance[i][k]));
                        lo[k] = e.mean[i][k] - sd;
                        hi[k] = e.mean[i][k] + sd;
                    }
                    ImPlot::PlotShaded((e.speciesIds[i] + " +/- 1 SD").c_str(), e.times.data(),
                                       lo.data(), hi.data(), static_cast<int>(e.times.size()));
                    ImPlot::PlotLine(e.speciesIds[i].c_str(), e.times.data(), e.mean[i].data(),
                                     static_cast<int>(e.times.size()));
                }
                ImPlot::EndPlot();
            }
            ImGui::TextDisabled("%s, %d replicates, seed %llu, %lld events, %.3f s. %s",
                                e.solver.method.c_str(), e.replicates,
                                static_cast<unsigned long long>(e.solver.seed),
                                static_cast<long long>(e.solver.acceptedSteps),
                                e.solver.cpuSeconds, e.solver.note.c_str());
            if (ws.tauLeap)
                ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic),
                                   "Tau-leaping is APPROXIMATE: its bias is bounded by its own "
                                   "error-control parameter, not by the replicate count. Use the "
                                   "exact algorithm when the mean itself matters.");
        }
    }
    ImGui::Spacing();

    // -------------------------------------------------------- control analysis
    theme::sectionHeader("METABOLIC CONTROL ANALYSIS");
    static bool runMca = false;
    if (ImGui::Button(runMca ? "Hide control analysis" : "Run control analysis"))
        runMca = !runMca;
    if (runMca) {
        const ControlAnalysis ca = s.simulation->controlAnalysis(ws.spec, ws.horizon * 10.0);
        if (ca.fluxControlCoefficients.empty()) {
            for (const std::string& w : ca.warnings) {
                ImGui::PushTextWrapPos(0.0f);
                ImGui::TextColored(theme::verdictColor(3), "%s", w.c_str());
                ImGui::PopTextWrapPos();
            }
        } else {
            // The two residuals come FIRST: they are theorems, and a large residual
            // means the table below it is not trustworthy.
            statCard("SUMMATION RESIDUAL", sci(ca.summationResidual),
                     "flux CCs must sum to 1, concentration CCs to 0", 300.0f);
            ImGui::SameLine();
            statCard("CONNECTIVITY RESIDUAL", sci(ca.connectivityResidual),
                     "sum_k C^J_k * eps_k,j must be 0", 300.0f);
            if (ImGui::BeginTable("mca", static_cast<int>(ca.reactionIds.size() + 1),
                                  ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
                ImGui::TableSetupColumn("flux \\ perturbed", ImGuiTableColumnFlags_WidthFixed,
                                        150.0f);
                for (const std::string& id : ca.reactionIds) ImGui::TableSetupColumn(id.c_str());
                ImGui::TableHeadersRow();
                for (std::size_t r = 0; r < ca.fluxControlCoefficients.size(); ++r) {
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0);
                    ImGui::TextUnformatted(ca.reactionIds[r].c_str());
                    for (std::size_t k = 0; k < ca.fluxControlCoefficients[r].size(); ++k) {
                        ImGui::TableSetColumnIndex(static_cast<int>(k + 1));
                        ImGui::TextUnformatted(f2(ca.fluxControlCoefficients[r][k]).c_str());
                    }
                }
                ImGui::EndTable();
            }
            for (const std::string& w : ca.warnings)
                ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic), "%s", w.c_str());
        }
    }
    ImGui::Spacing();

    // ------------------------------------------------------- chemical kinetics
    theme::sectionHeader("CHEMICAL KINETICS FROM YOUR OWN DATA");
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed),
                       "This is where a pH window, a temperature ceiling and a shelf life come "
                       "from. They used to be produced from a functional-group flag count and "
                       "are now computed only from rate constants you measured. Three or more "
                       "distinct temperatures are required before any extrapolation is offered.");
    ImGui::PopTextWrapPos();

    if (ImGui::Button("Add temperature point"))
        ws.temperatureRates.emplace_back(298.15 + 10.0 * ws.temperatureRates.size(), 1e-3);
    ImGui::SameLine();
    if (ImGui::Button("Clear temperatures")) ws.temperatureRates.clear();
    if (!ws.temperatureRates.empty() &&
        ImGui::BeginTable("kin-t", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, 40.0f);
        ImGui::TableSetupColumn("temperature (K)", ImGuiTableColumnFlags_WidthFixed, 180.0f);
        ImGui::TableSetupColumn("rate constant k (1/time)");
        ImGui::TableHeadersRow();
        for (std::size_t i = 0; i < ws.temperatureRates.size(); ++i) {
            ImGui::TableNextRow();
            ImGui::PushID(static_cast<int>(1000 + i));
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%zu", i + 1);
            ImGui::TableSetColumnIndex(1);
            ImGui::SetNextItemWidth(-1.0f);
            ImGui::InputDouble("##T", &ws.temperatureRates[i].first, 0.0, 0.0, "%.2f");
            ImGui::TableSetColumnIndex(2);
            ImGui::SetNextItemWidth(-1.0f);
            ImGui::InputDouble("##k", &ws.temperatureRates[i].second, 0.0, 0.0, "%.6g");
            ImGui::PopID();
        }
        ImGui::EndTable();
    }
    if (ws.temperatureRates.size() >= 2) {
        std::vector<double> temps, ks;
        for (const auto& [t, k] : ws.temperatureRates) {
            temps.push_back(t);
            ks.push_back(k);
        }
        const KineticsFit fit = s.simulation->arrhenius(temps, ks);
        drawQuantity("Pre-exponential A", fit.preExponential);
        drawQuantity("Activation energy Ea", fit.activationEnergy);
        drawQuantity("Enthalpy of activation dH*", fit.enthalpyOfActivation);
        drawQuantity("Entropy of activation dS*", fit.entropyOfActivation);
        drawQuantity("k at 25 C", fit.predictedRateAt25C);
        if (fit.extrapolationSupported) {
            ImGui::Text("95%% prediction interval on k(25 C): %s to %s",
                        sci(fit.predictionIntervalLow).c_str(),
                        sci(fit.predictionIntervalHigh).c_str());
            ImGui::SetNextItemWidth(140.0f);
            ImGui::InputDouble("storage (C)", &ws.storageTemperatureC, 1.0, 5.0, "%.1f");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(140.0f);
            ImGui::InputDouble("fraction lost", &ws.fractionLost, 0.01, 0.05, "%.3f");
        } else {
            ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed),
                               "Two temperatures fit a straight line EXACTLY and say nothing "
                               "about how wrong it is: no extrapolation, and no shelf life.");
        }
        ImGui::Text("R-squared %.6f", fit.rSquared);
        for (const std::string& a : fit.assumptions) ImGui::BulletText("assumes: %s", a.c_str());
        for (const std::string& w : fit.warnings)
            ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic), "%s", w.c_str());
        // The dH*/dS* joint ellipse, because the two are strongly correlated and two
        // separate error bars overstate what the experiment determined.
        if (fit.confidenceEllipse.size() > 2) {
            std::vector<double> ex, ey;
            for (const auto& [h, e] : fit.confidenceEllipse) {
                ex.push_back(h);
                ey.push_back(e);
            }
            if (ImPlot::BeginPlot("##kin-ellipse", ImVec2(-1, 200))) {
                ImPlot::SetupAxes("dH* (kJ/mol)", "dS* (J/(mol*K))");
                ImPlot::PlotLine("95% joint confidence region", ex.data(), ey.data(),
                                 static_cast<int>(ex.size()));
                const double hx = fit.enthalpyOfActivation.value;
                const double hy = fit.entropyOfActivation.value;
                ImPlot::PlotScatter("fit", &hx, &hy, 1);
                ImPlot::EndPlot();
            }
        }
    } else {
        ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed),
                           "No Arrhenius fit: needs rate constants at two or more temperatures "
                           "(three for an extrapolation).");
    }

    ImGui::Spacing();
    if (ImGui::Button("Add pH point")) ws.phRates.emplace_back(1.0 + ws.phRates.size(), 1e-4);
    ImGui::SameLine();
    if (ImGui::Button("Clear pH points")) ws.phRates.clear();
    if (!ws.phRates.empty() &&
        ImGui::BeginTable("kin-ph", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, 40.0f);
        ImGui::TableSetupColumn("pH", ImGuiTableColumnFlags_WidthFixed, 180.0f);
        ImGui::TableSetupColumn("observed k (1/time)");
        ImGui::TableHeadersRow();
        for (std::size_t i = 0; i < ws.phRates.size(); ++i) {
            ImGui::TableNextRow();
            ImGui::PushID(static_cast<int>(2000 + i));
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%zu", i + 1);
            ImGui::TableSetColumnIndex(1);
            ImGui::SetNextItemWidth(-1.0f);
            ImGui::InputDouble("##ph", &ws.phRates[i].first, 0.0, 0.0, "%.2f");
            ImGui::TableSetColumnIndex(2);
            ImGui::SetNextItemWidth(-1.0f);
            ImGui::InputDouble("##kobs", &ws.phRates[i].second, 0.0, 0.0, "%.6g");
            ImGui::PopID();
        }
        ImGui::EndTable();
    }
    if (ws.phRates.size() >= 3) {
        std::vector<double> phs, ks;
        for (const auto& [p, k] : ws.phRates) {
            phs.push_back(p);
            ks.push_back(k);
        }
        const PhRateProfile pr = s.simulation->phRate(phs, ks);
        drawQuantity("kH (acid-catalysed)", pr.kAcid);
        drawQuantity("k0 (pH-independent)", pr.kNeutral);
        drawQuantity("kOH (base-catalysed)", pr.kBase);
        drawQuantity("pH of minimum degradation", pr.minimumPh);
        drawQuantity("rate at the minimum", pr.minimumRate);
        ImGui::Text("R-squared %.6f", pr.rSquared);
        for (const std::string& w : pr.warnings)
            ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic), "%s", w.c_str());
        if (ImPlot::BeginPlot("##kin-ph-plot", ImVec2(-1, 200))) {
            ImPlot::SetupAxes("pH", "log10 k_obs");
            std::vector<double> ly;
            for (double k : ks) ly.push_back(std::log10(std::max(k, 1e-300)));
            ImPlot::PlotScatter("measured", phs.data(), ly.data(), static_cast<int>(phs.size()),
                                markerSpec(ImPlotMarker_Circle, 5.0f));
            if (pr.minimumPh.provenance != Provenance::NotComputed) {
                const double x = pr.minimumPh.value;
                ImPlot::PlotInfLines("closed-form minimum", &x, 1);
            }
            ImPlot::EndPlot();
        }
    } else {
        ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed),
                           "No pH-rate profile: the model has three parameters (kH, k0, kOH) and "
                           "needs at least three pH values.");
    }
}

// ---------------------------------------------------------------- Metabolic flux
void flux(AppShell& shell) {
    Services& s = shell.services();
    NetworkWorkspace& ws = networkWorkspace();

    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::Model),
                       "A flux distribution is a property of the stoichiometry, the bounds and "
                       "the objective below - all three of which you chose. A computed growth "
                       "rate is never a measurement of an organism.");
    ImGui::PopTextWrapPos();
    ImGui::Spacing();

    if (!s.flux) {
        // The feature is compiled out, and saying which switch turns it on is more
        // useful than an empty panel.
        ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed),
                           "Constraint-based flux is not in this build. Configure with "
                           "-DBIOCAD_ENABLE_FBA=ON to enable it.");
        return;
    }

    theme::sectionHeader("MASS AND CHARGE BALANCE (THE GATE)");
    const FluxSolution check = s.flux->balance(ws.spec);
    statCard("BALANCE", check.massBalanced ? "PASS" : "FAIL",
             check.massBalanced ? "every reaction conserves elements and charge"
                                : "fba() is refused until this passes",
             320.0f);
    if (!check.unbalancedReactions.empty() &&
        ImGui::BeginTable("flux-bal", 1, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Unbalanced or unchecked reaction");
        ImGui::TableHeadersRow();
        for (const std::string& u : check.unbalancedReactions) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextWrapped("%s", u.c_str());
        }
        ImGui::EndTable();
    }
    ImGui::Spacing();

    theme::sectionHeader("OBJECTIVE AND BOUNDS");
    static char objective[128] = "";
    ImGui::SetNextItemWidth(220.0f);
    ImGui::InputText("Objective reaction", objective, sizeof(objective));
    static std::vector<FluxBound> bounds;
    if (ImGui::Button("Add bound") && !ws.spec.reactions.empty())
        bounds.push_back({ws.spec.reactions.front().id, 0.0, 10.0});
    ImGui::SameLine();
    if (ImGui::Button("Clear bounds")) bounds.clear();
    for (std::size_t i = 0; i < bounds.size(); ++i) {
        ImGui::PushID(static_cast<int>(3000 + i));
        char id[128];
        std::snprintf(id, sizeof(id), "%s", bounds[i].reactionId.c_str());
        ImGui::SetNextItemWidth(140.0f);
        if (ImGui::InputText("reaction", id, sizeof(id))) bounds[i].reactionId = id;
        ImGui::SameLine();
        ImGui::SetNextItemWidth(110.0f);
        ImGui::InputDouble("lower", &bounds[i].lower, 0.0, 0.0, "%.4g");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(110.0f);
        ImGui::InputDouble("upper", &bounds[i].upper, 0.0, 0.0, "%.4g");
        ImGui::PopID();
    }
    ImGui::Spacing();

    static int mode = 0;
    static const char* kModes[] = {"FBA", "parsimonious FBA", "FVA", "single deletions",
                                   "double deletions"};
    ImGui::SetNextItemWidth(240.0f);
    ImGui::Combo("Analysis", &mode, kModes, 5);
    static double fraction = 1.0;
    if (mode == 2) {
        ImGui::SameLine();
        ImGui::SetNextItemWidth(140.0f);
        ImGui::InputDouble("objective fraction", &fraction, 0.05, 0.1, "%.2f");
    }
    if (!ImGui::Button("Solve")) {
        ImGui::TextDisabled("Press Solve to run the LP.");
        return;
    }

    const std::string obj = objective;
    if (mode == 0 || mode == 1) {
        const FluxSolution r = mode == 0 ? s.flux->fba(ws.spec, obj, bounds)
                                         : s.flux->parsimonious(ws.spec, obj, bounds);
        ImGui::Text("status: %s", r.solverStatus.c_str());
        ImGui::Text("objective %s = %.10g", r.objectiveReactionId.c_str(), r.objectiveValue);
        for (const std::string& w : r.warnings) {
            ImGui::PushTextWrapPos(0.0f);
            ImGui::TextColored(theme::verdictColor(3), "%s", w.c_str());
            ImGui::PopTextWrapPos();
        }
        // The bounds are printed BESIDE the fluxes, always: a flux without the medium
        // it was allowed is not interpretable.
        if (!r.fluxes.empty() &&
            ImGui::BeginTable("flux-v", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
            ImGui::TableSetupColumn("Reaction", ImGuiTableColumnFlags_WidthFixed, 120.0f);
            ImGui::TableSetupColumn("Flux", ImGuiTableColumnFlags_WidthFixed, 140.0f);
            ImGui::TableSetupColumn("Lower bound", ImGuiTableColumnFlags_WidthFixed, 140.0f);
            ImGui::TableSetupColumn("Upper bound");
            ImGui::TableHeadersRow();
            for (std::size_t i = 0; i < r.fluxes.size(); ++i) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::TextUnformatted(r.reactionIds[i].c_str());
                ImGui::TableSetColumnIndex(1);
                ImGui::TextUnformatted(f2(r.fluxes[i]).c_str());
                if (i < r.exchangeBounds.size()) {
                    ImGui::TableSetColumnIndex(2);
                    ImGui::TextUnformatted(f2(r.exchangeBounds[i].lower).c_str());
                    ImGui::TableSetColumnIndex(3);
                    ImGui::TextUnformatted(f2(r.exchangeBounds[i].upper).c_str());
                }
            }
            ImGui::EndTable();
        }
        return;
    }

    const std::vector<FluxRange> ranges =
        mode == 2 ? s.flux->fva(ws.spec, obj, bounds, fraction)
                  : s.flux->deletions(ws.spec, obj, bounds, mode == 3 ? 1 : 2);
    if (ranges.empty()) {
        ImGui::TextColored(theme::verdictColor(3),
                           "No result: the LP was infeasible, unbounded, or the balance gate "
                           "refused the network.");
        return;
    }
    if (ImGui::BeginTable("flux-range", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Reaction(s)", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn(mode == 2 ? "Minimum" : "Objective with the deletion",
                                ImGuiTableColumnFlags_WidthFixed, 220.0f);
        ImGui::TableSetupColumn(mode == 2 ? "Maximum" : "Undeleted optimum",
                                ImGuiTableColumnFlags_WidthFixed, 200.0f);
        ImGui::TableHeadersRow();
        for (const FluxRange& r : ranges) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextUnformatted(r.reactionId.c_str());
            ImGui::TableSetColumnIndex(1);
            ImGui::TextUnformatted(f2(r.minimum).c_str());
            ImGui::TableSetColumnIndex(2);
            ImGui::TextUnformatted(f2(r.maximum).c_str());
        }
        ImGui::EndTable();
    }
}

// ------------------------------------------------------------ Pathway enrichment
void enrichment(AppShell& shell) {
    Services& s = shell.services();
    if (!s.enrichment) return;

    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::Measured),
                       "The background is not optional. A hypergeometric p-value is a function "
                       "of the universe the query was drawn from, so paste the identifiers your "
                       "experiment could actually have detected - not 'all genes'.");
    ImGui::PopTextWrapPos();
    ImGui::Spacing();

    static char pack[128] = "reactome-human.gmt";
    static char queryText[4096] = "HK1\nPFKL\nALDOA\nGAPDH\nPGK1\nENO1\nPKM\nLDHA\n";
    static char backgroundText[16384] = "";
    ImGui::SetNextItemWidth(280.0f);
    ImGui::InputText("Gene-set pack", pack, sizeof(pack));
    ImGui::Text("Query identifiers (one per line)");
    ImGui::InputTextMultiline("##enr-query", queryText, sizeof(queryText), ImVec2(-1, 120));
    ImGui::Text("Background identifiers (one per line) - REQUIRED");
    ImGui::InputTextMultiline("##enr-bg", backgroundText, sizeof(backgroundText), ImVec2(-1, 120));

    auto split = [](const char* text) {
        std::vector<std::string> out;
        std::istringstream in(text);
        std::string line;
        while (std::getline(in, line)) {
            while (!line.empty() && (line.back() == '\r' || line.back() == ' ')) line.pop_back();
            if (!line.empty()) out.push_back(line);
        }
        return out;
    };

    static bool run = false;
    if (ImGui::Button("Run enrichment")) run = true;
    if (run) {
        const EnrichmentReport r =
            s.enrichment->enrich(split(queryText), split(backgroundText), pack);
        ImGui::TextDisabled("method: %s   |   database release: %s", r.method.c_str(),
                            r.databaseRelease.c_str());
        for (const std::string& w : r.warnings) {
            ImGui::PushTextWrapPos(0.0f);
            ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic), "%s", w.c_str());
            ImGui::PopTextWrapPos();
        }
        if (ImGui::BeginTable("enr", 6, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
            ImGui::TableSetupColumn("Pathway", ImGuiTableColumnFlags_WidthFixed, 130.0f);
            ImGui::TableSetupColumn("Name", ImGuiTableColumnFlags_WidthStretch);
            ImGui::TableSetupColumn("Overlap", ImGuiTableColumnFlags_WidthFixed, 90.0f);
            ImGui::TableSetupColumn("p", ImGuiTableColumnFlags_WidthFixed, 110.0f);
            ImGui::TableSetupColumn("q (BH)", ImGuiTableColumnFlags_WidthFixed, 110.0f);
            ImGui::TableSetupColumn("Fold", ImGuiTableColumnFlags_WidthFixed, 80.0f);
            ImGui::TableHeadersRow();
            for (const EnrichmentHit& h : r.hits) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(h.pathwayId.c_str());
                ImGui::TableSetColumnIndex(1); ImGui::TextUnformatted(h.pathwayName.c_str());
                ImGui::TableSetColumnIndex(2);
                ImGui::Text("%d / %d", h.inSetAndPathway, h.pathwaySize);
                ImGui::TableSetColumnIndex(3); ImGui::TextUnformatted(sci(h.pValue).c_str());
                // The q-value is the number to read, so it is coloured and the raw p
                // is left plain: a raw p over a thousand tests is not a result.
                ImGui::TableSetColumnIndex(4);
                ImGui::TextColored(theme::provenanceColor(h.qValue <= 0.05 ? Provenance::Measured
                                                                          : Provenance::Heuristic),
                                   "%s", sci(h.qValue).c_str());
                ImGui::TableSetColumnIndex(5); ImGui::TextUnformatted(f2(h.foldEnrichment).c_str());
            }
            ImGui::EndTable();
        }
    }
    ImGui::Spacing();

    // ------------------------------------------------------------ graph metrics
    theme::sectionHeader("TOPOLOGY OF THE LOADED REACTION NETWORK");
    ImGui::TextDisabled("Edges are induced by shared reactions in the Reaction Network panel's "
                        "network; each carries the reaction id as its evidence.");
    const GraphMetrics g = s.enrichment->graph(networkAsEdges(networkWorkspace().spec));
    if (g.nodes.empty()) {
        ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed),
                           "No edges with stated evidence: nothing to measure.");
        return;
    }
    statCard("NODES", std::to_string(g.nodes.size()), "species", 140.0f);
    ImGui::SameLine();
    statCard("EDGES", std::to_string(g.edges.size()), "with evidence", 150.0f);
    ImGui::SameLine();
    statCard("COMPONENTS", std::to_string(g.componentCount), "connected", 160.0f);
    ImGui::SameLine();
    statCard("MODULARITY", f2(g.modularity), "Louvain partition", 180.0f);
    if (ImGui::BeginTable("graph", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Node", ImGuiTableColumnFlags_WidthFixed, 140.0f);
        ImGui::TableSetupColumn("Degree", ImGuiTableColumnFlags_WidthFixed, 100.0f);
        ImGui::TableSetupColumn("Betweenness (Brandes)", ImGuiTableColumnFlags_WidthFixed, 200.0f);
        ImGui::TableSetupColumn("Community");
        ImGui::TableHeadersRow();
        for (std::size_t i = 0; i < g.nodes.size(); ++i) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(g.nodes[i].c_str());
            ImGui::TableSetColumnIndex(1); ImGui::Text("%d", g.degree[i]);
            ImGui::TableSetColumnIndex(2);
            ImGui::TextUnformatted(f2(g.betweenness[i]).c_str());
            ImGui::TableSetColumnIndex(3); ImGui::Text("%d", g.community[i]);
        }
        ImGui::EndTable();
    }
    for (const std::string& w : g.warnings)
        ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic), "%s", w.c_str());
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

    theme::sectionHeader("Jurisdiction");
    ImGui::TextUnformatted(r.jurisdiction.c_str());
    ImGui::Spacing();
    theme::sectionHeader("Substantial Similarity");
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

    theme::sectionHeader("Target");
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
                // of just this receptor (any of the 29 receptor presets, not only headlines).
                ImGui::SameLine();
                const std::string lbl = "Provision " + st.dockTarget;
                if (ImGui::Button(lbl.c_str()) && shell.provisionTarget(st.dockTarget))
                    userProvision = true;
            }
            ImGui::SameLine();
            ImGui::TextDisabled("%s", prov.status().c_str());
        }
        ImGui::TextDisabled("Downloads vina.exe (size-checked) + prepares receptor PDBQTs from RCSB "
                            "under %%APPDATA%%/BioCAD/runtime. Headline = DAT/NET/SERT/TAAR1; any "
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
        d.fromEngine() ? ("Best affinity " + f2(d.bestAffinity()) + " kcal/mol at " + st.dockTarget +
                  " (docked with " + d.engine + ").")
               : ("Estimated affinity " + f2(d.bestAffinity()) + " kcal/mol at " + st.dockTarget +
                  " (" + d.engine + " - structure-descriptor model, not a docked score).");

    theme::sectionHeader("Docking Result");

    // Headline metric cards. The affinity card is deliberately large + high-contrast so
    // the primary number reads at a glance; companions summarise pose count / provenance.
    statCard("BEST AFFINITY", hasResult ? f2(d.bestAffinity()) : "--",
             "kcal/mol  (more negative = stronger)", 300.0f, 120.0f, 2.7f);
    ImGui::SameLine();
    statCard("POSES", hasResult ? std::to_string(d.poses.size()) : "--",
             "ranked binding modes", 150.0f, 120.0f, 2.0f);
    ImGui::SameLine();
    statCard("ENGINE", hasResult ? (d.fromEngine() ? "REAL" : "ESTIMATE") : "--",
             d.fromEngine() ? "docked score" : "descriptor model", 160.0f, 120.0f, 1.6f);
    if (hasResult && d.fromEngine()) {
        ImGui::SameLine();
        statCard("CONFIDENCE", d.converged ? "HIGH" : "MODERATE",
                 d.converged ? "search converged" : "budget reached", 175.0f, 120.0f, 1.6f);
    }

    if (hasResult) {
        // The tier and the engine travel with the number itself, not in a tooltip.
        // A real dock is Provenance::Model (a constructed pose, not a measurement);
        // the descriptor fallback is Heuristic and therefore carries no unit.
        drawQuantity("Best affinity",
                     d.fromEngine()
                         ? makeQuantity(d.bestAffinity(), "kcal/mol", d.affinitySpread,
                                        Provenance::Model, d.engine)
                         : makeQuantity(d.bestAffinity(), "", 0.0, Provenance::Heuristic,
                                        "descriptor estimate - rank ordering only"));
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
            row("Result type", d.fromEngine() ? "real engine dock" : "descriptor estimate (not a docked score)");
            row("Target", st.dockTarget);
            if (const ReceptorTarget* preset = docking::findPreset(st.dockTarget)) {
                row("Receptor PDB", preset->pdb.empty() ? "(box only)" : preset->pdb);
                char box[96];
                std::snprintf(box, sizeof box, "%.0f x %.0f x %.0f A  @ (%.1f, %.1f, %.1f)",
                              preset->box.sx, preset->box.sy, preset->box.sz,
                              preset->box.cx, preset->box.cy, preset->box.cz);
                row("Search box", box);
            }
            if (d.fromEngine()) {
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
        theme::sectionHeader("Affinity By Pose");
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
    theme::sectionHeader("3D Pose / Ligand");
    if (hasResult) {
        ImGui::SetNextItemWidth(160);
        ImGui::SliderInt("Pose##dock", &poseSel, 0, static_cast<int>(d.poses.size()) - 1);
        ImGui::SameLine();
        ImGui::TextDisabled("engine: %s%s", d.engine.c_str(),
                            d.fromEngine() ? "" : "  (descriptor estimate, not a docked score)");
    }
    static std::string dockKey;
    static chem::Conformer dockConf;
    static chem::Conformer dockPocket;  // receptor binding-pocket overlay (real docks only)
    static ViewerUiState dockViewUi;
    // Data key (rebuilds the geometry) DOES include the pose, but the camera key passed
    // to the viewer below does NOT - so flipping between poses re-renders the new pose
    // without snapping the camera back to a re-fit (the poses share a coordinate frame).
    const std::string dataKey = m.id + "|" + st.dockTarget + "|" + std::to_string(poseSel) +
                                (hasResult ? "|r" : "|p") + (d.fromEngine() ? "|e" : "");
    if (dataKey != dockKey) {
        dockKey = dataKey;
        dockConf = chem::Conformer{};
        dockPocket = chem::Conformer{};
        if (poseSel < static_cast<int>(d.poses.size()) && !d.poses[poseSel].ligand.empty())
            dockConf = d.poses[poseSel].ligand;
        else if (auto parsed = chem::parsePerceived(m.smiles))
            dockConf = chem::embed3D(*parsed);
        // Overlay the receptor pocket only for a REAL dock (the pose shares the
        // receptor's coordinate frame); the descriptor estimate is not receptor-aligned.
        if (d.fromEngine() && !dockConf.empty()) {
            if (const ReceptorTarget* preset = docking::findPreset(st.dockTarget)) {
                const auto rp = docking::locatePreparedReceptor(preset->id);
                if (rp.ready) dockPocket = loadReceptorPocket(rp.path, dockConf, 5.5);
            }
        }
    }
    // Camera-fit key: stable across pose changes (molecule + target + render mode only).
    const std::string camKey = m.id + "|" + st.dockTarget + (d.fromEngine() ? "|e" : "|p");
    molViewer3D(shell, dockConf, camKey, dockViewUi, 280.0f,
                dockPocket.empty() ? nullptr : &dockPocket);
    if (!dockPocket.empty())
        ImGui::TextDisabled("Receptor binding pocket shown (muted); toggle with the Receptor checkbox.");
    // Screenshot/automation hook (mirrors BIOCAD_PANEL/TARGET): scroll the pose viewer
    // into view for capture tooling. Harmless in normal use.
    static const bool kScroll3d = std::getenv("BIOCAD_DOCK_SCROLL3D") != nullptr;
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

    theme::sectionHeader("Target");
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
                        "%%APPDATA%%/BioCAD/cache. Analysis only - no synthesis content.");
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
    ImGui::TextDisabled("%zu run(s) - persisted to SQLite under %%APPDATA%%/BioCAD/biocad.db",
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
    (void)shell;
    const auto& report = docking::targetPackReport();

    ImGui::TextWrapped(
        "The compound and target catalog is DATA, not code. Built-in packs ship as "
        "assets/packs/*.json beside the executable; your own packs go in "
        "%%APPDATA%%/BioCAD/packs and override a built-in pack with the same id. Receptor "
        "PDBQTs are prepared on demand into %%APPDATA%%/BioCAD/runtime/receptors.");
    ImGui::Spacing();
    if (ImGui::Button("Reload packs")) docking::reloadTargetPacks();
    ImGui::SameLine();
    ImGui::TextDisabled("re-reads every pack from disk without restarting");
    ImGui::Spacing();

    // Load failures are the whole point of this panel: a pack that silently
    // vanished is indistinguishable from a broken application.
    if (!report.errors.empty()) {
        theme::sectionHeader("LOAD ERRORS");
        ImGui::PushStyleColor(ImGuiCol_Text, theme::verdictColor(3));
        for (const auto& e : report.errors) ImGui::TextWrapped("%s", e.c_str());
        ImGui::PopStyleColor();
        ImGui::Spacing();
    }

    theme::sectionHeader("LOADED PACKS");
    if (report.packs.empty()) {
        ImGui::TextDisabled("No packs loaded - the application has no catalog.");
        return;
    }
    if (ImGui::BeginTable("packs", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                                          ImGuiTableFlags_SizingStretchProp)) {
        ImGui::TableSetupColumn("Pack", ImGuiTableColumnFlags_WidthFixed, 190.0f);
        ImGui::TableSetupColumn("Origin", ImGuiTableColumnFlags_WidthFixed, 80.0f);
        ImGui::TableSetupColumn("Compounds", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Targets", ImGuiTableColumnFlags_WidthFixed, 80.0f);
        ImGui::TableSetupColumn("Source");
        ImGui::TableHeadersRow();
        for (const auto& pk : report.packs) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(pk.id.c_str());
            ImGui::TableSetColumnIndex(1);
            ImGui::TextDisabled("%s", pk.builtin ? "built-in" : "user");
            ImGui::TableSetColumnIndex(2); ImGui::Text("%zu", pk.compounds.size());
            ImGui::TableSetColumnIndex(3); ImGui::Text("%zu", pk.targets.size());
            ImGui::TableSetColumnIndex(4);
            ImGui::TextDisabled("%s", pk.sourcePath.c_str());
        }
        ImGui::EndTable();
    }
    ImGui::Spacing();

    theme::sectionHeader("TARGETS");
    ImGui::TextDisabled(
        "A target without a binding-site box is a coverage gap, listed honestly rather "
        "than filled in with an invented site. Only boxed targets are dockable.");
    if (ImGui::BeginTable("packtargets", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                                                ImGuiTableFlags_SizingStretchProp)) {
        ImGui::TableSetupColumn("Id", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Name");
        ImGui::TableSetupColumn("PDB", ImGuiTableColumnFlags_WidthFixed, 60.0f);
        ImGui::TableSetupColumn("Box", ImGuiTableColumnFlags_WidthFixed, 110.0f);
        ImGui::TableSetupColumn("Tags", ImGuiTableColumnFlags_WidthFixed, 150.0f);
        ImGui::TableHeadersRow();
        for (const auto& pk : report.packs) {
            for (const auto& t : pk.targets) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(t.target.id.c_str());
                ImGui::TableSetColumnIndex(1); ImGui::TextUnformatted(t.target.name.c_str());
                ImGui::TableSetColumnIndex(2);
                ImGui::TextDisabled("%s", t.target.pdb.empty() ? "-" : t.target.pdb.c_str());
                ImGui::TableSetColumnIndex(3);
                ImGui::TextColored(theme::provenanceColor(t.hasBox ? Provenance::Model
                                                                   : Provenance::NotComputed),
                                   "%s", t.hasBox ? "dockable" : "no box");
                ImGui::TableSetColumnIndex(4);
                std::string tags;
                for (const auto& tag : t.panels) {
                    if (!tags.empty()) tags += ", ";
                    tags += tag;
                }
                if (t.headline) tags = tags.empty() ? "headline" : "headline, " + tags;
                ImGui::TextDisabled("%s", tags.c_str());
            }
        }
        ImGui::EndTable();
    }
}

// -------------------------------------------------------------------- Settings
void settings(AppShell& shell) {
    theme::sectionHeader("AI Assistant");
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
    theme::sectionHeader("API Key (encrypted at rest via Windows DPAPI)");
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

    theme::sectionHeader("Behavior");
    bool autop = shell.autopilot();
    if (ImGui::Checkbox("Autopilot - run navigate/highlight tools automatically", &autop))
        shell.setAutopilot(autop);
    ImGui::Spacing();

    theme::sectionHeader("Compute (docking engine)");
    int compute = shell.computeMode();
    bool computeChanged = false;
    computeChanged |= ImGui::RadioButton("Auto", &compute, 0); ImGui::SameLine();
    computeChanged |= ImGui::RadioButton("GPU", &compute, 1); ImGui::SameLine();
    computeChanged |= ImGui::RadioButton("CPU", &compute, 2);
    if (computeChanged) shell.setComputeMode(compute);
#ifdef BIOCAD_HAVE_CUDA
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

    theme::sectionHeader("Storage");
    ImGui::TextWrapped("All state lives under %%APPDATA%%/BioCAD (db, artifacts, runtime, presets, logs).");
    ImGui::Spacing();

    theme::sectionHeader("Runtime (self-provisioned components)");
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

    // Defaults resolved by compound id, not by index: a pack edit reorders the
    // library, and a raw index would silently start comparing something else.
    static int slot[3] = {-1, -1, -1};
    static bool slotInit = false;
    if (!slotInit) {
        slotInit = true;
        const char* wanted[3] = {"caffeine", "theobromine", "acetaminophen"};
        for (int i = 0; i < 3; ++i) {
            for (int k = 0; k < static_cast<int>(lib.size()); ++k) {
                if (lib[k].id == wanted[i]) { slot[i] = k; break; }
            }
            if (slot[i] < 0 && i < static_cast<int>(lib.size())) slot[i] = i;
        }
    }
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
        rowD("Shelf life", [](const Row& r) { return quantityShort(r.st.shelfLife); });
        rowD("Oral F%", [](const Row& r) { return f0(r.ab.bioavailabilityPct) + "%"; });
        rowD("HIA%", [](const Row& r) { return f0(r.ab.hiaPct) + "%"; });
        rowD("logBB", [](const Row& r) { return f2(r.ab.logBB) + (r.ab.cnsPenetrant ? " (BBB+)" : ""); });
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0); ImGui::TextDisabled("ADMET");  // inline table cell label, not a section header
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
// The sketch canvas is a hand-built graph, so it has no parser run behind it: the
// implicit hydrogens have to be assigned explicitly, and then the ONE canonical
// writer produces the string. This used to be 76 lines of bespoke SMILES emission
// living in a panel file, which meant a drawn structure and a typed structure could
// produce different strings for the same molecule - so nothing downstream could
// compare, cache or deduplicate them.
std::string sketchToSmiles(const Sketch& sk) {
    const int n = static_cast<int>(sk.atoms.size());
    if (n == 0) return "";

    chem::Molecule m;
    m.atoms.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) m.atoms[static_cast<std::size_t>(i)].z = sk.atoms[i].z;

    for (const auto& b : sk.bonds) {
        if (b.a < 0 || b.a >= n || b.b < 0 || b.b >= n || b.a == b.b) continue;
        chem::Bond nb;
        nb.a = b.a;
        nb.b = b.b;
        nb.aromatic = b.aromatic;
        nb.order = b.aromatic ? 1.5 : static_cast<double>(b.order);
        const int idx = static_cast<int>(m.bonds.size());
        m.bonds.push_back(nb);
        m.atoms[static_cast<std::size_t>(b.a)].nbr.push_back(b.b);
        m.atoms[static_cast<std::size_t>(b.b)].nbr.push_back(b.a);
        m.atoms[static_cast<std::size_t>(b.a)].bonds.push_back(idx);
        m.atoms[static_cast<std::size_t>(b.b)].bonds.push_back(idx);
    }

    chem::assignImplicitHydrogens(m);
    const chem::RingInfo rings = chem::perceiveRings(m);
    chem::annotateRings(m, rings);
    chem::perceiveAromaticity(m, rings);
    return chem::canonicalSmiles(m);
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

    // Screenshot/automation hook (mirrors BIOCAD_PANEL/TARGET): BIOCAD_ANALOG_DRAW
    // opens straight to the sketcher seeded with a phenethylamine core. Harmless in
    // normal use; only the capture tooling sets it.
    static const bool kAutoDraw = std::getenv("BIOCAD_ANALOG_DRAW") != nullptr;
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
                    theme::sectionHeader("Properties");
                    ImGui::SetNextItemWidth(200); ImGui::SliderFloat("logP", &logP, -3.0f, 5.0f, "%.2f");
                    ImGui::SetNextItemWidth(200); ImGui::SliderFloat("TPSA", &tpsa, 0.0f, 150.0f, "%.0f");
                    ImGui::SetNextItemWidth(200); ImGui::SliderFloat("MW", &mw, 80.0f, 400.0f, "%.0f");
                    ImGui::SetNextItemWidth(200); ImGui::SliderInt("H-bond donors", &hbd, 0, 6);
                    ImGui::SetNextItemWidth(200); ImGui::SliderInt("H-bond acceptors", &hba, 0, 10);
                    ImGui::SetNextItemWidth(200); ImGui::SliderInt("Rotatable bonds", &rot, 0, 12);
                    ImGui::TableSetColumnIndex(1);
                    theme::sectionHeader("Functional Groups");
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
            if (auto pm = chem::parsePerceived(c.smiles)) { c.formula = chem::molecularFormula(*pm); }
            haveGraph = true;
        } else {
            c.smiles = sketch.smiles;
            c.drugClass = "User-drawn structure";
            if (auto pm = chem::parsePerceived(c.smiles); pm && !pm->empty()) {
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
                if (auto pm = chem::parsePerceived(c.smiles)) aeConf = chem::embed3D(*pm);
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
                if (auto pm = chem::parsePerceived(act.smiles)) {
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
        statCard("STABILITY", f0(st.overallScore) + "/100", "rank ordering, no units", 190.0f);
        ImGui::SameLine();
        statCard("ORAL F", f0(ab.bioavailabilityPct) + "%",
                 ab.cnsPenetrant ? "BBB-permeant" : "low BBB partition", 160.0f);
        ImGui::SameLine();
        statCard("HIA", f0(ab.hiaPct) + "%", "intestinal abs.", 150.0f);
        ImGui::SameLine();
        {
            const std::string aeEpCount = std::to_string(ad.endpoints.size()) + " endpoints";
            const ImU32 aeAdmetCol = ImGui::ColorConvertFloat4ToU32(
                theme::verdictColor(static_cast<int>(ad.overall)));
            theme::metricCard("ADMET", verdictLabel(ad.overall), aeEpCount.c_str(),
                              aeAdmetCol, 170.0f, 104.0f);
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
        theme::sectionHeader("Predicted Byproducts / Interactions");
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

// ---------------------------------------------------------------------------
// DNA / RNA workbench.
//
// Two rules govern this panel. First, an off-target count never appears without
// the scope it was counted in: the scope statement is rendered ABOVE the guide
// table, not in a tooltip, because a number a reader can see while its scope is
// hidden is a number that lies. Second, there is nothing here that orders,
// quotes, or sends a sequence anywhere: the only outputs are FASTA and GenBank
// text in a read-only box the user can copy.
// ---------------------------------------------------------------------------
void nucleicAcid(AppShell& shell) {
    Services& s = shell.services();
    if (!s.nucleicAcid) return;

    // pUC19's polylinker region: small, real, and enough to exercise every track.
    static char text[65536] =
        ">demo Synthetic test insert (paste your own FASTA or GenBank here)\n"
        "GAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCACTGGCCGTCGTTTTACAA\n"
        "CGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATA\n"
        "GCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCCTGATGCGGTATTTTCT\n";
    static bool  circular = false;
    static char  enzymeText[512] = "EcoRI, BamHI, HindIII";
    static int   geneticCode = 1;
    static int   minOrfAa = 30;

    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::Measured), "%s", nucleicScopeNote());
    ImGui::PopTextWrapPos();
    ImGui::Spacing();

    theme::sectionHeader("SEQUENCE (FASTA OR GENBANK)");
    ImGui::InputTextMultiline("##nuc-text", text, sizeof(text), ImVec2(-1, 120));
    ImGui::Checkbox("Treat as circular", &circular);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputInt("transl_table", &geneticCode);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputInt("min ORF aa", &minOrfAa);

    const auto parsed = s.nucleicAcid->parse(text);
    if (!parsed) {
        ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed),
                           "That text is neither FASTA nor GenBank, so nothing below can be "
                           "computed.");
        return;
    }
    NucRecord rec = *parsed;
    if (circular) rec.circular = true;
    const int len = static_cast<int>(rec.sequence.size());
    if (len == 0) {
        ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed),
                           "The record parsed but carries no sequence.");
        return;
    }

    statCard("LENGTH", std::to_string(len) + (rec.kind == NucKind::Rna ? " nt" : " bp"),
             rec.circular ? "circular" : "linear", 190.0f);
    ImGui::SameLine();
    statCard("GC", f2(bio::gcPercent(rec.sequence)) + "%", "G+C / unambiguous bases", 190.0f);
    ImGui::SameLine();
    statCard("FEATURES", std::to_string(rec.features.size()), "from the GenBank table", 180.0f);
    ImGui::SameLine();
    statCard("ID", rec.id.empty() ? "<none>" : rec.id, rec.description.c_str(), 240.0f);
    for (const auto& w : rec.warnings)
        ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic), "%s", w.c_str());
    ImGui::Spacing();

    // ------------------------------------------------------------ feature track
    theme::sectionHeader("FEATURE TRACK");
    if (rec.features.empty()) {
        ImGui::TextDisabled("No features. A FASTA input has none by definition.");
    } else {
        const float w = ImGui::GetContentRegionAvail().x - 8.0f;
        const float rowH = 18.0f;
        const int rows = static_cast<int>(std::min<std::size_t>(rec.features.size(), 12));
        const ImVec2 origin = ImGui::GetCursorScreenPos();
        ImDrawList* dl = ImGui::GetWindowDrawList();
        for (int i = 0; i < rows; ++i) {
            const NucFeature& f = rec.features[static_cast<std::size_t>(i)];
            const float y = origin.y + rowH * static_cast<float>(i);
            dl->AddLine(ImVec2(origin.x, y + rowH * 0.5f), ImVec2(origin.x + w, y + rowH * 0.5f),
                        theme::kTextDim, 1.0f);
            for (const auto& part : f.parts) {
                const float x0 = origin.x + w * static_cast<float>(part.first) / static_cast<float>(len);
                const float x1 = origin.x + w * static_cast<float>(part.second) / static_cast<float>(len);
                const ImU32 col = ImGui::ColorConvertFloat4ToU32(
                    theme::provenanceColor(f.strand == Strand::Forward ? Provenance::Measured
                                                                      : Provenance::Model));
                dl->AddRectFilled(ImVec2(x0, y + 3.0f), ImVec2(std::max(x1, x0 + 2.0f), y + rowH - 3.0f),
                                  col, 2.0f);
            }
            std::string label = f.type;
            for (const auto& q : f.qualifiers)
                if (q.first == "gene" || q.first == "product") { label += " " + q.second; break; }
            dl->AddText(ImVec2(origin.x + 4.0f, y + 2.0f), theme::kTextHi, label.c_str());
        }
        ImGui::Dummy(ImVec2(w, rowH * static_cast<float>(rows) + 4.0f));
    }
    ImGui::Spacing();

    // -------------------------------------------------- restriction map and gel
    theme::sectionHeader("RESTRICTION MAP AND GEL SCHEMATIC");
    ImGui::SetNextItemWidth(-1);
    ImGui::InputText("##nuc-enz", enzymeText, sizeof(enzymeText));
    std::vector<std::string> enzymes;
    {
        std::string cur;
        for (const char* p = enzymeText; ; ++p) {
            if (*p == ',' || *p == '\0') {
                std::size_t a = cur.find_first_not_of(" \t");
                std::size_t b = cur.find_last_not_of(" \t");
                if (a != std::string::npos) enzymes.push_back(cur.substr(a, b - a + 1));
                cur.clear();
                if (*p == '\0') break;
            } else {
                cur += *p;
            }
        }
    }
    const RestrictionDigest dig = s.nucleicAcid->digest(rec, enzymes);
    for (const auto& wn : dig.warnings)
        ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic), "%s", wn.c_str());
    if (dig.sites.empty()) {
        ImGui::TextDisabled("No cut sites for those enzymes in this sequence.");
    } else {
        const float w = ImGui::GetContentRegionAvail().x - 8.0f;
        const ImVec2 origin = ImGui::GetCursorScreenPos();
        ImDrawList* dl = ImGui::GetWindowDrawList();
        dl->AddLine(ImVec2(origin.x, origin.y + 20.0f), ImVec2(origin.x + w, origin.y + 20.0f),
                    theme::kTextDim, 2.0f);
        for (const auto& site : dig.sites) {
            const float x = origin.x + w * static_cast<float>(site.position) / static_cast<float>(len);
            dl->AddLine(ImVec2(x, origin.y + 6.0f), ImVec2(x, origin.y + 34.0f),
                        ImGui::ColorConvertFloat4ToU32(theme::provenanceColor(Provenance::Measured)),
                        1.5f);
            dl->AddText(ImVec2(x + 2.0f, origin.y + 34.0f), theme::kTextHi, site.enzyme.c_str());
        }
        ImGui::Dummy(ImVec2(w, 56.0f));

        // Gel schematic: one lane, band position by log10(length), which is what a
        // real agarose gel's mobility approximates. Lengths are labelled, because
        // reading a size off a picture of a gel is exactly the mistake to avoid.
        if (!dig.fragmentLengths.empty()) {
            const float gw = 120.0f, gh = 220.0f;
            const ImVec2 g = ImGui::GetCursorScreenPos();
            dl->AddRectFilled(g, ImVec2(g.x + gw, g.y + gh), IM_COL32(18, 20, 26, 255), 3.0f);
            const double hi = std::log10(static_cast<double>(len));
            const double lo = std::log10(20.0);
            for (int fl : dig.fragmentLengths) {
                const double t = (hi - std::log10(std::max(20.0, static_cast<double>(fl)))) /
                                 std::max(1e-9, hi - lo);
                const float y = g.y + 12.0f + static_cast<float>(t) * (gh - 24.0f);
                dl->AddRectFilled(ImVec2(g.x + 10.0f, y - 2.0f), ImVec2(g.x + gw - 10.0f, y + 2.0f),
                                  theme::kTextHi, 1.0f);
                dl->AddText(ImVec2(g.x + gw + 6.0f, y - 7.0f), theme::kTextDim,
                            (std::to_string(fl) + " bp").c_str());
            }
            ImGui::Dummy(ImVec2(gw + 90.0f, gh + 6.0f));
            int total = 0;
            for (int fl : dig.fragmentLengths) total += fl;
            ImGui::TextColored(theme::provenanceColor(total == len ? Provenance::Measured
                                                                  : Provenance::Heuristic),
                               "%zu fragments, sum %d bp vs sequence %d bp.",
                               dig.fragmentLengths.size(), total, len);
        }
    }
    ImGui::Spacing();

    // ------------------------------------------ six-frame ruler and ORF listing
    theme::sectionHeader("SIX-FRAME TRANSLATION RULER");
    const TranslationResult tr = s.nucleicAcid->translate(rec, geneticCode, minOrfAa);
    for (const auto& wn : tr.warnings)
        ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic), "%s", wn.c_str());
    static int rulerStart = 0;
    ImGui::SetNextItemWidth(-120.0f);
    ImGui::SliderInt("window start (bp)", &rulerStart, 0, std::max(0, len - 60));
    {
        const int a = std::clamp(rulerStart, 0, std::max(0, len - 1));
        const int b = std::min(len, a + 60);
        ImGui::TextColored(theme::provenanceColor(Provenance::Measured), "%5d %s", a + 1,
                           rec.sequence.substr(static_cast<std::size_t>(a),
                                               static_cast<std::size_t>(b - a)).c_str());
        // Each frame's amino acids are spaced three columns apart so a residue sits
        // over its own codon: a ruler whose letters do not line up is decoration.
        for (std::size_t fi = 0; fi < tr.frames.size(); ++fi) {
            const bool reverse = fi >= 3;
            const int frame = static_cast<int>(fi % 3);
            std::string row(static_cast<std::size_t>(b - a), ' ');
            for (int i = a; i + 2 < b; ++i) {
                if (((i - frame) % 3) != 0) continue;
                const int aaIndex = reverse ? ((len - frame) / 3) - 1 - ((i - frame) / 3)
                                            : (i - frame) / 3;
                if (aaIndex < 0 || aaIndex >= static_cast<int>(tr.frames[fi].size())) continue;
                row[static_cast<std::size_t>(i - a)] = tr.frames[fi][static_cast<std::size_t>(aaIndex)];
            }
            ImGui::TextColored(theme::provenanceColor(reverse ? Provenance::Model
                                                             : Provenance::Predicted),
                               "%s%d   %s", reverse ? "-" : "+", frame + 1, row.c_str());
        }
        // The ruler relies on the theme's monospaced default font for alignment.
    }
    ImGui::Spacing();

    theme::sectionHeader("OPEN READING FRAMES");
    if (tr.orfs.empty()) {
        ImGui::TextDisabled("No ORF of at least %d amino acids in table %d.", minOrfAa,
                            tr.geneticCodeId);
    } else if (ImGui::BeginTable("nuc-orfs", 6, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders)) {
        ImGui::TableSetupColumn("strand");
        ImGui::TableSetupColumn("frame");
        ImGui::TableSetupColumn("begin");
        ImGui::TableSetupColumn("end");
        ImGui::TableSetupColumn("aa");
        ImGui::TableSetupColumn("protein (first 60)");
        ImGui::TableHeadersRow();
        for (const auto& o : tr.orfs) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextUnformatted(o.strand == Strand::Forward ? "+" : "-");
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%d", o.frame + 1);
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%d", o.begin + 1);
            ImGui::TableSetColumnIndex(3);
            ImGui::Text("%d", o.end);
            ImGui::TableSetColumnIndex(4);
            ImGui::Text("%zu%s", o.protein.size(), o.stopped ? "" : " (runs off the end)");
            ImGui::TableSetColumnIndex(5);
            ImGui::TextUnformatted(o.protein.substr(0, 60).c_str());
        }
        ImGui::EndTable();
    }
    ImGui::Spacing();

    // ------------------------------------------------ oligo thermodynamics
    theme::sectionHeader("OLIGO THERMODYNAMICS");
    static char  oligoText[256] = "GTAAAACGACGGCCAGT";
    static double naM = 0.05, mgM = 0.0, oligoM = 2.5e-7, dntpM = 0.0;
    ImGui::SetNextItemWidth(-1);
    ImGui::InputText("##nuc-oligo", oligoText, sizeof(oligoText));
    ImGui::SetNextItemWidth(130.0f);
    ImGui::InputDouble("Na+ (M)", &naM, 0.0, 0.0, "%.3f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(130.0f);
    ImGui::InputDouble("Mg2+ (M)", &mgM, 0.0, 0.0, "%.4f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(150.0f);
    ImGui::InputDouble("oligo (M)", &oligoM, 0.0, 0.0, "%.2e");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(150.0f);
    ImGui::InputDouble("dNTP (M)", &dntpM, 0.0, 0.0, "%.2e");
    {
        const OligoThermo t = s.nucleicAcid->oligo(oligoText, naM, mgM, oligoM, dntpM);
        statCard("Tm", f2(t.tm.value) + " C", "nearest-neighbour", 170.0f);
        ImGui::SameLine();
        statCard("GC", f2(t.gcPercent) + "%", "of the oligo", 150.0f);
        ImGui::SameLine();
        statCard("dG37", f2(t.deltaG37.value), "kcal/mol", 170.0f);
        ImGui::SameLine();
        statCard("LENGTH", std::to_string(t.sequence.size()), "nt", 130.0f);
        drawQuantity("Tm", t.tm);
        drawQuantity("dH", t.deltaH);
        drawQuantity("dS", t.deltaS);
        drawQuantity("dG37", t.deltaG37);
        // The salt and oligo concentration are part of the Tm, so they are printed
        // with it: a Tm quoted without them cannot be reproduced at a bench.
        ImGui::TextColored(theme::provenanceColor(Provenance::Measured),
                           "Conditions: Na+ %.3f M, Mg2+ %.4f M, oligo %.2e M, dNTP %.2e M",
                           t.naMolar, t.mgMolar, t.oligoMolar, t.dntpMolar);
        for (const auto& an : t.assumptions) ImGui::BulletText("%s", an.c_str());
        const std::vector<SecondaryStructure> ss = s.nucleicAcid->selfStructures(oligoText, naM);
        if (ss.empty()) {
            ImGui::TextDisabled("No hairpin or self-dimer above the reporting cutoff.");
        } else {
            for (const auto& st : ss) {
                ImGui::TextColored(theme::provenanceColor(st.deltaG37.value <= -6.0
                                                              ? Provenance::Heuristic
                                                              : Provenance::Model),
                                   "%s at %d: dG37 %.2f kcal/mol (1 M Na+ standard state)",
                                   st.kind.c_str(), st.position, st.deltaG37.value);
                if (!st.alignment.empty()) ImGui::TextUnformatted(st.alignment.c_str());
            }
        }
    }
    ImGui::Spacing();

    // ---------------------------------------------------------- primer design
    theme::sectionHeader("PRIMER DESIGN");
    static int   pBegin = 0, pEnd = 0;
    static double pTm = 60.0;
    if (pEnd == 0) pEnd = std::min(len, 500);
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputInt("product begin", &pBegin);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputInt("product end", &pEnd);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputDouble("target Tm (C)", &pTm, 0.0, 0.0, "%.1f");
    const std::vector<PrimerPair> pairs =
        s.nucleicAcid->designPrimers(rec, std::clamp(pBegin, 0, len),
                                     std::clamp(pEnd, 0, len), pTm);
    if (pairs.empty()) {
        ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed),
                           "No pair satisfies the limits for that interval. A pair that violates "
                           "a hairpin, self-dimer, cross-dimer, GC or Tm-difference limit is "
                           "absent rather than ranked low.");
    } else if (ImGui::BeginTable("nuc-primers", 6,
                                 ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders)) {
        ImGui::TableSetupColumn("forward");
        ImGui::TableSetupColumn("reverse");
        ImGui::TableSetupColumn("Tm F / R");
        ImGui::TableSetupColumn("|dTm|");
        ImGui::TableSetupColumn("product");
        ImGui::TableSetupColumn("worst structure dG37");
        ImGui::TableHeadersRow();
        for (const auto& p : pairs) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextUnformatted(p.forwardOligo.sequence.c_str());
            ImGui::TableSetColumnIndex(1);
            ImGui::TextUnformatted(p.reverseOligo.sequence.c_str());
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%.1f / %.1f C", p.forwardOligo.tm.value, p.reverseOligo.tm.value);
            ImGui::TableSetColumnIndex(3);
            ImGui::Text("%.2f C", p.tmDifference);
            ImGui::TableSetColumnIndex(4);
            ImGui::Text("%d-%d (%d bp)", p.productBegin + 1, p.productEnd, p.productLength);
            ImGui::TableSetColumnIndex(5);
            double worst = 0.0;
            for (const auto& l : p.liabilities) worst = std::min(worst, l.deltaG37.value);
            ImGui::Text("%.2f kcal/mol", worst);
        }
        ImGui::EndTable();
        for (const auto& wn : pairs.front().warnings) ImGui::TextWrapped("%s", wn.c_str());
    }
    ImGui::Spacing();

    // ------------------------------------------------------------ codon metrics
    theme::sectionHeader("CODON METRICS");
    static char cdsText[8192] = "ATGGCTAGCAAAGGTGAAGAACTGTTTACCGGTGTTGTTCCGATTCTGGTTGAACTGGATGGTGATGTTAACTAA";
    static char usageId[64] = "ecoli-k12";
    static char forbidden[256] = "GAATTC, GGATCC";
    ImGui::SetNextItemWidth(-1);
    ImGui::InputTextMultiline("##nuc-cds", cdsText, sizeof(cdsText), ImVec2(-1, 60));
    ImGui::SetNextItemWidth(200.0f);
    ImGui::InputText("usage table id", usageId, sizeof(usageId));
    ImGui::SameLine();
    ImGui::SetNextItemWidth(260.0f);
    ImGui::InputText("forbidden sites", forbidden, sizeof(forbidden));
    {
        const CodonMetrics cm = s.nucleicAcid->codonMetrics(cdsText, usageId);
        drawQuantity("CAI", cm.cai);
        drawQuantity("GC", cm.gcPercent);
        drawQuantity("GC3", cm.gc3Percent);
        ImGui::TextColored(theme::provenanceColor(Provenance::Measured), "Relative to: %s",
                           cm.usageTableName.c_str());
        for (const auto& wn : cm.warnings) ImGui::BulletText("%s", wn.c_str());
    }
    if (ImGui::TreeNode("Constraint-based codon optimization (not an expression prediction)")) {
        std::vector<std::string> forb;
        {
            std::string cur;
            for (const char* p = forbidden; ; ++p) {
                if (*p == ',' || *p == '\0') {
                    std::size_t a = cur.find_first_not_of(" \t");
                    std::size_t b = cur.find_last_not_of(" \t");
                    if (a != std::string::npos) forb.push_back(cur.substr(a, b - a + 1));
                    cur.clear();
                    if (*p == '\0') break;
                } else {
                    cur += *p;
                }
            }
        }
        const CodonOptimizationResult opt =
            s.nucleicAcid->optimizeCodons(cdsText, usageId, forb);
        ImGui::TextColored(theme::provenanceColor(opt.translationPreserved ? Provenance::Measured
                                                                          : Provenance::NotComputed),
                           "Translation preserved: %s", opt.translationPreserved ? "yes" : "NO");
        if (!opt.optimized.empty()) {
            ImGui::InputTextMultiline("##nuc-opt", const_cast<char*>(opt.optimized.c_str()),
                                      opt.optimized.size() + 1, ImVec2(-1, 60),
                                      ImGuiInputTextFlags_ReadOnly);
        }
        drawQuantity("CAI before", opt.before.cai);
        drawQuantity("CAI after", opt.after.cai);
        for (const auto& v : opt.remainingViolations)
            ImGui::TextColored(theme::verdictColor(3), "%s", v.c_str());
        for (const auto& an : opt.assumptions) ImGui::TextWrapped("%s", an.c_str());
        ImGui::TreePop();
    }
    ImGui::Spacing();

    // ------------------------------------------------------------- guide search
    theme::sectionHeader("CRISPR GUIDE SEARCH");
    static char pamText[16] = "NGG";
    static char refText[65536] = "";
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputText("PAM (IUPAC)", pamText, sizeof(pamText));
    ImGui::TextUnformatted("Reference to count off-targets in (FASTA; empty = the sequence above)");
    ImGui::InputTextMultiline("##nuc-ref", refText, sizeof(refText), ImVec2(-1, 60));
    NucRecord ref = rec;
    if (refText[0] != '\0') {
        if (const auto p = s.nucleicAcid->parse(refText)) {
            ref = *p;
        } else {
            ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic),
                               "The reference text did not parse; the sequence above was searched "
                               "instead, and the scope statement below says so.");
        }
    }
    const GuideSearchResult gs = s.nucleicAcid->findGuides(rec, ref, pamText);
    // THE SCOPE COMES FIRST, ALWAYS. Every number in the table below is a count
    // inside this reference and nothing more.
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::Measured), "%s",
                       gs.scopeStatement.c_str());
    ImGui::PopTextWrapPos();
    statCard("BASES SEARCHED", std::to_string(gs.basesSearched), gs.referenceName.c_str(), 230.0f);
    ImGui::SameLine();
    statCard("GENOME-WIDE CLAIM", gs.genomeWideClaimPossible ? "not verifiable" : "impossible",
             "see the scope statement", 250.0f);
    ImGui::SameLine();
    statCard("GUIDES", std::to_string(gs.guides.size()), "PAM-adjacent sites in the target",
             240.0f);
    if (gs.guides.empty()) {
        ImGui::TextDisabled("No %s-adjacent protospacer in the target sequence.", pamText);
    } else if (ImGui::BeginTable("nuc-guides", 7,
                                 ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders |
                                     ImGuiTableFlags_ScrollY,
                                 ImVec2(0, 240))) {
        ImGui::TableSetupColumn("protospacer");
        ImGui::TableSetupColumn("PAM");
        ImGui::TableSetupColumn("pos");
        ImGui::TableSetupColumn("strand");
        ImGui::TableSetupColumn("GC");
        ImGui::TableSetupColumn("off-target 0 / 1 / 2 mm (in this reference only)");
        ImGui::TableSetupColumn("notes");
        ImGui::TableHeadersRow();
        for (const auto& g : gs.guides) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextUnformatted(g.protospacer.c_str());
            ImGui::TableSetColumnIndex(1);
            ImGui::TextUnformatted(g.pam.c_str());
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%d", g.position + 1);
            ImGui::TableSetColumnIndex(3);
            ImGui::TextUnformatted(g.strand == Strand::Forward ? "+" : "-");
            ImGui::TableSetColumnIndex(4);
            ImGui::Text("%.0f%%", g.gcPercent);
            ImGui::TableSetColumnIndex(5);
            ImGui::TextColored(theme::provenanceColor(g.exactOffTargets > 0 ? Provenance::Heuristic
                                                                           : Provenance::Measured),
                               "%d / %d / %d", g.exactOffTargets, g.oneMismatchOffTargets,
                               g.twoMismatchOffTargets);
            ImGui::TableSetColumnIndex(6);
            std::string notes;
            for (const auto& wn : g.warnings) notes += (notes.empty() ? "" : " ") + wn;
            ImGui::TextUnformatted(notes.c_str());
        }
        ImGui::EndTable();
    }
    for (const auto& wn : gs.warnings) ImGui::TextWrapped("%s", wn.c_str());
    ImGui::Spacing();

    // ----------------------------------------------------------------- export
    // FASTA and GenBank, in a read-only box, and that is the complete list of
    // export paths this panel has by design.
    theme::sectionHeader("EXPORT (FASTA / GENBANK ONLY)");
    static int which = 0;
    ImGui::RadioButton("FASTA", &which, 0);
    ImGui::SameLine();
    ImGui::RadioButton("GenBank", &which, 1);
    const std::string out =
        which == 0 ? s.nucleicAcid->toFasta(rec) : s.nucleicAcid->toGenBank(rec);
    ImGui::InputTextMultiline("##nuc-export", const_cast<char*>(out.c_str()), out.size() + 1,
                              ImVec2(-1, 140), ImGuiInputTextFlags_ReadOnly);
}

// ------------------------------------------------------------------ Antibody
// The header states the one thing a reader will otherwise get wrong: the closest
// germline set is a similarity result, not a species. Everything else in the panel
// is either exact arithmetic or an explicit refusal.
void antibody(AppShell& shell) {
    Services& s = shell.services();
    if (!s.biologics) return;

    // 1N8Z chain B and chain A (RCSB entry FASTA): a real, published Fab, so the
    // panel is never empty and the example is checkable.
    static char heavy[8192] =
        "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAY"
        "LQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS";
    static char light[8192] =
        "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQP"
        "EDFATYYCQQHYTTPPTFGQGTKVEIK";
    static int  scheme = 0;
    static int  disulfides = 4;
    static int  protease = 0;
    static int  missed = 1;
    static char structPath[1024] = "";
    static char chainsA[64] = "AB";
    static char chainsB[64] = "C";
    static std::optional<bio::Structure> loaded;
    static std::string loadError;
    static bool ranScan = false;
    static AlanineScanReport scan;

    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::Measured),
                       "IMGT is canonical here; Kabat and Chothia are table-driven VIEWS of it. "
                       "The closest germline set is a SEQUENCE-SIMILARITY result and is NEVER a "
                       "species identification - a camelid VHH scores best against a human IGHV3 "
                       "gene. There is no humanness score, no humanization suggestion, no "
                       "immunogenicity prediction, no affinity-maturation proposal, no aggregation "
                       "free energy, no viscosity, no expression titre and no shelf life in this "
                       "panel, and their absence is deliberate.");
    ImGui::PopTextWrapPos();
    ImGui::Spacing();

    theme::sectionHeader("V-DOMAIN SEQUENCES (ONE-LETTER)");
    ImGui::Text("Heavy / VHH");
    ImGui::InputTextMultiline("##ab-h", heavy, sizeof(heavy), ImVec2(-1, 56));
    ImGui::Text("Light (optional)");
    ImGui::InputTextMultiline("##ab-l", light, sizeof(light), ImVec2(-1, 56));

    const char* schemeNames[] = {"IMGT", "Kabat", "Chothia", "Martin", "AHo"};
    const NumberingScheme schemes[] = {NumberingScheme::Imgt, NumberingScheme::Kabat,
                                      NumberingScheme::Chothia, NumberingScheme::Martin,
                                      NumberingScheme::Aho};
    ImGui::SetNextItemWidth(160.0f);
    ImGui::Combo("Numbering scheme", &scheme, schemeNames, IM_ARRAYSIZE(schemeNames));
    ImGui::Spacing();

    // Numbering aligns the sequence against every germline in the pack, which is far
    // too much work to redo 60 times a second. The result depends only on the two text
    // boxes, so it is cached on them.
    struct NumberedCache {
        std::string heavy, light;
        AbDomain    h, l;
    };
    static NumberedCache cache;
    if (cache.heavy != heavy || cache.light != light) {
        cache.heavy = heavy;
        cache.light = light;
        cache.h = s.biologics->number(heavy, NumberingScheme::Imgt);
        cache.l = *light ? s.biologics->number(light, NumberingScheme::Imgt) : AbDomain{};
    }
    const AbDomain& imgtH = cache.h;
    const AbDomain& imgtL = cache.l;
    const AbDomain shownH = s.biologics->convertScheme(imgtH, schemes[scheme]);
    const AbDomain shownL = imgtL.numbered ? s.biologics->convertScheme(imgtL, schemes[scheme])
                                           : imgtL;

    std::vector<std::string> chains;
    if (*heavy) chains.push_back(heavy);
    if (*light) chains.push_back(light);

    // Everything below either touches SASA or walks the formula table, so the whole
    // bundle is recomputed only when one of its inputs actually changes. A panel that
    // recomputes an interface every frame is a panel that cannot be opened.
    struct ResultCache {
        std::string signature;
        std::vector<SequenceLiability> liabH, liabL;
        DevelopabilityReport dev;
        MassLadder           ladder;
        PeptideMap           map;
        InterfaceReport      interface;
    };
    static ResultCache rc;
    const std::string protName = bio::proteaseNames()[static_cast<std::size_t>(
        std::clamp(protease, 0, static_cast<int>(bio::proteaseNames().size()) - 1))];
    const std::string signature = std::string(heavy) + "|" + light + "|" + structPath + "|" +
                                  chainsA + "|" + chainsB + "|" + std::to_string(disulfides) +
                                  "|" + protName + "|" + std::to_string(missed) + "|" +
                                  (loaded ? loaded->id + loaded->source : std::string("-"));
    if (rc.signature != signature) {
        rc.signature = signature;
        const bio::Structure* st = loaded ? &*loaded : nullptr;
        rc.liabH = s.biologics->liabilities(imgtH, st);
        rc.liabL = imgtL.numbered ? s.biologics->liabilities(imgtL, st)
                                  : std::vector<SequenceLiability>{};
        rc.dev = s.biologics->developability(chains, st);
        rc.ladder = s.biologics->massLadder(chains, disulfides);
        rc.map = s.biologics->digest(heavy, protName, missed);
        rc.interface = st ? s.biologics->interfaceOf(*st, chainsA, chainsB) : InterfaceReport{};
    }

    auto chainName = [](AbChainType c) {
        switch (c) {
            case AbChainType::HeavyVh:      return "VH";
            case AbChainType::LightVKappa:  return "VK";
            case AbChainType::LightVLambda: return "VL";
            case AbChainType::Vhh:          return "VHH";
            case AbChainType::TcrBeta:      return "TCR beta (REJECTED)";
            case AbChainType::TcrAlpha:     return "TCR alpha (REJECTED)";
            case AbChainType::Unknown:      break;
        }
        return "unknown";
    };

    auto domainBlock = [&](const char* title, const AbDomain& imgt, const AbDomain& shown,
                          const std::vector<SequenceLiability>& liabs) {
        if (imgt.sequence.empty()) return;
        theme::sectionHeader(title);
        if (!imgt.numbered) {
            ImGui::TextColored(theme::verdictColor(3),
                               "NOT NUMBERED - no positions are reported at all, because plausible "
                               "wrong numbering is worse than none:");
            for (const auto& f : imgt.anchorFailures)
                ImGui::TextWrapped("  - %s", f.c_str());
            ImGui::TextDisabled("closest germline set %s (%.1f bits) - a similarity result, not a "
                                "species", imgt.closestGermlineSet.c_str(), imgt.bestBitScore);
            return;
        }
        statCard("CHAIN", chainName(imgt.chain), "from the germline gene class", 200.0f);
        ImGui::SameLine();
        statCard("CDR1/2/3", std::to_string(imgt.cdrLengths[0]) + "/" +
                                 std::to_string(imgt.cdrLengths[1]) + "/" +
                                 std::to_string(imgt.cdrLengths[2]),
                 "IMGT lengths", 180.0f);
        ImGui::SameLine();
        statCard("GERMLINE", imgt.closestGermlineSet,
                 (std::string("best ") + f2(imgt.bestBitScore) + " bits, runner-up " +
                  imgt.runnerUpGermlineSet + " " + f2(imgt.runnerUpBitScore))
                     .c_str(),
                 320.0f);
        ImGui::TextDisabled("V gene %s, J gene %s. Germline identity is NOT species "
                            "identification.", imgt.vGene.c_str(), imgt.jGene.c_str());
        for (const auto& w : imgt.warnings) ImGui::TextWrapped("%s", w.c_str());
        if (!shown.numbered) {
            ImGui::TextColored(theme::verdictColor(3), "%s view REFUSED:", schemeNames[scheme]);
            for (const auto& w : shown.warnings) ImGui::TextWrapped("  - %s", w.c_str());
        }

        // The liability overlay: same table, extra column, so a flagged position is
        // seen next to its region rather than in a separate list.
        std::map<int, std::string> flagAt;
        std::map<int, bool> exposureKnown;
        for (const auto& f : liabs) {
            std::string& cell = flagAt[f.sequenceIndex];
            if (!cell.empty()) cell += ", ";
            cell += f.motif;
            exposureKnown[f.sequenceIndex] = f.exposureKnown;
        }

        if (ImGui::BeginTable("ab-numbering", 5,
                              ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                                  ImGuiTableFlags_ScrollY,
                              ImVec2(-1, 260))) {
            ImGui::TableSetupScrollFreeze(0, 1);
            ImGui::TableSetupColumn("IMGT", ImGuiTableColumnFlags_WidthFixed, 70.0f);
            ImGui::TableSetupColumn(schemeNames[scheme], ImGuiTableColumnFlags_WidthFixed, 80.0f);
            ImGui::TableSetupColumn("AA", ImGuiTableColumnFlags_WidthFixed, 40.0f);
            ImGui::TableSetupColumn("Region", ImGuiTableColumnFlags_WidthFixed, 80.0f);
            ImGui::TableSetupColumn("Liability motif");
            ImGui::TableHeadersRow();
            for (std::size_t i = 0; i < imgt.residues.size(); ++i) {
                const auto& r = imgt.residues[i];
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                const bool cdr = r.region.rfind("CDR", 0) == 0;
                ImGui::TextColored(cdr ? ImGui::ColorConvertU32ToFloat4(theme::kWarn)
                                       : ImGui::GetStyle().Colors[ImGuiCol_Text],
                                   "%s", bio::positionLabel(r).c_str());
                ImGui::TableSetColumnIndex(1);
                if (shown.numbered && i < shown.residues.size())
                    ImGui::TextUnformatted(bio::positionLabel(shown.residues[i]).c_str());
                else
                    ImGui::TextDisabled("-");
                ImGui::TableSetColumnIndex(2);
                ImGui::Text("%c", r.aminoAcid);
                ImGui::TableSetColumnIndex(3);
                ImGui::TextUnformatted(r.region.c_str());
                ImGui::TableSetColumnIndex(4);
                const auto f = flagAt.find(r.sequenceIndex);
                if (f == flagAt.end()) {
                    ImGui::TextDisabled("-");
                } else {
                    ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(theme::kWarn), "%s%s",
                                       f->second.c_str(),
                                       exposureKnown[r.sequenceIndex] ? "" : " (exposure unknown)");
                }
            }
            ImGui::EndTable();
        }
        ImGui::Spacing();
    };

    domainBlock("HEAVY / VHH DOMAIN", imgtH, shownH, rc.liabH);
    domainBlock("LIGHT DOMAIN", imgtL, shownL, rc.liabL);

    // ---------------------------------------------------------- developability
    theme::sectionHeader("DEVELOPABILITY DESCRIPTORS (SEQUENCE ARITHMETIC)");
    const DevelopabilityReport& dev = rc.dev;
    drawQuantity("Isoelectric point", dev.isoelectricPoint);
    drawQuantity("Net charge at pH 7.4", dev.netChargeAtPh74);
    drawQuantity("eps280", dev.extinctionCoefficient280);
    drawQuantity("GRAVY", dev.gravy);
    drawQuantity("Aliphatic index", dev.aliphaticIndex);
    drawQuantity("Instability index", dev.instabilityIndex);
    drawQuantity("Fv charge symmetry", dev.fvChargeSymmetry);
    drawQuantity("TAP-style PSH", dev.tapPsh);
    drawQuantity("TAP-style PPC", dev.tapPpc);
    drawQuantity("TAP-style PNC", dev.tapPnc);
    drawQuantity("TAP-style SFvCSP", dev.tapSfvcsp);
    ImGui::TextDisabled("Structure origin: %s",
                        dev.structureOrigin.empty() ? "(none stated)"
                                                    : dev.structureOrigin.c_str());
    for (const auto& a : dev.assumptions) ImGui::TextWrapped("%s", a.c_str());
    for (const auto& w : dev.warnings) ImGui::TextWrapped("%s", w.c_str());
    ImGui::Spacing();

    // ------------------------------------------------------------- mass ladder
    theme::sectionHeader("MASS LADDER");
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputInt("Disulfide bonds", &disulfides);
    const MassLadder& ladder = rc.ladder;
    if (ImGui::BeginTable("ab-mass", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Species");
        ImGui::TableSetupColumn("Average (Da)", ImGuiTableColumnFlags_WidthFixed, 130.0f);
        ImGui::TableSetupColumn("Monoisotopic (Da)", ImGuiTableColumnFlags_WidthFixed, 150.0f);
        ImGui::TableSetupColumn("SS", ImGuiTableColumnFlags_WidthFixed, 40.0f);
        ImGui::TableHeadersRow();
        for (const auto& e : ladder.entries) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextUnformatted(e.species.c_str());
            if (ImGui::IsItemHovered() && !e.note.empty()) ImGui::SetTooltip("%s", e.note.c_str());
            ImGui::TableSetColumnIndex(1);
            // The average mass leads, because above ~10 kDa it is the observable one.
            ImGui::TextColored(theme::provenanceColor(e.average.provenance), "%.3f",
                               e.average.value);
            ImGui::TableSetColumnIndex(2);
            ImGui::TextDisabled("%.4f", e.monoisotopic.value);
            ImGui::TableSetColumnIndex(3);
            ImGui::Text("%d", e.disulfides);
        }
        ImGui::EndTable();
    }
    drawQuantity("Resolving power needed for +0.984 vs 1.003", ladder.requiredResolvingPower);
    for (const auto& a : ladder.assumptions) ImGui::TextWrapped("%s", a.c_str());
    ImGui::Spacing();

    // ---------------------------------------------------------- peptide map
    theme::sectionHeader("PEPTIDE MAP");
    const auto proteases = bio::proteaseNames();
    std::vector<const char*> pnames;
    for (const auto& p : proteases) pnames.push_back(p.c_str());
    ImGui::SetNextItemWidth(160.0f);
    ImGui::Combo("Protease", &protease, pnames.data(), static_cast<int>(pnames.size()));
    ImGui::SameLine();
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputInt("Missed cleavages", &missed);
    const PeptideMap& map = rc.map;
    ImGui::Text("%zu peptides, %.1f%% coverage at zero missed cleavages", map.peptides.size(),
                map.coveragePct);
    if (ImGui::BeginTable("ab-peptides", 6,
                          ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                              ImGuiTableFlags_ScrollY,
                          ImVec2(-1, 200))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("Peptide");
        ImGui::TableSetupColumn("Range", ImGuiTableColumnFlags_WidthFixed, 100.0f);
        ImGui::TableSetupColumn("MC", ImGuiTableColumnFlags_WidthFixed, 40.0f);
        ImGui::TableSetupColumn("Mono (Da)", ImGuiTableColumnFlags_WidthFixed, 110.0f);
        ImGui::TableSetupColumn("b1 / y1", ImGuiTableColumnFlags_WidthFixed, 150.0f);
        ImGui::TableSetupColumn("Variable modifications");
        ImGui::TableHeadersRow();
        for (const auto& p : map.peptides) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextUnformatted(p.sequence.c_str());
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%d..%d", p.begin, p.end);
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%d", p.missedCleavages);
            ImGui::TableSetColumnIndex(3);
            ImGui::Text("%.4f", p.monoisotopic.value);
            ImGui::TableSetColumnIndex(4);
            ImGui::Text("%.3f / %.3f", p.bIons.empty() ? 0.0 : p.bIons.front(),
                        p.yIons.empty() ? 0.0 : p.yIons.front());
            ImGui::TableSetColumnIndex(5);
            std::string mods;
            for (const auto& m : p.modifications) mods += (mods.empty() ? "" : "; ") + m;
            if (mods.empty()) ImGui::TextDisabled("-");
            else ImGui::TextUnformatted(mods.c_str());
        }
        ImGui::EndTable();
    }
    for (const auto& w : map.warnings) ImGui::TextWrapped("%s", w.c_str());
    ImGui::Spacing();

    // -------------------------------------------------- interface + ala scan
    // Same load idiom as the Protein Structure panel: a LOCAL file, nothing fetched.
    theme::sectionHeader("INTERFACE (LOCAL .PDB / .CIF, NOTHING IS FETCHED)");
    ImGui::SetNextItemWidth(-140.0f);
    ImGui::InputTextWithHint("##ab-structpath", "Path to a complex .pdb / .cif file...",
                             structPath, sizeof(structPath));
    ImGui::SameLine();
    if (ImGui::Button("Load complex", ImVec2(120, 0))) {
        ranScan = false;
        loaded = s.structure ? s.structure->load(std::filesystem::path(structPath))
                             : std::optional<bio::Structure>{};
        loadError = loaded ? std::string()
                           : "Could not read '" + std::string(structPath) + "'.";
    }
    if (!loadError.empty()) ImGui::TextColored(theme::verdictColor(3), "%s", loadError.c_str());
    if (!loaded) {
        ImGui::TextDisabled("No complex loaded: interface contacts, BSA and the alanine scan are "
                            "NotComputed until one is.");
        return;
    }
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputText("Chains A", chainsA, sizeof(chainsA));
    ImGui::SameLine();
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputText("Chains B", chainsB, sizeof(chainsB));

    const InterfaceReport& ifr = rc.interface;
    drawQuantity("Buried surface area", ifr.buriedSurfaceArea);
    drawQuantity("SASA(A)", ifr.sasaA);
    drawQuantity("SASA(B)", ifr.sasaB);
    drawQuantity("SASA(AB)", ifr.sasaComplex);
    ImGui::TextWrapped("SASA parameters: %s", ifr.sasaParameters.c_str());
    ImGui::Text("Levy partition: %zu support, %zu core, %zu rim", ifr.supportResidues.size(),
                ifr.coreResidues.size(), ifr.rimResidues.size());
    ImGui::Text("Paratope %zu residues (%zu of them in CDRs), epitope %zu residues",
                ifr.paratope.size(), ifr.cdrContacts.size(), ifr.epitope.size());
    for (const auto& w : ifr.warnings) ImGui::TextWrapped("%s", w.c_str());
    if (ImGui::BeginTable("ab-contacts", 5,
                          ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                              ImGuiTableFlags_ScrollY,
                          ImVec2(-1, 200))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("Residue A");
        ImGui::TableSetupColumn("Residue B");
        ImGui::TableSetupColumn("Min dist", ImGuiTableColumnFlags_WidthFixed, 80.0f);
        ImGui::TableSetupColumn("Atom pairs", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Interactions");
        ImGui::TableHeadersRow();
        for (const auto& c : ifr.contacts) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextUnformatted(c.residueA.c_str());
            ImGui::TableSetColumnIndex(1);
            ImGui::TextUnformatted(c.residueB.c_str());
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%.2f A", c.minDistance);
            ImGui::TableSetColumnIndex(3);
            ImGui::Text("%d", c.atomContacts);
            ImGui::TableSetColumnIndex(4);
            std::string k;
            auto push = [&k](bool on, const char* n) { if (on) k += (k.empty() ? "" : ", ") + std::string(n); };
            push(c.hydrogenBond, "H-bond");
            push(c.saltBridge, "salt bridge");
            push(c.hydrophobic, "hydrophobic");
            push(c.piStacking, "pi-pi");
            push(c.cationPi, "cation-pi");
            push(c.disulfide, "disulfide");
            if (k.empty()) ImGui::TextDisabled("van der Waals only");
            else ImGui::TextUnformatted(k.c_str());
        }
        ImGui::EndTable();
    }
    ImGui::Spacing();

    theme::sectionHeader("GEOMETRIC ALANINE SCAN (NOT A DELTA-DELTA-G)");
    if (ImGui::Button("Run scan", ImVec2(140, 0))) {
        scan = s.biologics->alanineScan(*loaded, chainsA, chainsB);
        ranScan = true;
    }
    if (!ranScan) {
        ImGui::TextDisabled("Not run. It re-measures the interface once per truncated side chain, "
                            "so it costs a few seconds.");
        return;
    }
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(theme::kWarn), "%s", scan.disclaimer.c_str());
    ImGui::PopTextWrapPos();
    drawQuantity("Benchmark Spearman", scan.benchmarkSpearman);
    if (ImGui::BeginTable("ab-scan", 6,
                          ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                              ImGuiTableFlags_ScrollY,
                          ImVec2(-1, 220))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("Residue");
        ImGui::TableSetupColumn("Impact (unit-free)", ImGuiTableColumnFlags_WidthFixed, 140.0f);
        ImGui::TableSetupColumn("Lost BSA (A^2)", ImGuiTableColumnFlags_WidthFixed, 120.0f);
        ImGui::TableSetupColumn("Lost contacts", ImGuiTableColumnFlags_WidthFixed, 110.0f);
        ImGui::TableSetupColumn("Lost H-bonds", ImGuiTableColumnFlags_WidthFixed, 110.0f);
        ImGui::TableSetupColumn("Lost salt bridges", ImGuiTableColumnFlags_WidthFixed, 130.0f);
        ImGui::TableHeadersRow();
        for (const auto& p : scan.positions) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextUnformatted(p.residue.c_str());
            ImGui::TableSetColumnIndex(1);
            ImGui::TextColored(theme::provenanceColor(p.impact.provenance), "%.4f",
                               p.impact.value);
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%.1f", p.lostBuriedAreaA2);
            ImGui::TableSetColumnIndex(3);
            ImGui::Text("%d", p.lostContacts);
            ImGui::TableSetColumnIndex(4);
            ImGui::Text("%d", p.lostHydrogenBonds);
            ImGui::TableSetColumnIndex(5);
            ImGui::Text("%d", p.lostSaltBridges);
        }
        ImGui::EndTable();
    }
    for (const auto& a : scan.assumptions) ImGui::TextWrapped("%s", a.c_str());
}

// ---------------------------------------------------------------- Mechanism of action
// Everything here is RETRIEVED. The panel's job is to keep two facts visible: where
// each record came from, and that an empty list is a statement about the query and the
// source rather than about the compound.
namespace {

void warningsBlock(const std::vector<std::string>& warnings) {
    for (const std::string& w : warnings) {
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextColored(theme::provenanceColor(Provenance::Heuristic), "%s", w.c_str());
        ImGui::PopTextWrapPos();
    }
}

// One shared line: this build cannot retrieve anything, and that is not a result.
void offlineBanner(bool online) {
    ImGui::TextColored(theme::provenanceColor(online ? Provenance::Measured
                                                     : Provenance::NotComputed),
                       online ? "Network transport: available (science feature on)."
                              : "Network transport: NOT compiled in. Nothing below was retrieved, "
                                "and an empty result is an absence of retrieval.");
}

}  // namespace

void mechanism(AppShell& shell) {
    Services& s = shell.services();
    if (!s.mechanism) return;
    const Molecule m = shell.currentMolecule();

    offlineBanner(mechanismNetworkAvailable());
    ImGui::TextWrapped("A mechanism is retrieved with its reference. BioCAD never infers one from a "
                       "docking pose or a fingerprint, and the source's free-text mechanism is shown "
                       "verbatim - it is never parsed into a vocabulary.");
    ImGui::Spacing();

    static std::string lastId;
    static MechanismReport report;
    static bool ran = false;
    if (ImGui::Button("Retrieve mechanisms", ImVec2(180, 0)) || (ran && lastId != m.id)) {
        report = s.mechanism->mechanisms(m.id);
        lastId = m.id;
        ran = true;
    }
    ImGui::SameLine();
    ImGui::TextDisabled("compound: %s", m.name.c_str());
    if (!ran) {
        ImGui::TextDisabled("Not retrieved yet. Nothing is queried until you ask.");
        return;
    }

    statCard("RECORDS", std::to_string(report.entries.size()), "retrieved, with references", 200.0f);
    ImGui::SameLine();
    statCard("ATTEMPTED", report.retrievalAttempted ? "yes" : "no",
             report.networkAvailable ? "query formed and sent" : "no transport", 200.0f);
    ImGui::Spacing();
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed), "%s",
                       report.coverageNote.c_str());
    ImGui::PopTextWrapPos();
    warningsBlock(report.warnings);

    if (report.entries.empty()) return;
    if (ImGui::BeginTable("moa", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Target", ImGuiTableColumnFlags_WidthFixed, 200.0f);
        ImGui::TableSetupColumn("Accession", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Action type", ImGuiTableColumnFlags_WidthFixed, 130.0f);
        ImGui::TableSetupColumn("Mechanism (verbatim)", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("Modality", ImGuiTableColumnFlags_WidthFixed, 110.0f);
        ImGui::TableHeadersRow();
        for (const MechanismEntry& e : report.entries) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextColored(theme::provenanceColor(e.provenance), "%s", e.targetName.c_str());
            ImGui::TableSetColumnIndex(1);
            ImGui::TextUnformatted(e.targetAccession.c_str());
            ImGui::TableSetColumnIndex(2);
            ImGui::TextUnformatted(e.actionType.c_str());
            ImGui::TableSetColumnIndex(3);
            ImGui::TextWrapped("%s", e.freeTextMechanism.c_str());
            ImGui::TableSetColumnIndex(4);
            // BioCAD's own axis, and it stays "unknown" until an assay fit supplies it:
            // the source cannot express inhibition modality, so neither do we.
            ImGui::TextDisabled("%s",
                                e.modality == InhibitionModality::Unknown ? "unknown" : "set");
        }
        ImGui::EndTable();
    }
    ImGui::Spacing();
    theme::sectionHeader("REFERENCES (THE RECORD IS ONLY AS GOOD AS THESE)");
    for (const MechanismEntry& e : report.entries) {
        ImGui::TextDisabled("%s  [%s]  %s", e.targetName.c_str(), e.organism.c_str(),
                            e.source.c_str());
        for (const std::string& r : e.references) ImGui::BulletText("%s", r.c_str());
    }
}

// ----------------------------------------------------------------- Off-target panel
// The unscreened count is rendered FIRST and largest, because it is the honest
// headline: a 44-target panel with a handful screened leaves the unknown fraction
// dominant, and a view that leads with the hits inverts that.
void offTargetPanel(AppShell& shell) {
    Services& s = shell.services();
    if (!s.mechanism) return;
    const Molecule m = shell.currentMolecule();

    static char panelId[64] = "safetyscreen44";
    static double hergIc50Nm = 0.0, freeCmaxNm = 0.0;
    static char hergCite[256] = "";
    static PanelScreenReport report;
    static bool ran = false;

    ImGui::TextWrapped("Runs the docking module over every target in the named panel pack. There is "
                       "no composite safety score and no cross-target ranking: receptor "
                       "preparations and box volumes differ per row, so the scores are not on a "
                       "common scale.");
    ImGui::SetNextItemWidth(220.0f);
    ImGui::InputText("Panel id", panelId, sizeof(panelId));
    ImGui::TextDisabled("shipped: safetyscreen44, safetyscreen87, cipa-currents");

    theme::sectionHeader("MEASURED hERG INPUTS (OPTIONAL, USER-SUPPLIED ONLY)");
    ImGui::TextWrapped("The margin is measured IC50 / free Cmax. BioCAD does not predict a hERG "
                       "IC50 and derives no QT or TdP risk, so leaving these empty leaves the "
                       "margin not computed rather than estimated.");
    ImGui::SetNextItemWidth(160.0f);
    ImGui::InputDouble("hERG IC50 (nM)", &hergIc50Nm, 0.0, 0.0, "%.1f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(160.0f);
    ImGui::InputDouble("free Cmax (nM)", &freeCmaxNm, 0.0, 0.0, "%.1f");
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::InputText("IC50 source", hergCite, sizeof(hergCite));

    if (ImGui::Button("Screen panel", ImVec2(160, 0))) {
        // Exactly one implementation of the interface exists, so the measured hERG
        // inputs are handed to it directly rather than widened into the contract.
        if (auto* real = dynamic_cast<RealMechanism*>(s.mechanism)) {
            HergInput in;
            in.measuredIc50Molar = hergIc50Nm > 0.0 ? hergIc50Nm * 1e-9 : -1.0;
            in.freeCmaxMolar = freeCmaxNm > 0.0 ? freeCmaxNm * 1e-9 : -1.0;
            in.citation = hergCite;
            real->setHergInput(in);
        }
        report = s.mechanism->screenPanel(m, panelId);
        ran = true;
    }
    if (!ran) {
        ImGui::TextDisabled("Not run. Every row docks independently, so a full panel takes a while.");
        return;
    }

    statCard("NOT SCREENED", std::to_string(report.unscreened), "the number that matters", 220.0f);
    ImGui::SameLine();
    statCard("SCREENED", std::to_string(report.screened), "real engine poses only", 200.0f);
    ImGui::SameLine();
    statCard("PANEL SIZE", std::to_string(report.panelSize), "as the panel declares it", 200.0f);
    ImGui::Spacing();
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed), "%s",
                       report.coverageStatement.c_str());
    ImGui::PopTextWrapPos();
    drawQuantity("hERG margin (IC50 / free Cmax)", report.hergSafetyMargin);
    if (report.hergSafetyMargin.provenance != Provenance::NotComputed &&
        report.hergSafetyMargin.value < report.hergMarginFlagBelow)
        ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(theme::kWarn),
                          "Margin below %.0f-fold. That is a flag on the ratio of two "
                          "measurements, not a proarrhythmia prediction.",
                          report.hergMarginFlagBelow);
    warningsBlock(report.warnings);

    if (ImGui::BeginTable("panel", 5,
                          ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY,
                          ImVec2(-1, 300))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("Target", ImGuiTableColumnFlags_WidthFixed, 120.0f);
        ImGui::TableSetupColumn("Screened", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Affinity", ImGuiTableColumnFlags_WidthFixed, 150.0f);
        ImGui::TableSetupColumn("Receptor preparation", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("Box", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableHeadersRow();
        for (const PanelTargetResult& r : report.results) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextUnformatted(r.targetId.c_str());
            ImGui::TableSetColumnIndex(1);
            ImGui::TextColored(theme::provenanceColor(r.screened ? Provenance::Model
                                                                : Provenance::NotComputed),
                               "%s", r.screened ? "yes" : "no");
            ImGui::TableSetColumnIndex(2);
            ImGui::TextColored(theme::provenanceColor(r.affinity.provenance), "%s",
                               quantityShort(r.affinity).c_str());
            ImGui::TableSetColumnIndex(3);
            ImGui::TextWrapped("%s", r.screened ? r.receptorPreparation.c_str()
                                                : r.skipReason.c_str());
            ImGui::TableSetColumnIndex(4);
            ImGui::TextWrapped("%s", r.boxDefinition.c_str());
        }
        ImGui::EndTable();
    }
}

// ------------------------------------------------------------------ Pathway context
void pathwayContext(AppShell& shell) {
    Services& s = shell.services();
    if (!s.mechanism) return;

    static char accession[32] = "P29274";
    static PathwayContext ctx;
    static bool ran = false;

    offlineBanner(mechanismNetworkAvailable());
    ImGui::TextWrapped("Reactome pathway membership (CC0) for one UniProt accession, with its event "
                       "ancestors. There is no pathway impact score: no database supports "
                       "propagating a docking score through a pathway graph.");
    ImGui::SetNextItemWidth(160.0f);
    ImGui::InputText("UniProt accession", accession, sizeof(accession));
    if (ImGui::Button("Retrieve pathways", ImVec2(180, 0))) {
        ctx = s.mechanism->pathways(accession);
        ran = true;
    }
    if (!ran) {
        ImGui::TextDisabled("Not retrieved yet.");
        return;
    }
    ImGui::TextDisabled("source: %s", ctx.source.c_str());
    warningsBlock(ctx.warnings);
    for (const PathwayNode& n : ctx.pathways) {
        // Indent by ancestor depth so the flat mapping reads as the hierarchy it is.
        ImGui::Indent(12.0f * static_cast<float>(std::min<std::size_t>(n.ancestorIds.size(), 6)));
        ImGui::BulletText("%s  %s  [%s]", n.stableId.c_str(), n.name.c_str(), n.species.c_str());
        ImGui::SameLine();
        ImGui::TextDisabled("%s", n.url.c_str());
        ImGui::Unindent(12.0f * static_cast<float>(std::min<std::size_t>(n.ancestorIds.size(), 6)));
    }
    ImGui::Spacing();
    theme::sectionHeader("EXTERNAL PATHWAY BROWSERS");
    ImGui::TextWrapped("%s", keggDeepLinkNote());
    ImGui::TextDisabled("%s", keggPathwayUrl("hsa04726").c_str());
}

// -------------------------------------------------------------------- Stack checker
void stackCheck(AppShell& shell) {
    Services& s = shell.services();
    if (!s.mechanism) return;

    static char membersText[4096] = "caffeine\nfluvoxamine\ngrapefruit-juice\nst-johns-wort\n";
    static StackReport report;
    static bool ran = false;

    ImGui::TextWrapped("Enter every drug, supplement, food or habit in the stack, one per line. Each "
                       "result is a MECHANISM WITH A CITATION, never a severity score and never a "
                       "recommendation. A member this build does not know is listed as unknown "
                       "rather than silently ignored.");
    ImGui::InputTextMultiline("##stack-members", membersText, sizeof(membersText), ImVec2(-1, 120));
    if (ImGui::Button("Check stack", ImVec2(150, 0))) {
        std::vector<std::string> ids;
        std::istringstream in(membersText);
        std::string line;
        while (std::getline(in, line)) {
            while (!line.empty() && (line.back() == '\r' || line.back() == ' ')) line.pop_back();
            if (!line.empty()) ids.push_back(line);
        }
        report = s.mechanism->checkStack(ids);
        ran = true;
    }
    if (!ran) return;

    statCard("FLAGS", std::to_string(report.flags.size()), "mechanisms, not severities", 200.0f);
    ImGui::SameLine();
    statCard("UNKNOWN", std::to_string(report.unknownMembers.size()), "entered but not screened",
             210.0f);
    ImGui::Spacing();
    if (!report.unknownMembers.empty()) {
        std::string list;
        for (const std::string& u : report.unknownMembers) list += (list.empty() ? "" : ", ") + u;
        ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(theme::kWarn),
                          "NOT screened (unknown to the pack): %s", list.c_str());
    }
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(theme::provenanceColor(Provenance::NotComputed), "%s",
                       report.coverageNote.c_str());
    ImGui::PopTextWrapPos();
    warningsBlock(report.warnings);

    if (ImGui::BeginTable("stack", 4,
                          ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY,
                          ImVec2(-1, 280))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("Mechanism", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("Direction", ImGuiTableColumnFlags_WidthFixed, 220.0f);
        ImGui::TableSetupColumn("Evidence", ImGuiTableColumnFlags_WidthFixed, 180.0f);
        ImGui::TableSetupColumn("Citation", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableHeadersRow();
        for (const InteractionFlag& f : report.flags) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextColored(theme::provenanceColor(Provenance::Measured), "%s",
                               f.mechanism.c_str());
            ImGui::TableSetColumnIndex(1);
            ImGui::TextWrapped("%s", f.direction.c_str());
            ImGui::TableSetColumnIndex(2);
            ImGui::TextWrapped("%s", f.evidence.c_str());
            ImGui::TableSetColumnIndex(3);
            ImGui::TextWrapped("%s", f.citation.c_str());
        }
        ImGui::EndTable();
    }
    if (!report.flags.empty()) {
        ImGui::Spacing();
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(theme::kWarn), "%s",
                           report.flags.front().boundaryNote.c_str());
        ImGui::PopTextWrapPos();
    }

    ImGui::Spacing();
    theme::sectionHeader("PHARMACOGENOMIC NOTES (CONDITIONAL, CPIC CC0)");
    const Molecule m = shell.currentMolecule();
    const PharmacogenomicReport pgx = s.mechanism->pharmacogenomics(m.id);
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(ImGui::ColorConvertU32ToFloat4(theme::kWarn), "%s",
                       pgx.boundaryStatement.c_str());
    ImGui::PopTextWrapPos();
    if (pgx.notes.empty()) {
        warningsBlock(pgx.warnings);
        return;
    }
    if (ImGui::BeginTable("pgx", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Gene", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Phenotype", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Activity score", ImGuiTableColumnFlags_WidthFixed, 150.0f);
        ImGui::TableSetupColumn("Conditional implication", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableHeadersRow();
        for (const PharmacogenomicNote& n : pgx.notes) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextUnformatted(n.gene.c_str());
            ImGui::TableSetColumnIndex(1);
            ImGui::TextUnformatted(n.phenotype.c_str());
            ImGui::TableSetColumnIndex(2);
            ImGui::TextUnformatted(n.activityScoreBand.c_str());
            ImGui::TableSetColumnIndex(3);
            ImGui::TextWrapped("%s", n.implication.c_str());
        }
        ImGui::EndTable();
    }
    warningsBlock(pgx.warnings);
}

}  // namespace biocad::panels
