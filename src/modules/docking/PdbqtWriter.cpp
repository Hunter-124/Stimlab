#include "modules/docking/PdbqtWriter.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

namespace stimlab::docking {
namespace {

namespace chem = stimlab::chem;

// Pauling electronegativity for the elements StimLab's parser admits (others fall
// back to carbon). Drives the approximate bond-charge smear below.
double electronegativity(int z) {
    switch (z) {
        case 1:  return 2.20;  // H
        case 5:  return 2.04;  // B
        case 6:  return 2.55;  // C
        case 7:  return 3.04;  // N
        case 8:  return 3.44;  // O
        case 9:  return 3.98;  // F
        case 15: return 2.19;  // P
        case 16: return 2.58;  // S
        case 17: return 3.16;  // Cl
        case 35: return 2.96;  // Br
        case 53: return 2.66;  // I
        default: return 2.55;  // carbon-like default
    }
}

// Element symbol, AutoDock-style (two-letter elements keep their case). Used for
// the PDBQT element field and as the base of the atom-type code.
const char* elementSymbol(int z) {
    switch (z) {
        case 1:  return "H";
        case 5:  return "B";
        case 6:  return "C";
        case 7:  return "N";
        case 8:  return "O";
        case 9:  return "F";
        case 15: return "P";
        case 16: return "S";
        case 17: return "Cl";
        case 35: return "Br";
        case 53: return "I";
        default: return "C";
    }
}

// Build, for every position in the conformer, the list of bonded position indices.
std::vector<std::vector<int>> adjacency(const chem::Conformer& conf) {
    std::vector<std::vector<int>> adj(conf.size());
    for (const auto& b : conf.bonds) {
        if (b.first < 0 || b.second < 0) continue;
        if (b.first >= conf.size() || b.second >= conf.size()) continue;
        adj[b.first].push_back(b.second);
        adj[b.second].push_back(b.first);
    }
    return adj;
}

// A hydrogen position is "polar" (kept in a united-atom PDBQT as HD) iff it is
// bonded to N, O or S. Hydrogens on carbon are nonpolar and omitted.
bool isPolarHydrogen(const chem::Conformer& conf, const std::vector<std::vector<int>>& adj,
                     int posIndex) {
    if (conf.z[posIndex] != 1) return false;
    for (int nb : adj[posIndex]) {
        const int zn = conf.z[nb];
        if (zn == 7 || zn == 8 || zn == 16) return true;
    }
    return false;
}

}  // namespace

// Forward declarations of the public typing/charge helpers used by the line emitter.
std::string autodockAtomType(int z, bool aromatic, bool polarH);
double approxPartialCharge(const chem::Molecule& graph, const chem::Conformer& conf, int posIndex);

namespace {

// Format one fixed-width AutoDock PDBQT ATOM line for kept position `i` with the
// given 1-based `serial`. Shared by the rigid and flexible writers so both emit
// byte-identical atom records (only the surrounding tree differs).
std::string pdbqtAtomLine(const chem::Molecule& graph, const chem::Conformer& conf,
                          const std::vector<std::vector<int>>& adj, int i, int serial,
                          const std::string& res) {
    const int z = conf.z[i];
    const bool isH = (z == 1);
    const bool polar = isH && isPolarHydrogen(conf, adj, i);
    bool aromatic = false;
    if (i < conf.heavyCount && i < static_cast<int>(graph.atoms.size()))
        aromatic = graph.atoms[i].aromatic;
    const std::string adType = autodockAtomType(z, aromatic, polar);
    const double charge = approxPartialCharge(graph, conf, i);
    const chem::Vec3 p = conf.pos[i];
    const double x = std::isfinite(p.x) ? p.x : 0.0;
    const double y = std::isfinite(p.y) ? p.y : 0.0;
    const double zc = std::isfinite(p.z) ? p.z : 0.0;
    char atomName[8];
    std::snprintf(atomName, sizeof(atomName), "%.2s%d", elementSymbol(z), serial);
    char line[128];
    std::snprintf(line, sizeof(line),
                  "ATOM  %5d %-4.4s %-3.3s A%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%10.3f %-2.2s\n",
                  serial, atomName, res.c_str(), 1, x, y, zc, 1.00, 0.00, charge, adType.c_str());
    return line;
}

}  // namespace

std::string autodockAtomType(int z, bool aromatic, bool polarH) {
    switch (z) {
        case 1:  return polarH ? "HD" : "H";   // HD = polar (H-bonding) hydrogen
        case 6:  return aromatic ? "A" : "C";  // A = aromatic carbon
        case 7:  return aromatic ? "NA" : "N"; // NA carries the H-bond-acceptor flag
        case 8:  return "OA";                  // OA = H-bonding oxygen
        case 16: return aromatic ? "SA" : "S"; // SA = H-bonding sulfur
        case 9:  return "F";
        case 15: return "P";
        case 17: return "Cl";
        case 35: return "Br";
        case 53: return "I";
        case 5:  return "B";
        default: return "C";
    }
}

double approxPartialCharge(const chem::Molecule& graph, const chem::Conformer& conf,
                           int posIndex) {
    // Approximate partial charge: a small per-element formal/seed offset plus a
    // bond-by-bond electronegativity smear (the more electronegative partner pulls
    // negative charge). This is intentionally a single-pass approximation, NOT a
    // converged Gasteiger PEOE solve - Vina ranks tolerantly on approximate charges.
    const auto adj = adjacency(conf);
    const int z = conf.z[posIndex];
    double q = 0.0;

    // Seed from any formal charge carried on the heavy-atom graph (positions
    // [0, heavyCount) map 1:1 to graph.atoms; hydrogens have no graph entry).
    if (posIndex < conf.heavyCount && posIndex < static_cast<int>(graph.atoms.size())) {
        q += static_cast<double>(graph.atoms[posIndex].charge);
    }

    const double enSelf = electronegativity(z);
    for (int nb : adj[posIndex]) {
        const double enOther = electronegativity(conf.z[nb]);
        // Each bond shifts ~0.07 e per unit electronegativity difference toward the
        // more electronegative atom; capped so pathological graphs stay finite.
        double dq = 0.07 * (enOther - enSelf);
        if (dq > 0.30) dq = 0.30;
        if (dq < -0.30) dq = -0.30;
        q += dq;
    }
    // Clamp to a physically reasonable window and guard against non-finite input.
    if (!std::isfinite(q)) q = 0.0;
    if (q > 1.5) q = 1.5;
    if (q < -1.5) q = -1.5;
    return q;
}

PdbqtLigand writeRigidPdbqt(const chem::Molecule& graph, const chem::Conformer& conf,
                            const std::string& resName) {
    PdbqtLigand out{};
    if (conf.empty()) return out;

    const auto adj = adjacency(conf);
    std::string body;
    body.reserve(conf.size() * 80);

    // Residue label: PDB columns are 3 wide; truncate/pad defensively.
    const std::string res = resName.empty() ? std::string("LIG") : resName.substr(0, 3);

    int serial = 0;
    for (int i = 0; i < conf.size(); ++i) {
        const int z = conf.z[i];
        const bool isH = (z == 1);
        // United-atom PDBQT: keep heavy atoms and polar H only.
        if (isH && !isPolarHydrogen(conf, adj, i)) continue;
        ++serial;
        body += pdbqtAtomLine(graph, conf, adj, i, serial, res);
        out.serialToConf.push_back(i);
        ++out.atomCount;
        if (isH) ++out.polarH; else ++out.heavyCount;
    }

    // Assemble the rigid block. A single ROOT/ENDROOT with TORSDOF 0 is what keeps
    // Vina's tree parser happy (no flexible torsion tree to mis-read).
    std::string text;
    text.reserve(body.size() + 128);
    text += "REMARK  StimLab rigid ligand (approximate partial charges)\n";
    text += "REMARK  binding-affinity prediction only; not a synthesis artifact\n";
    text += "ROOT\n";
    text += body;
    text += "ENDROOT\n";
    text += "TORSDOF 0\n";
    out.text = std::move(text);
    out.torsions = 0;
    return out;
}

PdbqtLigand writeFlexiblePdbqt(const chem::Molecule& graph, const chem::Conformer& conf,
                               const std::string& resName) {
    PdbqtLigand out{};
    if (conf.empty()) return out;

    const auto adj = adjacency(conf);
    const std::string res = resName.empty() ? std::string("LIG") : resName.substr(0, 3);
    const int n = conf.size();

    // Kept atoms = heavy atoms + polar hydrogens (united-atom convention).
    std::vector<int> keptIndex(n, -1);
    std::vector<int> keptList;
    keptList.reserve(n);
    for (int i = 0; i < n; ++i) {
        const bool isH = (conf.z[i] == 1);
        if (isH && !isPolarHydrogen(conf, adj, i)) continue;
        keptIndex[i] = static_cast<int>(keptList.size());
        keptList.push_back(i);
    }
    if (keptList.empty()) return out;
    const int K = static_cast<int>(keptList.size());

    // ---- rotatable-bond perception (read from the heavy-atom graph) -----------
    auto heavyDeg = [&](int i) {
        return (i < static_cast<int>(graph.atoms.size())) ? static_cast<int>(graph.atoms[i].nbr.size()) : 0;
    };
    auto isCarbonyl = [&](int c) {
        if (c < 0 || c >= static_cast<int>(graph.atoms.size())) return false;
        for (int bi : graph.atoms[c].bonds) {
            const auto& bb = graph.bonds[bi];
            if (bb.order == 2.0 && graph.atoms[bb.other(c)].z == 8) return true;
        }
        return false;
    };
    auto graphBond = [&](int a, int b) -> const chem::Bond* {
        if (a < 0 || a >= static_cast<int>(graph.atoms.size())) return nullptr;
        for (int bi : graph.atoms[a].bonds) {
            const auto& bb = graph.bonds[bi];
            if (bb.other(a) == b) return &bb;
        }
        return nullptr;
    };
    auto isRotatable = [&](int a, int b) -> bool {
        if (a >= conf.heavyCount || b >= conf.heavyCount) return false;  // heavy-heavy only
        const chem::Bond* bb = graphBond(a, b);
        if (!bb) return false;
        if (bb->order != 1.0 || bb->aromatic || bb->inRing) return false;   // single, acyclic
        if (heavyDeg(a) < 2 || heavyDeg(b) < 2) return false;               // non-terminal
        // Amide C-N is planar (partial double-bond character): not a free torsion.
        if ((graph.atoms[a].z == 6 && graph.atoms[b].z == 7 && isCarbonyl(a)) ||
            (graph.atoms[b].z == 6 && graph.atoms[a].z == 7 && isCarbonyl(b)))
            return false;
        return true;
    };

    // ---- fragment the ligand on rotatable bonds (union-find over kept atoms) ---
    std::vector<int> uf(K);
    for (int i = 0; i < K; ++i) uf[i] = i;
    std::function<int(int)> findRoot = [&](int x) {
        while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; }
        return x;
    };
    auto unite = [&](int a, int b) { uf[findRoot(a)] = findRoot(b); };

    struct RotEdge { int a, b; };  // conformer heavy-atom indices
    std::vector<RotEdge> rotEdges;
    for (const auto& bp : conf.bonds) {
        const int u = bp.first, v = bp.second;
        if (u < 0 || v < 0 || u >= n || v >= n) continue;
        if (keptIndex[u] < 0 || keptIndex[v] < 0) continue;  // bond to a dropped nonpolar H
        if (isRotatable(u, v)) { rotEdges.push_back({u, v}); continue; }  // split: do not unite
        unite(keptIndex[u], keptIndex[v]);
    }
    const int nRot = static_cast<int>(rotEdges.size());

    std::unordered_map<int, std::vector<int>> frags;  // fragment root -> kept conf indices
    for (int ki = 0; ki < K; ++ki) frags[findRoot(ki)].push_back(keptList[ki]);
    int rootFrag = findRoot(0);
    size_t bestSize = 0;
    for (const auto& kv : frags)
        if (kv.second.size() > bestSize) { bestSize = kv.second.size(); rootFrag = kv.first; }

    auto fragOf = [&](int confIdx) { return findRoot(keptIndex[confIdx]); };

    // ---- serialise the AutoDock torsion tree ----------------------------------
    std::vector<int> serial(n, -1);
    int next = 1;
    std::string text;
    text.reserve(K * 80 + 128);

    auto writeAtom = [&](int confIdx) {
        serial[confIdx] = next++;
        text += pdbqtAtomLine(graph, conf, adj, confIdx, serial[confIdx], res);
        out.serialToConf.push_back(confIdx);
        if (conf.z[confIdx] == 1) ++out.polarH; else ++out.heavyCount;
        ++out.atomCount;
    };
    auto writeAtoms = [&](int frag, int entry) {
        if (entry >= 0) writeAtom(entry);
        for (int ci : frags[frag])
            if (ci != entry) writeAtom(ci);
    };
    std::function<void(int, int)> writeBranches = [&](int frag, int parentFrag) {
        for (const auto& e : rotEdges) {
            const int fa = fragOf(e.a), fb = fragOf(e.b);
            int a = -1, b = -1, child = -1;
            if (fa == frag && fb != frag && fb != parentFrag) { a = e.a; b = e.b; child = fb; }
            else if (fb == frag && fa != frag && fa != parentFrag) { a = e.b; b = e.a; child = fa; }
            else continue;
            // The child's entry atom `b` will be the very next serial assigned.
            text += "BRANCH " + std::to_string(serial[a]) + " " + std::to_string(next) + "\n";
            writeAtoms(child, b);
            writeBranches(child, frag);
            text += "ENDBRANCH " + std::to_string(serial[a]) + " " + std::to_string(serial[b]) + "\n";
        }
    };

    std::string header;
    header += "REMARK  StimLab flexible ligand (" + std::to_string(nRot) +
              " active torsions; approximate partial charges)\n";
    header += "REMARK  binding-affinity prediction only; not a synthesis artifact\n";
    text = header + "ROOT\n";
    writeAtoms(rootFrag, -1);
    text += "ENDROOT\n";
    writeBranches(rootFrag, -1);
    text += "TORSDOF " + std::to_string(nRot) + "\n";

    out.text = std::move(text);
    out.torsions = nRot;
    return out;
}

}  // namespace stimlab::docking
