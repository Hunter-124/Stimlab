#include "modules/docking/PdbqtWriter.h"

#include <array>
#include <cmath>
#include <cstdio>
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
    out.atomCount = 0;
    out.heavyCount = 0;
    out.polarH = 0;
    if (conf.empty()) return out;

    const auto adj = adjacency(conf);
    std::string body;
    body.reserve(conf.size() * 80);

    // Residue label: PDB columns are 3 wide; truncate/pad defensively.
    std::string res = resName.empty() ? std::string("LIG") : resName.substr(0, 3);

    int serial = 0;
    char line[128];
    for (int i = 0; i < conf.size(); ++i) {
        const int z = conf.z[i];
        const bool isH = (z == 1);
        const bool polar = isH && isPolarHydrogen(conf, adj, i);
        // United-atom PDBQT: keep heavy atoms and polar H only.
        if (isH && !polar) continue;

        bool aromatic = false;
        if (i < conf.heavyCount && i < static_cast<int>(graph.atoms.size()))
            aromatic = graph.atoms[i].aromatic;

        const std::string adType = autodockAtomType(z, aromatic, polar);
        const double charge = approxPartialCharge(graph, conf, i);
        const chem::Vec3 p = conf.pos[i];
        const double x = std::isfinite(p.x) ? p.x : 0.0;
        const double y = std::isfinite(p.y) ? p.y : 0.0;
        const double zc = std::isfinite(p.z) ? p.z : 0.0;

        ++serial;
        // PDB ATOM record name: element symbol + position-local index (cosmetic).
        char atomName[8];
        std::snprintf(atomName, sizeof(atomName), "%.2s%d", elementSymbol(z), serial);

        // Fixed-width AutoDock PDBQT ATOM line. Columns mirror the PDB spec; the
        // trailing partial charge (%6.3f) and right-justified 2-char type are the
        // PDBQT-specific extension fields Vina reads.
        std::snprintf(line, sizeof(line),
                      "ATOM  %5d %-4.4s %-3.3s A%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%10.3f %-2.2s\n",
                      serial, atomName, res.c_str(), 1, x, y, zc, 1.00, 0.00, charge,
                      adType.c_str());
        body += line;

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
    return out;
}

}  // namespace stimlab::docking
