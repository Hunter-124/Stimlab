#include "chem/Analysis.h"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace biocad::chem {

namespace {

bool dbondTo(const Molecule& m, int atom, int z) {
    for (int bi : m.atoms[atom].bonds) {
        const Bond& b = m.bonds[bi];
        if (b.order == 2.0 && m.atoms[b.other(atom)].z == z) return true;
    }
    return false;
}

bool isCarbonyl(const Molecule& m, int c) {
    return m.atoms[c].z == 6 && dbondTo(m, c, 8);
}

bool aromaticNeighbor(const Molecule& m, int atom) {
    for (int nb : m.atoms[atom].nbr)
        if (m.atoms[nb].aromatic) return true;
    return false;
}

// Does aromatic carbon `c` carry a hydroxyl (-OH) substituent?
bool ringCarbonHasOH(const Molecule& m, int c) {
    if (!(m.atoms[c].z == 6 && m.atoms[c].aromatic)) return false;
    for (int nb : m.atoms[c].nbr) {
        const Atom& o = m.atoms[nb];
        if (o.z == 8 && !o.aromatic && o.degree() == 1 && o.totalH() >= 1) return true;
    }
    return false;
}

std::uint64_t mix(std::uint64_t h, std::uint64_t v) {
    h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return h;
}

}  // namespace

FunctionalGroups detectGroups(const Molecule& m) {
    FunctionalGroups g;
    const int n = static_cast<int>(m.atoms.size());

    for (int i = 0; i < n; ++i) {
        const Atom& a = m.atoms[i];
        if (a.aromatic) g.aromaticRing = true;
        if (a.z == 9 || a.z == 17 || a.z == 35 || a.z == 53) g.halogen = true;

        // Carbonyl-centered groups.
        if (isCarbonyl(m, i)) {
            bool singleO = false, singleOwithH = false, toN = false;
            int arylC = 0, alkylC = 0, hOnC = a.totalH();
            for (int nb : a.nbr) {
                const Atom& x = m.atoms[nb];
                // skip the =O itself (double bond)
                bool isCarbonylO = false;
                for (int bi : a.bonds)
                    if (m.bonds[bi].other(i) == nb && m.bonds[bi].order == 2.0) isCarbonylO = true;
                if (isCarbonylO) continue;
                if (x.z == 8) { singleO = true; if (x.totalH() >= 1) singleOwithH = true; }
                else if (x.z == 7) toN = true;
                else if (x.z == 6) { if (x.aromatic) ++arylC; else ++alkylC; }
            }
            if (singleOwithH) g.carboxylicAcid = true;
            else if (singleO) g.ester = true;
            if (toN) g.amide = true;
            if (!singleO && !toN) {
                if (arylC >= 1) g.arylKetone = true;
                if (arylC + alkylC >= 2) g.ketone = true;
                if (hOnC >= 1 && (arylC + alkylC) == 1) g.aldehyde = true;
            }
        }

        // Ether / methylenedioxy (O with two single C neighbors).
        if (a.z == 8 && !a.aromatic && a.degree() == 2 && a.totalH() == 0) {
            bool carbonylAdjacent = false;
            for (int nb : a.nbr) if (isCarbonyl(m, nb)) carbonylAdjacent = true;
            if (!carbonylAdjacent) g.ether = true;
        }

        // Nitrogen classes.
        if (a.z == 7 && !a.aromatic) {
            bool amideN = false;
            for (int nb : a.nbr) if (isCarbonyl(m, nb)) amideN = true;
            if (dbondTo(m, i, 8) && a.charge >= 0 && a.degree() >= 2) g.nitro =
                (a.degree() == 3);  // crude N(=O)-O
            if (!amideN) {
                const int h = a.totalH();
                if (a.degree() == 1 && h >= 2) g.primaryAmine = true;
                else if (a.degree() == 2 && h >= 1) g.secondaryAmine = true;
                else if (a.degree() == 3 && h == 0) g.tertiaryAmine = true;
                if (a.totalH() >= 0 && a.degree() >= 1) g.basicAmine = true;
            }
        }
        if (a.z == 7) {
            // nitrile
            for (int bi : a.bonds) if (m.bonds[bi].order == 3.0) g.nitrile = true;
        }

        // Sulfur oxidation state.
        if (a.z == 16 && dbondTo(m, i, 8)) {
            int doubleO = 0;
            for (int bi : a.bonds)
                if (m.bonds[bi].order == 2.0 && m.atoms[m.bonds[bi].other(i)].z == 8) ++doubleO;
            if (doubleO >= 2) g.sulfone = true; else g.sulfoxide = true;
        }
    }

    // Phenol & catechol.
    for (int i = 0; i < n; ++i) {
        if (ringCarbonHasOH(m, i)) {
            g.phenol = true;
            for (int nb : m.atoms[i].nbr)
                if (m.atoms[nb].aromatic && ringCarbonHasOH(m, nb)) g.catechol = true;
        }
    }

    // Methylenedioxy: a CH2 carbon bonded to two ether oxygens that touch an arene.
    for (int i = 0; i < n; ++i) {
        const Atom& a = m.atoms[i];
        if (a.z != 6 || a.aromatic) continue;
        int oCount = 0; bool arene = false;
        for (int nb : a.nbr) {
            const Atom& o = m.atoms[nb];
            if (o.z == 8 && o.degree() == 2) { ++oCount; if (aromaticNeighbor(m, nb)) arene = true; }
        }
        if (oCount >= 2 && arene) g.methylenedioxy = true;
    }

    // Phenethylamine core: Ar-C-C-N (non-amide amine two carbons from an arene).
    for (int nIdx = 0; nIdx < n && !g.phenethylamine; ++nIdx) {
        const Atom& na = m.atoms[nIdx];
        if (na.z != 7 || na.aromatic) continue;
        bool amideN = false;
        for (int nb : na.nbr) if (isCarbonyl(m, nb)) amideN = true;
        if (amideN) continue;
        for (int c1 : na.nbr) {
            if (m.atoms[c1].z != 6 || m.atoms[c1].aromatic) continue;
            for (int c2 : m.atoms[c1].nbr) {
                if (c2 == nIdx || m.atoms[c2].z != 6) continue;
                if (m.atoms[c2].aromatic) { g.phenethylamine = true; break; }
                for (int c3 : m.atoms[c2].nbr)
                    if (c3 != c1 && m.atoms[c3].z == 6 && m.atoms[c3].aromatic) g.phenethylamine = true;
            }
        }
    }
    g.catecholamine = g.catechol && g.phenethylamine;
    return g;
}

std::vector<std::uint32_t> morganFingerprint(const Molecule& m, int radius) {
    const int n = static_cast<int>(m.atoms.size());
    std::vector<std::uint64_t> inv(n);
    for (int i = 0; i < n; ++i) {
        const Atom& a = m.atoms[i];
        std::uint64_t h = 1469598103934665603ULL;
        h = mix(h, static_cast<std::uint64_t>(a.z));
        h = mix(h, static_cast<std::uint64_t>(a.degree()));
        h = mix(h, static_cast<std::uint64_t>(a.totalH()));
        h = mix(h, static_cast<std::uint64_t>(a.charge + 8));
        h = mix(h, static_cast<std::uint64_t>(a.aromatic ? 1 : 0));
        h = mix(h, static_cast<std::uint64_t>(a.inRing ? 1 : 0));
        inv[i] = h;
    }

    std::vector<std::uint32_t> bits;
    auto emit = [&](std::uint64_t v) { bits.push_back(static_cast<std::uint32_t>(v ^ (v >> 32))); };
    for (int i = 0; i < n; ++i) emit(inv[i]);

    for (int r = 0; r < radius; ++r) {
        std::vector<std::uint64_t> next(n);
        for (int i = 0; i < n; ++i) {
            std::vector<std::uint64_t> env;
            for (int bi : m.atoms[i].bonds) {
                const Bond& b = m.bonds[bi];
                env.push_back(mix(static_cast<std::uint64_t>(b.order * 2), inv[b.other(i)]));
            }
            std::sort(env.begin(), env.end());
            std::uint64_t h = mix(inv[i], static_cast<std::uint64_t>(r + 1));
            for (std::uint64_t e : env) h = mix(h, e);
            next[i] = h;
        }
        inv.swap(next);
        for (int i = 0; i < n; ++i) emit(inv[i]);
    }

    std::sort(bits.begin(), bits.end());
    bits.erase(std::unique(bits.begin(), bits.end()), bits.end());
    return bits;
}

double tanimoto(const std::vector<std::uint32_t>& a, const std::vector<std::uint32_t>& b) {
    if (a.empty() && b.empty()) return 1.0;
    int inter = 0;
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] == b[j]) { ++inter; ++i; ++j; }
        else if (a[i] < b[j]) ++i;
        else ++j;
    }
    const int uni = static_cast<int>(a.size() + b.size()) - inter;
    return uni > 0 ? static_cast<double>(inter) / uni : 0.0;
}

}  // namespace biocad::chem
