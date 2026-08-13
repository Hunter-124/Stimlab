#include "chem/Descriptors.h"

#include <algorithm>
#include <map>
#include <numeric>
#include <vector>

namespace biocad::chem {

namespace {

bool hasBondOrder(const Molecule& m, int atom, double order) {
    for (int bi : m.atoms[atom].bonds)
        if (m.bonds[bi].order == order) return true;
    return false;
}

bool isDoubleBondedTo(const Molecule& m, int atom, int z) {
    for (int bi : m.atoms[atom].bonds) {
        const Bond& b = m.bonds[bi];
        if (b.order == 2.0 && m.atoms[b.other(atom)].z == z) return true;
    }
    return false;
}

// Is C-N bond an amide (the C carries a carbonyl =O)?
bool isAmideBond(const Molecule& m, const Bond& b) {
    int c = -1, n = -1;
    if (m.atoms[b.a].z == 6 && m.atoms[b.b].z == 7) { c = b.a; n = b.b; }
    else if (m.atoms[b.b].z == 6 && m.atoms[b.a].z == 7) { c = b.b; n = b.a; }
    if (c < 0) return false;
    (void)n;
    return isDoubleBondedTo(m, c, 8);
}

}  // namespace

int heavyAtomCount(const Molecule& m) { return static_cast<int>(m.atoms.size()); }

int formalCharge(const Molecule& m) {
    int q = 0;
    for (const auto& a : m.atoms) q += a.charge;
    return q;
}

double molecularWeight(const Molecule& m) {
    double w = 0.0;
    for (const auto& a : m.atoms) w += atomicMass(a.z) + a.totalH() * atomicMass(1);
    return w;
}

std::string molecularFormula(const Molecule& m) {
    std::map<std::string, int> counts;
    int h = 0;
    for (const auto& a : m.atoms) {
        counts[symbolByZ(a.z)]++;
        h += a.totalH();
    }
    if (h > 0) counts["H"] += h;

    auto emit = [](const std::string& sym, int n) {
        return n > 1 ? sym + std::to_string(n) : sym;
    };
    std::string out;
    // Hill order: C first, then H, then the rest alphabetically.
    if (counts.count("C")) { out += emit("C", counts["C"]); counts.erase("C"); }
    if (counts.count("H")) { out += emit("H", counts["H"]); counts.erase("H"); }
    for (const auto& [sym, n] : counts) out += emit(sym, n);
    return out;
}

int hbdCount(const Molecule& m) {
    int d = 0;
    for (const auto& a : m.atoms)
        if ((a.z == 7 || a.z == 8) && a.totalH() > 0) ++d;
    return d;
}

int hbaCount(const Molecule& m) {
    int n = 0;
    for (const auto& a : m.atoms)
        if (a.z == 7 || a.z == 8) ++n;
    return n;
}

int rotatableBondCount(const Molecule& m) {
    int n = 0;
    for (const auto& b : m.bonds) {
        if (b.order != 1.0 || b.aromatic || b.inRing) continue;
        if (m.atoms[b.a].degree() < 2 || m.atoms[b.b].degree() < 2) continue;
        if (isAmideBond(m, b)) continue;
        ++n;
    }
    return n;
}

int ringCount(const Molecule& m) {
    // Cyclomatic number = E - V + C (number of independent cycles = SSSR size).
    const int v = static_cast<int>(m.atoms.size());
    const int e = static_cast<int>(m.bonds.size());
    std::vector<int> comp(v, -1);
    int c = 0;
    for (int i = 0; i < v; ++i) {
        if (comp[i] != -1) continue;
        std::vector<int> stack{i};
        comp[i] = c;
        while (!stack.empty()) {
            const int u = stack.back(); stack.pop_back();
            for (int nb : m.atoms[u].nbr)
                if (comp[nb] == -1) { comp[nb] = c; stack.push_back(nb); }
        }
        ++c;
    }
    return e - v + c;
}

int aromaticAtomCount(const Molecule& m) {
    int n = 0;
    for (const auto& a : m.atoms) if (a.aromatic) ++n;
    return n;
}

double fractionCsp3(const Molecule& m) {
    int carbons = 0, sp3 = 0;
    for (size_t i = 0; i < m.atoms.size(); ++i) {
        const Atom& a = m.atoms[i];
        if (a.z != 6) continue;
        ++carbons;
        bool isSp3 = !a.aromatic;
        for (int bi : a.bonds)
            if (m.bonds[bi].order > 1.0) { isSp3 = false; break; }
        if (isSp3) ++sp3;
    }
    return carbons ? static_cast<double>(sp3) / carbons : 0.0;
}

// ----------------------------------------------------------------- Ertl TPSA
double tpsa(const Molecule& m) {
    double sum = 0.0;
    for (size_t i = 0; i < m.atoms.size(); ++i) {
        const Atom& a = m.atoms[i];
        const int deg = a.degree();
        const int h = a.totalH();
        const bool dbl = hasBondOrder(m, static_cast<int>(i), 2.0);
        const bool trp = hasBondOrder(m, static_cast<int>(i), 3.0);

        if (a.z == 7) {  // nitrogen
            if (a.aromatic) {
                if (a.charge == 0 && deg == 2 && h == 0) sum += 12.89;       // pyridine
                else if (a.charge == 0 && deg == 2 && h == 1) sum += 15.79;  // pyrrole NH
                else if (a.charge == 0 && deg == 3 && h == 0) sum += 4.41;   // bridgehead / N-subst
                else if (a.charge == 1 && deg == 3 && h == 0) sum += 4.10;
                else if (a.charge == 1 && deg == 3 && h == 1) sum += 14.14;
                else sum += 12.89;
            } else if (a.charge == 0) {
                if (trp) sum += 23.79;                                       // nitrile
                else if (deg == 1 && h == 2) sum += 26.02;                   // -NH2
                else if (deg == 1 && h == 1 && dbl) sum += 23.85;            // =NH
                else if (deg == 1 && h == 0 && dbl) sum += 23.85;
                else if (deg == 2 && h == 1) sum += 12.03;                   // -NH-
                else if (deg == 2 && h == 0 && dbl) sum += 12.36;            // =N-
                else if (deg == 2 && h == 0) sum += 12.03;
                else if (deg == 3 && h == 0 && dbl) sum += 11.68;            // e.g. amide-conjugated
                else if (deg == 3 && h == 0) sum += 3.24;                    // tertiary amine
                else if (deg == 3 && h == 1) sum += 3.01;
                else sum += 3.24;
            } else if (a.charge > 0) {
                if (deg == 4 && h == 0) sum += 0.00;
                else if (deg == 3 && h == 1) sum += 4.44;
                else if (deg == 2 && h == 0) sum += 13.97;
                else sum += 4.36;
            } else {
                sum += 12.03;
            }
        } else if (a.z == 8) {  // oxygen
            if (a.aromatic) sum += 13.14;                                    // furan-type
            else if (a.charge < 0) sum += 23.06;                             // O-
            else if (deg == 1 && h == 1) sum += 20.23;                       // -OH
            else if (deg == 1 && h == 0 && dbl) sum += 17.07;               // =O
            else if (deg == 2 && h == 0) sum += 9.23;                        // -O-
            else sum += 9.23;
        } else if (a.z == 16) {  // sulfur
            const bool dblToO = isDoubleBondedTo(m, static_cast<int>(i), 8);
            if (a.aromatic) sum += 28.24;
            else if (deg == 3 && dblToO) sum += 19.21;                       // sulfoxide
            else if (deg == 4 && dblToO) sum += 8.38;                        // sulfone
            else if (deg == 1 && dbl) sum += 32.09;                          // =S
            else if (deg == 2) sum += 25.30;                                 // -S-
            else sum += 25.30;
        }
    }
    return sum;
}

// ------------------------------------------------ Wildman-Crippen logP
// Atom-additive octanol-water logP after Wildman & Crippen (J. Chem. Inf.
// Comput. Sci. 1999, 39, 868): every atom is classified into an atom type by
// its element, hybridization, aromaticity, attached heteroatoms and hydrogen
// environment, then assigned the published per-type contribution. Hydrogens are
// folded into their heavy-atom carrier (H1/H2/H3 contributions). Validated to a
// tolerance against reference logP for the library in tests/test_chem.cpp.
namespace {

bool isHeteroZ(int z) {
    return z == 7 || z == 8 || z == 15 || z == 16 || z == 9 || z == 17 || z == 35 || z == 53;
}

// Contribution of the hydrogens carried by a heavy atom (WC H1/H2/H3 types).
double hydrogenContribution(int carrierZ, int nH) {
    if (nH <= 0) return 0.0;
    double per;
    switch (carrierZ) {
        case 6:  per = 0.1230; break;   // H1: hydrocarbon H
        case 7:  per = 0.2142; break;   // H3: H on nitrogen
        case 8:  per = -0.2677; break;  // H2: H on oxygen (alcohol/phenol/acid)
        case 16: per = 0.2142; break;   // thiol H ~ H3
        default: per = 0.1125; break;   // HS: default
    }
    return per * nH;
}

double carbonContribution(const Molecule& m, int i) {
    const Atom& a = m.atoms[i];
    const int H = a.totalH();
    int cN = 0, hetN = 0;
    for (int nb : a.nbr) {
        if (m.atoms[nb].z == 6) ++cN;
        else if (isHeteroZ(m.atoms[nb].z)) ++hetN;
    }
    bool dbl = false, trp = false, dblToHetero = false;
    for (int bi : a.bonds) {
        const Bond& b = m.bonds[bi];
        const int o = m.atoms[b.other(i)].z;
        if (b.order == 2.0) { dbl = true; if (o != 6) dblToHetero = true; }
        if (b.order == 3.0) trp = true;
    }
    if (a.aromatic) {
        if (H >= 1) return 0.1581;                    // C18: aromatic CH
        if (hetN > 0) return 0.2955;                  // C23/24: aromatic C-heteroatom
        return 0.1360;                                // C21: aromatic C-C(sub)
    }
    if (trp) return 0.0017;                           // C7: sp carbon
    if (dbl) return dblToHetero ? -0.2783 : 0.1551;   // C5 (=hetero) / C6 (=C)
    if (hetN > 0) return H >= 2 ? -0.2035 : -0.2051;  // C3 / C4: sp3 bonded to heteroatom
    return (cN <= 2 && H >= 1) ? 0.1441 : 0.0000;     // C1 / C2: sp3 all-carbon
}

double nitrogenContribution(const Molecule& m, int i) {
    const Atom& a = m.atoms[i];
    const int H = a.totalH();
    bool trp = false, dbl = false;
    for (int bi : a.bonds) {
        if (m.bonds[bi].order == 3.0) trp = true;
        if (m.bonds[bi].order == 2.0) dbl = true;
    }
    bool amide = false, arylAttached = false;
    for (int nb : a.nbr) {
        if (m.atoms[nb].aromatic) arylAttached = true;
        if (m.atoms[nb].z == 6)
            for (int bi : m.atoms[nb].bonds) {
                const Bond& b = m.bonds[bi];
                if (b.order == 2.0 && m.atoms[b.other(nb)].z == 8) amide = true;
            }
    }
    if (trp) return -0.3396;                           // N: nitrile
    if (a.aromatic) return H >= 1 ? -0.3239 : -0.3396; // N11/12: pyrrole NH / pyridine n
    if (amide) return -0.4458;                         // amide N
    if (dbl) return 0.08387;                           // N5: imine
    if (a.charge > 0) return -0.3187;                  // ammonium ~ tertiary
    if (arylAttached) return H >= 2 ? -1.0270 : (H == 1 ? -0.5188 : -0.3187);  // aniline N3/N4/N7
    return H >= 2 ? -1.0190 : (H == 1 ? -0.7096 : -0.3187);                    // amine N1/N2/N7
}

double oxygenContribution(const Molecule& m, int i) {
    const Atom& a = m.atoms[i];
    const int H = a.totalH();
    bool dbl = false;
    for (int bi : a.bonds)
        if (m.bonds[bi].order == 2.0) dbl = true;
    if (a.aromatic) return 0.1552;                     // O1: aromatic O (furan)
    if (a.charge < 0) return -0.3339;                  // O6: O-
    if (dbl) return -0.1188;                           // carbonyl =O (ketone/amide/acid/ester)
    if (H >= 1) return -0.2893;                         // O2: hydroxyl (alcohol/phenol/acid OH)
    return -0.0684;                                     // O3: ether -O- (incl. methylenedioxy)
}

}  // namespace

double crippenLogP(const Molecule& m) {
    double logp = 0.0;
    for (size_t i = 0; i < m.atoms.size(); ++i) {
        const Atom& a = m.atoms[i];
        const int idx = static_cast<int>(i);
        switch (a.z) {
            case 6:  logp += carbonContribution(m, idx); break;
            case 7:  logp += nitrogenContribution(m, idx); break;
            case 8:  logp += oxygenContribution(m, idx); break;
            case 9:  logp += 0.4202; break;   // F
            case 17: logp += 0.6895; break;   // Cl
            case 35: logp += 0.8456; break;   // Br
            case 53: logp += 0.8857; break;   // I
            case 16: logp += isDoubleBondedTo(m, idx, 8) ? -0.0024 : 0.6482; break;  // S=O / thioether
            case 15: logp += 0.8612; break;   // P
            default: break;
        }
        logp += hydrogenContribution(a.z, a.totalH());
    }
    return logp;
}

}  // namespace biocad::chem
