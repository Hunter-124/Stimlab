#include "chem/Descriptors.h"

#include <algorithm>
#include <map>
#include <numeric>
#include <vector>

namespace stimlab::chem {

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

// ------------------------------------------------------ Crippen-style logP
// Pragmatic atom-contribution estimate (not full Wildman-Crippen typing);
// validated to a tolerance in tests. Swappable for RDKit's MolLogP later.
double crippenLogP(const Molecule& m) {
    double logp = 0.0;
    for (size_t i = 0; i < m.atoms.size(); ++i) {
        const Atom& a = m.atoms[i];
        switch (a.z) {
            case 6: {  // carbon
                bool boundToHetero = false;
                for (int nb : a.nbr)
                    if (m.atoms[nb].z != 1 && m.atoms[nb].z != 6) boundToHetero = true;
                double cval = a.aromatic ? 0.290 : 0.205;
                if (boundToHetero) cval -= 0.151;
                logp += cval;
                logp += a.totalH() * 0.110;  // nonpolar H
                break;
            }
            case 7:  // nitrogen
                logp += a.aromatic ? -0.40 : -0.90;
                logp += a.totalH() * (-0.20);
                break;
            case 8:  // oxygen
                logp += isDoubleBondedTo(m, static_cast<int>(i), 6) ? -0.20 : -0.38;
                logp += a.totalH() * (-0.25);
                break;
            case 9:  logp += 0.40; break;   // F
            case 17: logp += 0.64; break;   // Cl
            case 35: logp += 0.84; break;   // Br
            case 53: logp += 1.00; break;   // I
            case 16: logp += isDoubleBondedTo(m, static_cast<int>(i), 8) ? -0.40 : 0.62; break;  // S
            case 15: logp += 0.30; break;   // P
            default: break;
        }
    }
    return logp;
}

}  // namespace stimlab::chem
