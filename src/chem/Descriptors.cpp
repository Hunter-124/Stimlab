#include "chem/Descriptors.h"

#include "chem/Crippen.h"
#include "chem/Smarts.h"

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
// The published contributions are keyed on aromaticity, so a caller that hands
// over a raw parse of a Kekule SMILES gets a different area for the same
// molecule (caffeine: 56.22 from CN1C=NC2=C1... vs 60.26 perceived). Every one
// of the ~15 call sites did exactly that, so perception happens here rather than
// being a rule callers must remember - the same decision chem::crippen() makes.
static double tpsaPerceived(const Molecule& m, bool includeSulfurAndPhosphorus) {
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
        } else if (a.z == 16 && includeSulfurAndPhosphorus) {  // sulfur
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

// Ertl's paper tabulates sulfur and phosphorus contributions, but every published
// TPSA THRESHOLD - Veber's <=140 A^2 for oral bioavailability, the <=90 A^2 rule of
// thumb for CNS penetration - was derived on the N,O-only sum, so that is the
// default here. Including S and P is a different number, not a better one: on
// famotidine the two conventions differ by 62 A^2, which straddles the Veber cutoff.
// Verified against RDKit 2026.03.5 over all 69 shipped library compounds.
double tpsa(const Molecule& m) {
    return tpsaPerceived(prepareMolecule(m).mol, /*includeSulfurAndPhosphorus=*/false);
}

double tpsaIncludingSulfurAndPhosphorus(const Molecule& m) {
    return tpsaPerceived(prepareMolecule(m).mol, /*includeSulfurAndPhosphorus=*/true);
}

// ------------------------------------------------ Wildman-Crippen logP
// WHY THIS IS A ONE-LINE DELEGATION NOW: the implementation that lived here was
// an ad-hoc re-derivation of the method - it invented its own decision tree over
// degree/heteroatom counts instead of the paper's class table, and it folded all
// hydrogens on a heavy atom into one contribution keyed only on the carrier's
// element. Measured against experiment on a 14-compound reference set it
// diverged from the published method by up to 1.35 log units on a single
// compound (caffeine, ibuprofen and dopamine were the worst offenders). The
// faithful implementation is chem::Crippen: 106 ordered SMARTS entries covering
// the method's 67 heavy-atom classes (several classes need more than one pattern)
// are data in assets/packs/descriptors/crippen.json, and the four hydrogen
// classes are derived in Crippen.cpp because this graph has no explicit H atoms.
//
// The symbol is kept because ~15 call sites use it. It returns 0.0 when the
// descriptor pack cannot be loaded; a caller that must distinguish "zero logP"
// from "no parameters" calls chem::crippen() directly and reads CrippenResult::ok.
double crippenLogP(const Molecule& m) { return crippen(m).logP; }

}  // namespace biocad::chem
