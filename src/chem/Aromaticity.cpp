#include "chem/Aromaticity.h"

#include <algorithm>
#include <optional>
#include <vector>

namespace biocad::chem {

namespace {

// Ring bonds are looked up through a per-ring mask so "is this double bond
// inside the ring I am counting?" is an O(1) question.
using Mask = std::vector<char>;

// A ring bond that carries pi density: a real double bond (Kekule input) or a
// bond the parser marked aromatic (order 1.5, lowercase SMILES input).
bool piRingBond(const Molecule& m, int atom, const Mask& ringBond) {
    for (int bi : m.atoms[static_cast<std::size_t>(atom)].bonds) {
        if (!ringBond[static_cast<std::size_t>(bi)]) continue;
        if (m.bonds[static_cast<std::size_t>(bi)].order >= 1.5) return true;
    }
    return false;
}

// A localised C=X / N=X double bond inside the ring. Heteroatoms must use this
// strict form: for lowercase input every ring bond has order 1.5, and treating
// that as "double" would make pyrrole's [nH] a pyridine-type 1-electron donor.
bool strictDoubleInRing(const Molecule& m, int atom, const Mask& ringBond) {
    for (int bi : m.atoms[static_cast<std::size_t>(atom)].bonds) {
        if (!ringBond[static_cast<std::size_t>(bi)]) continue;
        if (m.bonds[static_cast<std::size_t>(bi)].order == 2.0) return true;
    }
    return false;
}

// Partner of a localised double bond that leaves the ring, or -1.
int exocyclicDoublePartner(const Molecule& m, int atom, const Mask& ringBond) {
    for (int bi : m.atoms[static_cast<std::size_t>(atom)].bonds) {
        if (ringBond[static_cast<std::size_t>(bi)]) continue;
        const Bond& b = m.bonds[static_cast<std::size_t>(bi)];
        if (b.order == 2.0) return b.other(atom);
    }
    return -1;
}

// Pi-electron contribution of one ring atom, or nullopt when the atom cannot be
// sp2 and therefore disqualifies the whole ring. The table this implements, with
// the reason for every row, is in Aromaticity.h - keep the two in step.
std::optional<int> contribution(const Molecule& m, int i, const Mask& ringBond,
                                const std::vector<char>& atomAromatic) {
    const Atom& a = m.atoms[static_cast<std::size_t>(i)];
    switch (a.z) {
        case 6: {  // carbon
            // The exocyclic double bond is tested FIRST, and deliberately so:
            // lowercase input gives every ring bond order 1.5, so asking
            // "does a ring pi bond exist?" first would score a lowercase
            // carbonyl carbon (caffeine's c(=O), 2-pyridone) as 1 electron and
            // make perception depend on the input spelling. A carbon can only
            // carry one double bond anyway, so this ordering loses nothing.
            const int exo = exocyclicDoublePartner(m, i, ringBond);
            if (exo >= 0) {
                const int pz = m.atoms[static_cast<std::size_t>(exo)].z;
                if (pz == 6) {
                    // Fused-bond conjugation: only counts once the partner ring
                    // is known aromatic. This is the fixed-point rule.
                    return atomAromatic[static_cast<std::size_t>(exo)] ? 1 : 0;
                }
                return 0;  // carbonyl / thiocarbonyl / imine carbon: sp2, 0 e-
            }
            if (piRingBond(m, i, ringBond)) return 1;
            if (a.charge == -1) return 2;  // cyclopentadienyl-type carbanion
            if (a.charge == 1) return 0;   // tropylium-type carbocation
            return std::nullopt;           // sp3: kills the ring
        }
        case 7:
        case 15: {  // nitrogen, phosphorus
            if (strictDoubleInRing(m, i, ringBond)) return 1;  // pyridine-type
            if (a.charge == -1) return 2;                      // azolate anion
            // Cationic nitrogen: a pyridinium, thiazolium or isoquinolinium N has
            // spent its lone pair on the fourth sigma bond (or on a proton), so it
            // donates one p electron like pyridine, NOT two like pyrrole. Tested
            // before the pyrrole-type row because an N-methylpyridinium also has
            // three sigma connections, and lowercase input gives every ring bond
            // order 1.5, so the formal charge is the only thing that tells the two
            // apart. Without this row N-methylated cofactors (NMN/NAD+, thiamine,
            // berberine) come out non-aromatic and every aromatic rule pack misses
            // them.
            if (a.charge == 1 && piRingBond(m, i, ringBond)) return 1;
            if (a.totalH() >= 1 || a.degree() == 3) return 2;   // pyrrole-type
            if (piRingBond(m, i, ringBond)) return 1;  // lowercase aromatic n
            return std::nullopt;
        }
        case 8:
        case 16:
        case 34: {  // oxygen, sulfur, selenium
            if (strictDoubleInRing(m, i, ringBond))
                return a.charge == 1 ? std::optional<int>(1)  // pyrylium-type
                                     : std::nullopt;
            if (a.degree() == 2) return 2;  // furan / thiophene lone pair
            return std::nullopt;
        }
        case 5:  // boron: empty p orbital
            return a.degree() <= 3 ? std::optional<int>(0) : std::nullopt;
        default:
            return std::nullopt;
    }
}

}  // namespace

void perceiveAromaticity(Molecule& m, const RingInfo& info) {
    for (Atom& a : m.atoms) a.aromatic = false;
    for (Bond& b : m.bonds) b.aromatic = false;
    if (info.count() == 0) return;

    const std::size_t nRings = info.count();
    std::vector<Mask> ringBond(nRings, Mask(m.bonds.size(), 0));
    for (std::size_t r = 0; r < nRings; ++r)
        for (int bi : info.bondRings[r]) ringBond[r][static_cast<std::size_t>(bi)] = 1;

    std::vector<char> ringArom(nRings, 0);
    std::vector<char> atomArom(m.atoms.size(), 0);

    // Monotone fixed point: a ring may only ever be promoted to aromatic, so the
    // loop cannot oscillate and is bounded by the ring count.
    for (std::size_t pass = 0; pass <= nRings; ++pass) {
        bool changed = false;
        for (std::size_t r = 0; r < nRings; ++r) {
            if (ringArom[r]) continue;
            int total = 0;
            bool ok = true;
            for (int i : info.atomRings[r]) {
                const std::optional<int> e = contribution(m, i, ringBond[r], atomArom);
                if (!e) {
                    ok = false;
                    break;
                }
                total += *e;
            }
            // Huckel 4n + 2, n >= 0.
            if (!ok || total < 2 || (total - 2) % 4 != 0) continue;
            ringArom[r] = 1;
            for (int i : info.atomRings[r]) atomArom[static_cast<std::size_t>(i)] = 1;
            changed = true;
        }
        if (!changed) break;
    }

    for (std::size_t r = 0; r < nRings; ++r) {
        if (!ringArom[r]) continue;
        for (int i : info.atomRings[r]) m.atoms[static_cast<std::size_t>(i)].aromatic = true;
        for (int bi : info.bondRings[r]) m.bonds[static_cast<std::size_t>(bi)].aromatic = true;
    }
}

RingInfo perceiveRingsAndAromaticity(Molecule& m) {
    RingInfo info = perceiveRings(m);
    annotateRings(m, info);
    perceiveAromaticity(m, info);
    return info;
}

void normalizeAromaticBondOrders(Molecule& m) {
    for (Bond& b : m.bonds)
        if (b.aromatic) b.order = 1.5;
}

}  // namespace biocad::chem
