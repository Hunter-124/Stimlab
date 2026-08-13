// Ring perception (SSSR) and graph-based aromaticity perception.
//
// The load-bearing property under test is INPUT-FORM INDEPENDENCE: a lowercase
// aromatic SMILES and its Kekule spelling describe the same graph and must
// therefore perceive identically. Before Rings/Aromaticity existed the engine
// trusted the input, so C1=CC=CC=C1 was not benzene as far as any downstream
// rule was concerned. Every expectation here is the value the documented model
// in Aromaticity.h actually yields (verified by execution), including the cases
// the model deliberately gets "wrong" (2-pyridone aromatic, azulene not).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <string>
#include <vector>

#include "chem/Aromaticity.h"
#include "chem/Rings.h"
#include "chem/Smiles.h"

using namespace biocad::chem;

namespace {

struct Perceived {
    Molecule m;
    RingInfo info;
};

Perceived perceive(const std::string& smi) {
    auto parsed = parseSmiles(smi);
    REQUIRE(parsed.has_value());
    Perceived p;
    p.m = *parsed;
    p.info = perceiveRingsAndAromaticity(p.m);
    return p;
}

int aromaticAtoms(const Molecule& m) {
    return static_cast<int>(std::count_if(m.atoms.begin(), m.atoms.end(),
                                          [](const Atom& a) { return a.aromatic; }));
}
int aromaticBonds(const Molecule& m) {
    return static_cast<int>(std::count_if(m.bonds.begin(), m.bonds.end(),
                                          [](const Bond& b) { return b.aromatic; }));
}
std::vector<int> ringSizes(const RingInfo& info) {
    std::vector<int> sizes;
    for (const auto& r : info.atomRings) sizes.push_back(static_cast<int>(r.size()));
    std::sort(sizes.begin(), sizes.end());
    return sizes;
}

// Euler: |bonds| - |atoms| + |components|. The SSSR must contain exactly this
// many rings, which is the invariant that makes the greedy GF(2) search correct.
int cycleSpaceDim(const Molecule& m) {
    std::vector<char> seen(m.atoms.size(), 0);
    int comps = 0;
    std::vector<int> stack;
    for (std::size_t s = 0; s < m.atoms.size(); ++s) {
        if (seen[s]) continue;
        ++comps;
        stack.push_back(static_cast<int>(s));
        seen[s] = 1;
        while (!stack.empty()) {
            const int u = stack.back();
            stack.pop_back();
            for (int v : m.atoms[static_cast<std::size_t>(u)].nbr) {
                if (seen[static_cast<std::size_t>(v)]) continue;
                seen[static_cast<std::size_t>(v)] = 1;
                stack.push_back(v);
            }
        }
    }
    const int dim = static_cast<int>(m.bonds.size()) - static_cast<int>(m.atoms.size()) + comps;
    return dim < 0 ? 0 : dim;
}

}  // namespace

TEST_CASE("SSSR ring count always equals the Euler cycle-space dimension", "[chem][rings]") {
    const char* smiles[] = {
        "CC.CC", "CCO", "C1CCCCC1", "c1ccccc1", "c1ccc2ccccc2c1", "C1=CC=C2C=CC=CC2=C1",
        "C1C2CC3CC1CC(C2)C3",        // adamantane
        "C12C3C4C1C5C4C3C25",        // cubane
        "C1CCC2(C1)CCCCC2",          // spiro[4.5]decane
        "c1ccc2cc3ccccc3cc2c1",      // anthracene
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "c1ccccc1.C1CCCCC1",
    };
    for (const char* s : smiles) {
        Perceived p = perceive(s);
        INFO("smiles " << s);
        REQUIRE(static_cast<int>(p.info.count()) == cycleSpaceDim(p.m));
    }
}

TEST_CASE("bondRings walks the same cycle as atomRings", "[chem][rings]") {
    // Downstream code (SMARTS ring bonds, cartoon/geometry consumers) relies on
    // bondRings[i][k] joining atomRings[i][k] to atomRings[i][k+1].
    for (const char* s : {"c1ccc2[nH]ccc2c1", "C12C3C4C1C5C4C3C25", "C1CCC2(C1)CCCCC2"}) {
        Perceived p = perceive(s);
        INFO("smiles " << s);
        for (std::size_t r = 0; r < p.info.count(); ++r) {
            const auto& ats = p.info.atomRings[r];
            const auto& bds = p.info.bondRings[r];
            REQUIRE(ats.size() == bds.size());
            for (std::size_t k = 0; k < ats.size(); ++k) {
                const Bond& b = p.m.bonds[static_cast<std::size_t>(bds[k])];
                const int u = ats[k];
                const int v = ats[(k + 1) % ats.size()];
                REQUIRE(((b.a == u && b.b == v) || (b.a == v && b.b == u)));
            }
        }
    }
}

TEST_CASE("Benzene is perceived identically from aromatic and Kekule input", "[chem][aromatic]") {
    for (const char* s : {"c1ccccc1", "C1=CC=CC=C1"}) {
        Perceived p = perceive(s);
        INFO("smiles " << s);
        REQUIRE(p.info.count() == 1);
        REQUIRE(p.info.atomRings[0].size() == 6);
        REQUIRE(aromaticAtoms(p.m) == 6);
        REQUIRE(aromaticBonds(p.m) == 6);
    }
}

TEST_CASE("Cyclohexane has one ring and no aromaticity", "[chem][aromatic]") {
    Perceived p = perceive("C1CCCCC1");
    REQUIRE(p.info.count() == 1);
    REQUIRE(p.info.atomRings[0].size() == 6);
    REQUIRE(aromaticAtoms(p.m) == 0);
}

TEST_CASE("Fused systems resolve to a fixed point", "[chem][aromatic]") {
    SECTION("naphthalene, both spellings: 2 rings, 10 aromatic atoms") {
        // The Kekule form is the reason the fixed-point pass exists: the second
        // ring's fusion carbons point their double bonds into the first ring, so
        // in isolation it counts only 4 electrons.
        for (const char* s : {"c1ccc2ccccc2c1", "C1=CC=C2C=CC=CC2=C1"}) {
            Perceived p = perceive(s);
            INFO("smiles " << s);
            REQUIRE(p.info.count() == 2);
            REQUIRE(ringSizes(p.info) == std::vector<int>{6, 6});
            REQUIRE(aromaticAtoms(p.m) == 10);
            REQUIRE(aromaticBonds(p.m) == 11);
        }
    }
    SECTION("indole, both spellings: both rings aromatic, 9 atoms") {
        for (const char* s : {"c1ccc2[nH]ccc2c1", "C1=CC=C2NC=CC2=C1"}) {
            Perceived p = perceive(s);
            INFO("smiles " << s);
            REQUIRE(p.info.count() == 2);
            REQUIRE(ringSizes(p.info) == std::vector<int>{5, 6});
            REQUIRE(aromaticAtoms(p.m) == 9);
        }
    }
    SECTION("anthracene: 3 rings, 14 aromatic atoms") {
        Perceived p = perceive("c1ccc2cc3ccccc3cc2c1");
        REQUIRE(p.info.count() == 3);
        REQUIRE(aromaticAtoms(p.m) == 14);
    }
}

TEST_CASE("Heteroaromatics follow the documented electron table", "[chem][aromatic]") {
    struct Case {
        const char* smiles;
        int aromatic;
    };
    // Each pair is (aromatic spelling, Kekule spelling) where both exist.
    const Case cases[] = {
        {"c1ccncc1", 6},    {"C1=CC=NC=C1", 6},  // pyridine N: lone pair in sigma plane, 1 e-
        {"c1cc[nH]c1", 5},  {"C1=CC=CN1", 5},    // pyrrole N-H: lone pair donates 2 e-
        {"Cn1cccc1", 5},                          // N-methylpyrrole: 3 sigma bonds, 2 e-
        {"c1ccoc1", 5},     {"C1=CC=CO1", 5},    // furan O: one lone pair donates 2 e-
        {"c1ccsc1", 5},                           // thiophene S: 2 e-
        {"c1cnc[nH]1", 5},  {"C1=CN=CN1", 5},    // imidazole: 3*1 + 1 + 2 = 6 e-
        {"c1ncc2[nH]cnc2n1", 9},                  // purine: both rings aromatic
        {"[CH-]1C=CC=C1", 5},                     // Cp anion: 2 + 4*1 = 6 e-
    };
    for (const Case& c : cases) {
        Perceived p = perceive(c.smiles);
        INFO("smiles " << c.smiles);
        REQUIRE(aromaticAtoms(p.m) == c.aromatic);
    }
}

TEST_CASE("Non-aromatic rings are rejected for the documented reason", "[chem][aromatic]") {
    SECTION("cyclopentadiene: one sp3 CH2 disqualifies the ring") {
        Perceived p = perceive("C1=CCC=C1");
        REQUIRE(p.info.count() == 1);
        REQUIRE(aromaticAtoms(p.m) == 0);
    }
    SECTION("cyclobutadiene: 4 electrons is not 4n+2") {
        Perceived p = perceive("C1=CC=C1");
        REQUIRE(aromaticAtoms(p.m) == 0);
    }
    SECTION("para-benzoquinone: two exocyclic C=O give 0 e- each, total 4") {
        Perceived p = perceive("O=C1C=CC(=O)C=C1");
        REQUIRE(p.info.count() == 1);
        REQUIRE(aromaticAtoms(p.m) == 0);
    }
    SECTION("azulene: the 10-electron system is on the perimeter, not on an SSSR ring") {
        // Documented model limitation, asserted so it cannot change silently.
        Perceived p = perceive("c1ccc2cccc2cc1");
        REQUIRE(ringSizes(p.info) == std::vector<int>{5, 7});
        REQUIRE(aromaticAtoms(p.m) == 0);
    }
}

TEST_CASE("2-pyridone is reported aromatic by this model", "[chem][aromatic]") {
    // Exocyclic C=O contributes 0, the amide N contributes 2, the four CH carbons
    // contribute 1 each -> 6 e-. RDKit's default model agrees; a chemist calling
    // this an amide is not wrong. Asserted so the choice stays deliberate.
    Perceived p = perceive("O=C1C=CC=CN1");
    REQUIRE(p.info.count() == 1);
    REQUIRE(aromaticAtoms(p.m) == 6);
}

TEST_CASE("Caffeine perceives identically from four different spellings", "[chem][aromatic]") {
    // 9 ring atoms (5-ring 5 + 6-ring 6 - 2 shared fusion atoms) and 10 distinct
    // ring bonds (5 + 6 - 1 shared). The two exocyclic C=O carbons are ring
    // members that contribute 0 electrons, which is exactly why the
    // pyrimidinedione ring reaches 6 and not 8.
    for (const char* s : {"CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
                          "Cn1c(=O)c2c(ncn2C)n(C)c1=O", "n1cn(C)c2C(=O)N(C)C(=O)N(C)c12"}) {
        Perceived p = perceive(s);
        INFO("smiles " << s);
        REQUIRE(p.info.count() == 2);
        REQUIRE(ringSizes(p.info) == std::vector<int>{5, 6});
        REQUIRE(aromaticAtoms(p.m) == 9);
        REQUIRE(aromaticBonds(p.m) == 10);
    }
}

TEST_CASE("Biphenyl keeps its inter-ring single bond non-aromatic", "[chem][aromatic]") {
    Perceived p = perceive("c1ccccc1-c1ccccc1");
    REQUIRE(p.info.count() == 2);
    REQUIRE(aromaticAtoms(p.m) == 12);
    REQUIRE(aromaticBonds(p.m) == 12);  // 6 + 6; the biaryl bond is in no ring
}

TEST_CASE("Cage and spiro systems get the Euler ring count and terminate", "[chem][rings]") {
    SECTION("cubane: 12 bonds, 8 atoms -> 5 rings, not its 6 faces") {
        Perceived p = perceive("C12C3C4C1C5C4C3C25");
        REQUIRE(p.m.bonds.size() == 12);
        REQUIRE(p.info.count() == 5);
        REQUIRE(ringSizes(p.info) == std::vector<int>{4, 4, 4, 4, 4});
    }
    SECTION("adamantane: 3 six-membered rings") {
        Perceived p = perceive("C1C2CC3CC1CC(C2)C3");
        REQUIRE(ringSizes(p.info) == std::vector<int>{6, 6, 6});
    }
    SECTION("spiro[4.5]decane: a 5-ring and a 6-ring sharing one atom") {
        Perceived p = perceive("C1CCC2(C1)CCCCC2");
        REQUIRE(ringSizes(p.info) == std::vector<int>{5, 6});
    }
}

TEST_CASE("Disconnected input has the right cycle space", "[chem][rings]") {
    SECTION("two acyclic fragments: no rings, no crash") {
        Perceived p = perceive("CC.CC");
        REQUIRE(p.info.count() == 0);
        REQUIRE(aromaticAtoms(p.m) == 0);
    }
    SECTION("one aromatic and one saturated fragment") {
        Perceived p = perceive("c1ccccc1.C1CCCCC1");
        REQUIRE(p.info.count() == 2);
        REQUIRE(aromaticAtoms(p.m) == 6);
    }
}

TEST_CASE("Ring-size queries answer the SMARTS primitives", "[chem][rings]") {
    Perceived p = perceive("c1ccc2[nH]ccc2c1");  // indole
    int fusion = -1, sixOnly = -1;
    for (int i = 0; i < static_cast<int>(p.m.atoms.size()); ++i) {
        if (ringCountOf(p.info, i) == 2) fusion = i;
        else if (ringCountOf(p.info, i) == 1 && inRingOfSize(p.info, i, 6)) sixOnly = i;
    }
    REQUIRE(fusion >= 0);
    REQUIRE(sixOnly >= 0);
    REQUIRE(ringSizeOf(p.info, fusion) == 5);  // smallest ring containing the atom
    REQUIRE(inRingOfSize(p.info, fusion, 5));
    REQUIRE(inRingOfSize(p.info, fusion, 6));
    REQUIRE(ringSizeOf(p.info, sixOnly) == 6);
    REQUIRE_FALSE(inRingOfSize(p.info, sixOnly, 5));

    // annotateRings replaced the parser's crude bridge flag on every ring member.
    for (const auto& ring : p.info.atomRings)
        for (int a : ring) REQUIRE(p.m.atoms[static_cast<std::size_t>(a)].inRing);

    // A side-chain atom is in no ring at all.
    Perceived tol = perceive("Cc1ccccc1");
    REQUIRE(ringSizeOf(tol.info, 0) == 0);
    REQUIRE(ringCountOf(tol.info, 0) == 0);
    REQUIRE_FALSE(tol.m.atoms[0].inRing);
}

TEST_CASE("Perception is deterministic for a non-unique SSSR", "[chem][rings]") {
    // Cubane's SSSR must discard one of six equivalent faces; which one is
    // arbitrary, so the choice is pinned to be reproducible instead.
    Perceived a = perceive("C12C3C4C1C5C4C3C25");
    Perceived b = perceive("C12C3C4C1C5C4C3C25");
    REQUIRE(a.info.atomRings == b.info.atomRings);
    REQUIRE(a.info.bondRings == b.info.bondRings);
}

TEST_CASE("normalizeAromaticBondOrders converges Kekule and aromatic input",
          "[chem][aromatic]") {
    Perceived kek = perceive("C1=CC=CC=C1");
    Perceived low = perceive("c1ccccc1");
    normalizeAromaticBondOrders(kek.m);
    normalizeAromaticBondOrders(low.m);
    for (const Bond& b : kek.m.bonds) REQUIRE(b.order == 1.5);
    for (const Bond& b : low.m.bonds) REQUIRE(b.order == 1.5);
}
