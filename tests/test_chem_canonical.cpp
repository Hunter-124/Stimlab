#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "chem/Aromaticity.h"
#include "chem/Canonical.h"
#include "chem/Rings.h"
#include "chem/Smiles.h"

using namespace biocad::chem;

namespace {

Molecule mustParse(const std::string& smi) {
    auto m = parseSmiles(smi);
    REQUIRE(m.has_value());
    return *m;
}

// Parse and perceive: aromaticity must come from the GRAPH, not from whether the
// author happened to type lowercase, or Kekule and aromatic input can never agree.
Molecule perceived(const std::string& smi) {
    Molecule m = mustParse(smi);
    const RingInfo rings = perceiveRings(m);
    annotateRings(m, rings);
    perceiveAromaticity(m, rings);
    return m;
}

// Relabel atoms by `p` and shuffle every adjacency list: the same graph presented
// in a different atom order. Canonicalisation must not notice.
Molecule permute(const Molecule& m, const std::vector<int>& p, std::mt19937& rng) {
    const int n = static_cast<int>(m.atoms.size());
    Molecule o;
    o.atoms.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        Atom a = m.atoms[static_cast<std::size_t>(i)];
        a.nbr.clear();
        a.bonds.clear();
        o.atoms[static_cast<std::size_t>(p[static_cast<std::size_t>(i)])] = a;
    }
    std::vector<int> order(m.bonds.size());
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), rng);
    for (int bi : order) {
        Bond b = m.bonds[static_cast<std::size_t>(bi)];
        b.a = p[static_cast<std::size_t>(b.a)];
        b.b = p[static_cast<std::size_t>(b.b)];
        const int idx = static_cast<int>(o.bonds.size());
        o.bonds.push_back(b);
        o.atoms[static_cast<std::size_t>(b.a)].nbr.push_back(b.b);
        o.atoms[static_cast<std::size_t>(b.a)].bonds.push_back(idx);
        o.atoms[static_cast<std::size_t>(b.b)].nbr.push_back(b.a);
        o.atoms[static_cast<std::size_t>(b.b)].bonds.push_back(idx);
    }
    return o;
}

struct Case {
    const char* name;
    std::vector<const char*> smiles;  // different inputs, same molecule
};

const std::vector<Case>& cases() {
    static const std::vector<Case> kCases = {
        {"caffeine", {"Cn1cnc2c1c(=O)n(C)c(=O)n2C", "Cn1c(=O)c2c(ncn2C)n(C)c1=O"}},
        {"cocaine", {"CN1C2CCC1C(C(=O)OC)C(OC(=O)c1ccccc1)C2",
                     "COC(=O)C1C2CCC(N2C)CC1OC(=O)c1ccccc1"}},
        {"mdma", {"CC(NC)Cc1ccc2OCOc2c1", "CNC(C)Cc1ccc2c(c1)OCO2"}},
        {"ibuprofen", {"CC(C)Cc1ccc(cc1)C(C)C(=O)O", "CC(C(=O)O)c1ccc(CC(C)C)cc1",
                       "OC(=O)C(C)c1ccc(cc1)CC(C)C"}},
        {"naphthalene", {"c1ccc2ccccc2c1", "c1ccc2c(c1)cccc2", "C1=CC2=CC=CC=C2C=C1"}},
        {"pyridine", {"n1ccccc1", "c1ccncc1", "c1ccccn1", "C1=CC=NC=C1"}},
        {"benzene", {"c1ccccc1", "C1=CC=CC=C1"}},
        {"glycine_zwitterion", {"[NH3+]CC(=O)[O-]", "[O-]C(=O)C[NH3+]", "C([NH3+])C(=O)[O-]"}},
        {"acetic_acid", {"CC(=O)O", "OC(C)=O", "OC(=O)C"}},
        {"toluene", {"Cc1ccccc1", "c1ccccc1C"}},
        {"phenol", {"Oc1ccccc1", "c1cc(O)ccc1"}},
        {"amphetamine", {"CC(N)Cc1ccccc1", "NC(C)Cc1ccccc1"}},
        {"aspirin", {"CC(=O)Oc1ccccc1C(=O)O", "O=C(O)c1ccccc1OC(C)=O"}},
        {"nicotine", {"CN1CCCC1c1cccnc1", "C1(c2cccnc2)CCCN1C"}},
        {"cyclohexane", {"C1CCCCC1"}},
        {"ethylbenzene_mix", {"CC.c1ccccc1", "c1ccccc1.CC"}},
    };
    return kCases;
}

}  // namespace

// THE load-bearing test. If different spellings of one molecule disagree, the
// writer is not canonical and nothing downstream (cache keys, dedup, identity)
// is trustworthy.
TEST_CASE("Canonical SMILES is independent of how the molecule was written",
          "[chem][canonical]") {
    for (const Case& c : cases()) {
        const std::string ref = canonicalSmiles(perceived(c.smiles.front()));
        REQUIRE_FALSE(ref.empty());
        for (const char* in : c.smiles) {
            INFO(c.name << " input " << in);
            CHECK(canonicalSmiles(perceived(in)) == ref);
        }
    }
}

// The same graph, atoms and adjacency lists shuffled. This isolates the writer
// from the parser: any residual dependence on input atom order shows up here.
TEST_CASE("Canonical SMILES is independent of atom ordering", "[chem][canonical]") {
    std::mt19937 rng(20260813u);
    for (const Case& c : cases()) {
        const Molecule m = perceived(c.smiles.front());
        const std::string ref = canonicalSmiles(m);
        for (int trial = 0; trial < 25; ++trial) {
            std::vector<int> p(m.atoms.size());
            std::iota(p.begin(), p.end(), 0);
            std::shuffle(p.begin(), p.end(), rng);
            INFO(c.name << " permutation trial " << trial);
            CHECK(canonicalSmiles(permute(m, p, rng)) == ref);
        }
    }
}

TEST_CASE("Canonical SMILES round-trips through the parser", "[chem][canonical]") {
    for (const Case& c : cases()) {
        for (const char* in : c.smiles) {
            const std::string once = canonicalSmiles(perceived(in));
            INFO(c.name << " input " << in << " canonical " << once);
            auto again = parseSmiles(once);
            REQUIRE(again.has_value());
            RingInfo rings = perceiveRings(*again);
            annotateRings(*again, rings);
            perceiveAromaticity(*again, rings);
            CHECK(canonicalSmiles(*again) == once);
        }
    }
}

TEST_CASE("Distinct molecules get distinct canonical strings and hashes",
          "[chem][canonical]") {
    std::map<std::string, std::string> byString;
    std::map<std::uint64_t, std::string> byHash;
    for (const Case& c : cases()) {
        const Molecule m = perceived(c.smiles.front());
        const std::string cs = canonicalSmiles(m);
        const std::uint64_t h = graphHash(m);

        auto its = byString.find(cs);
        INFO(c.name << " canonical " << cs);
        CHECK(its == byString.end());
        auto ith = byHash.find(h);
        CHECK(ith == byHash.end());

        byString[cs] = c.name;
        byHash[h] = c.name;
    }
    CHECK(byString.size() == cases().size());
    CHECK(byHash.size() == cases().size());
}

// A disconnected input must not depend on which component the author typed first;
// components are ordered by their own canonical string.
TEST_CASE("Disconnected components are ordered canonically", "[chem][canonical]") {
    const std::string a = canonicalSmiles(perceived("CC.c1ccccc1"));
    const std::string b = canonicalSmiles(perceived("c1ccccc1.CC"));
    CHECK(a == b);
    CHECK(a.find('.') != std::string::npos);
}

// canonicalRanks is a total order, not merely a refined partition.
TEST_CASE("Canonical ranks are a permutation of 0..n-1", "[chem][canonical]") {
    for (const Case& c : cases()) {
        const Molecule m = perceived(c.smiles.front());
        std::vector<int> r = canonicalRanks(m);
        REQUIRE(r.size() == m.atoms.size());
        std::sort(r.begin(), r.end());
        for (std::size_t i = 0; i < r.size(); ++i) {
            INFO(c.name << " rank slot " << i);
            REQUIRE(r[i] == static_cast<int>(i));
        }
    }
}

// The hash folds in atom and bond counts, so a mismatch in either cannot collide
// even before the string is considered.
TEST_CASE("graphHash separates same-string-length neighbours", "[chem][canonical]") {
    const std::uint64_t benzene = graphHash(perceived("c1ccccc1"));
    const std::uint64_t pyridine = graphHash(perceived("c1ccncc1"));
    const std::uint64_t cyclohexane = graphHash(perceived("C1CCCCC1"));
    CHECK(benzene != pyridine);
    CHECK(benzene != cyclohexane);
    CHECK(pyridine != cyclohexane);
    CHECK(graphHash(perceived("c1ccccc1")) == benzene);  // stable across calls
}
