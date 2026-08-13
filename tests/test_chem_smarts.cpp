// Tests for chem/Smarts: the SMARTS parser and the VF2 monomorphism matcher.
// The cases here are the patterns real rule sets use - a structural-alert or
// biotransformation pack is only as trustworthy as these counts.
#include <catch2/catch_test_macros.hpp>

#include <string>

#include "chem/Smarts.h"
#include "chem/Smiles.h"

using namespace biocad::chem;

namespace {

// Ring and aromaticity perception are prerequisites of the ring/aromatic
// primitives, so every fixture goes through prepareMolecule.
PreparedMolecule prepared(const std::string& smiles) {
    auto m = parseSmiles(smiles);
    REQUIRE(m.has_value());
    return prepareMolecule(*m);
}

std::size_t hits(const std::string& smarts, const PreparedMolecule& pm) {
    std::string err;
    auto p = parseSmarts(smarts, &err);
    REQUIRE(err.empty());
    REQUIRE(p.has_value());
    return matchAll(*p, pm.mol, pm.rings, 10000).size();
}

}  // namespace

TEST_CASE("SMARTS hydroxyl oxygen is not a carbonyl oxygen", "[chem][smarts]") {
    CHECK(hits("[OX2H]", prepared("CCO")) == 1);            // ethanol
    CHECK(hits("[OX2H]", prepared("c1ccc(O)cc1")) == 1);    // phenol
    CHECK(hits("[OX2H]", prepared("CC(=O)C")) == 0);        // acetone
}

TEST_CASE("SMARTS carboxylic acid distinguishes an acid from its ester", "[chem][smarts]") {
    const auto acid = prepared("CC(C)Cc1ccc(cc1)C(C)C(=O)O");         // ibuprofen
    const auto ester = prepared("CC(C)Cc1ccc(cc1)C(C)C(=O)OC");       // methyl ester
    CHECK(hits("[CX3](=O)[OX2H1]", acid) == 1);
    CHECK(hits("[CX3](=O)[OX2H1]", ester) == 0);
}

TEST_CASE("SMARTS aromatic ring matching uses the matched-atom-set convention",
          "[chem][smarts]") {
    // Matches are distinct by SET of matched atoms, so the six automorphisms of
    // benzene collapse to one match and naphthalene reports its two six-cycles.
    CHECK(hits("c1ccccc1", prepared("c1ccccc1")) == 1);
    CHECK(hits("c1ccccc1", prepared("c1ccc2ccccc2c1")) == 2);
    CHECK(hits("c1ccccc1", prepared("Cc1ccccc1")) == 1);
    CHECK(hits("c1ccccc1", prepared("C1CCCCC1")) == 0);
    // Aromaticity is perceived from the graph, so Kekule uppercase input works.
    CHECK(hits("c1ccccc1", prepared("C1=CC=CC=C1")) == 1);
}

TEST_CASE("SMARTS primary/secondary non-amide amine alert", "[chem][smarts]") {
    // The single most important case: ';' , ',' and recursive $() together.
    const std::string alert = "[NX3;H2,H1;!$(NC=O)]";
    CHECK(hits(alert, prepared("CC(N)Cc1ccccc1")) == 1);      // amphetamine
    CHECK(hits(alert, prepared("CNC(C)Cc1ccccc1")) == 1);     // methamphetamine
    CHECK(hits(alert, prepared("CC(=O)Nc1ccc(O)cc1")) == 0);  // acetaminophen anilide
}

TEST_CASE("SMARTS any-bond and atomic-number primitives", "[chem][smarts]") {
    CHECK(hits("[#6]~[#7]", prepared("CC(N)Cc1ccccc1")) == 1);
    CHECK(hits("[#6]~[#7]", prepared("CC(=O)Nc1ccc(O)cc1")) == 2);
    CHECK(hits("[#6]~[#7]", prepared("CCO")) == 0);
}

TEST_CASE("SMARTS ring primitives", "[chem][smarts]") {
    CHECK(hits("[r5]", prepared("C1CCCC1")) == 5);
    CHECK(hits("[r5]", prepared("C1CCCCC1")) == 0);
    CHECK(hits("[R2]", prepared("c1ccc2ccccc2c1")) == 2);  // naphthalene fusion atoms
    CHECK(hits("[R2]", prepared("c1ccccc1")) == 0);
    CHECK(hits("[!R]", prepared("Cc1ccccc1")) == 1);       // the methyl only
    CHECK(hits("[!R]", prepared("c1ccccc1")) == 0);
    CHECK(hits("[!R]", prepared("CCO")) == 3);
    // Caffeine's fused 5/6 system: both fusion atoms are in two rings.
    CHECK(hits("[r5;R2]", prepared("Cn1cnc2c1c(=O)n(C)c(=O)n2C")) == 2);
}

TEST_CASE("SMARTS charge, degree and bond expressions", "[chem][smarts]") {
    const auto acetate = prepared("CC(=O)[O-]");
    CHECK(hits("[O-]", acetate) == 1);
    CHECK(hits("[OX1-1]", acetate) == 1);
    CHECK(hits("[CD3]", acetate) == 1);
    CHECK(hits("C-,=O", acetate) == 2);   // either bond order
    CHECK(hits("C!=O", acetate) == 1);    // the single-bonded oxygen only
    const auto toluene = prepared("c1ccccc1C");
    CHECK(hits("c@c", toluene) == 6);     // ring bonds
    CHECK(hits("c!@C", toluene) == 1);    // the exocyclic bond
}

TEST_CASE("A malformed SMARTS is an error that names itself", "[chem][smarts]") {
    for (const std::string bad : {"[C", "C(", "[C&]", "[^3]", "[13C]", "[Q]", "C1CC", ""}) {
        std::string err;
        const auto p = parseSmarts(bad, &err);
        CHECK_FALSE(p.has_value());
        CHECK_FALSE(err.empty());
    }
    // Unsupported primitives are refused by name rather than silently ignored,
    // because a silently dropped primitive matches the wrong atoms.
    std::string err;
    CHECK_FALSE(parseSmarts("[C^3]", &err).has_value());
    CHECK(err.find("hybridisation") != std::string::npos);
}

TEST_CASE("SMARTS atom maps are parsed and retained", "[chem][smarts]") {
    std::string err;
    const auto p = parseSmarts("[C:1][O:2]", &err);
    REQUIRE(p.has_value());
    REQUIRE(p->atoms.size() == 2);
    CHECK(p->atoms[0].map == 1);
    CHECK(p->atoms[1].map == 2);
}

TEST_CASE("matchAll honours its limit and matches() short-circuits", "[chem][smarts]") {
    const auto naph = prepared("c1ccc2ccccc2c1");
    const auto p = parseSmarts("c", nullptr);
    REQUIRE(p.has_value());
    CHECK(matchAll(*p, naph.mol, naph.rings, 3).size() == 3);
    CHECK(matchAll(*p, naph.mol, naph.rings, 100).size() == 10);
    CHECK(matches(*p, naph.mol, naph.rings));
    const auto none = parseSmarts("[Br]", nullptr);
    REQUIRE(none.has_value());
    CHECK_FALSE(matches(*none, naph.mol, naph.rings));
}
