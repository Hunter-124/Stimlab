#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <optional>
#include <string>

#include "chem/Analysis.h"
#include "chem/Descriptors.h"
#include "chem/Smiles.h"

using namespace stimlab::chem;
using Catch::Matchers::WithinAbs;

namespace {
Molecule mustParse(const std::string& smi) {
    auto m = parseSmiles(smi);
    REQUIRE(m.has_value());
    return *m;
}
}  // namespace

// The strongest parser/H-model validation: computed Hill formula must match the
// known molecular formula for a broad spread of real drug structures.
TEST_CASE("Molecular formula matches reference for the library", "[chem][formula]") {
    struct Case { const char* smiles; const char* formula; };
    const Case cases[] = {
        {"CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "C8H10N4O2"},        // caffeine
        {"CC(N)Cc1ccccc1", "C9H13N"},                          // amphetamine
        {"CNC(C)Cc1ccccc1", "C10H15N"},                        // methamphetamine
        {"CNC(C)Cc1ccc2OCOc2c1", "C11H15NO2"},                 // MDMA
        {"COC(=O)C(c1ccccc1)C1CCCCN1", "C14H19NO2"},           // methylphenidate
        {"COC(=O)C1C(OC(=O)c2ccccc2)CC3CCC1N3C", "C17H21NO4"}, // cocaine
        {"NCCc1ccc(O)c(O)c1", "C8H11NO2"},                     // dopamine
        {"CN1CCCC1c1cccnc1", "C10H14N2"},                      // nicotine
        {"CC(=O)Nc1ccc(O)cc1", "C8H9NO2"},                     // acetaminophen
        {"Cn1cnc2c1c(=O)[nH]c(=O)n2C", "C7H8N4O2"},            // theobromine
        {"NC(=O)CS(=O)C(c1ccccc1)c1ccccc1", "C15H15NO2S"},     // modafinil
        {"CC(N)Cc1ccc(F)cc1", "C9H12FN"},                      // 4-fluoroamphetamine
        {"NCCc1ccccc1", "C8H11N"},                             // phenethylamine
        {"CC(C)(N)Cc1ccccc1", "C10H15N"},                      // phentermine
    };
    for (const auto& c : cases) {
        const auto m = mustParse(c.smiles);
        INFO(c.smiles);
        REQUIRE(molecularFormula(m) == c.formula);
    }
}

TEST_CASE("Molecular weight matches reference", "[chem][mw]") {
    REQUIRE_THAT(molecularWeight(mustParse("CC(N)Cc1ccccc1")), WithinAbs(135.21, 0.1));
    REQUIRE_THAT(molecularWeight(mustParse("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")), WithinAbs(194.19, 0.1));
    REQUIRE_THAT(molecularWeight(mustParse("NCCc1ccc(O)c(O)c1")), WithinAbs(153.18, 0.1));
}

TEST_CASE("H-bond donor/acceptor counts (Lipinski)", "[chem][lipinski]") {
    const auto amp = mustParse("CC(N)Cc1ccccc1");
    REQUIRE(hbdCount(amp) == 1);
    REQUIRE(hbaCount(amp) == 1);

    const auto dopa = mustParse("NCCc1ccc(O)c(O)c1");
    REQUIRE(hbdCount(dopa) == 3);   // NH2 + 2 phenol OH
    REQUIRE(hbaCount(dopa) == 3);   // 1 N + 2 O

    const auto caf = mustParse("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
    REQUIRE(hbaCount(caf) == 6);    // 4 N + 2 O
}

TEST_CASE("Rotatable bonds (Veber)", "[chem][rotatable]") {
    REQUIRE(rotatableBondCount(mustParse("CC(N)Cc1ccccc1")) == 2);   // amphetamine
    REQUIRE(rotatableBondCount(mustParse("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")) == 0);  // caffeine
}

TEST_CASE("Ertl TPSA matches reference within tolerance", "[chem][tpsa]") {
    REQUIRE_THAT(tpsa(mustParse("CC(N)Cc1ccccc1")), WithinAbs(26.02, 0.5));        // amphetamine
    REQUIRE_THAT(tpsa(mustParse("NCCc1ccc(O)c(O)c1")), WithinAbs(66.48, 0.5));     // dopamine
    REQUIRE_THAT(tpsa(mustParse("CC(=O)Nc1ccc(O)cc1")), WithinAbs(49.33, 0.5));    // acetaminophen
    REQUIRE_THAT(tpsa(mustParse("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")), WithinAbs(58.44, 8.0));  // caffeine
}

TEST_CASE("Ring perception", "[chem][rings]") {
    REQUIRE(ringCount(mustParse("CC(N)Cc1ccccc1")) == 1);            // one benzene
    REQUIRE(ringCount(mustParse("CNC(C)Cc1ccc2OCOc2c1")) == 2);      // MDMA: arene + dioxole
    REQUIRE(ringCount(mustParse("CCC")) == 0);
}

TEST_CASE("Functional-group perception", "[chem][groups]") {
    const auto amp = detectGroups(mustParse("CC(N)Cc1ccccc1"));
    REQUIRE(amp.phenethylamine);
    REQUIRE(amp.primaryAmine);
    REQUIRE(amp.aromaticRing);
    REQUIRE_FALSE(amp.ester);

    REQUIRE(detectGroups(mustParse("COC(=O)C1C(OC(=O)c2ccccc2)CC3CCC1N3C")).ester);     // cocaine
    REQUIRE(detectGroups(mustParse("NCCc1ccc(O)c(O)c1")).catechol);                      // dopamine
    REQUIRE(detectGroups(mustParse("CC(=O)Nc1ccc(O)cc1")).amide);                        // acetaminophen
    REQUIRE(detectGroups(mustParse("CNC(C)Cc1ccc2OCOc2c1")).methylenedioxy);            // MDMA
    REQUIRE(detectGroups(mustParse("CC(N)C(=O)c1ccccc1")).arylKetone);                   // cathinone
    REQUIRE(detectGroups(mustParse("NC(=O)CS(=O)C(c1ccccc1)c1ccccc1")).sulfoxide);       // modafinil
}

TEST_CASE("Morgan/Tanimoto similarity is sane", "[chem][fingerprint]") {
    const auto amp = morganFingerprint(mustParse("CC(N)Cc1ccccc1"));
    const auto meth = morganFingerprint(mustParse("CNC(C)Cc1ccccc1"));
    const auto caf = morganFingerprint(mustParse("CN1C=NC2=C1C(=O)N(C(=O)N2C)C"));

    REQUIRE(tanimoto(amp, amp) == 1.0);
    const double close = tanimoto(amp, meth);
    const double far = tanimoto(amp, caf);
    REQUIRE(close > far);
    REQUIRE(close > 0.3);
    REQUIRE(far < 0.3);
}
