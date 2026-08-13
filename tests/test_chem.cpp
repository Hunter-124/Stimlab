#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <iterator>
#include <optional>
#include <string>

#include "chem/Analysis.h"
#include "chem/Descriptors.h"
#include "chem/Smiles.h"

using namespace biocad::chem;
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
    // Kekule input, deliberately: tpsa() perceives aromaticity itself, so this is
    // the same 60.26 the aromatic spelling gives. The 1.82 gap to the published
    // 58.44 is the fused-imidazole nitrogen, not a perception bug - which is why
    // the tolerance can now be 2.0 instead of the 8.0 an unperceived 56.22 needed.
    REQUIRE_THAT(tpsa(mustParse("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")), WithinAbs(58.44, 2.0));  // caffeine
}

TEST_CASE("Ring perception", "[chem][rings]") {
    REQUIRE(ringCount(mustParse("CC(N)Cc1ccccc1")) == 1);            // one benzene
    REQUIRE(ringCount(mustParse("CNC(C)Cc1ccc2OCOc2c1")) == 2);      // MDMA: arene + dioxole
    REQUIRE(ringCount(mustParse("CCC")) == 0);
}

// Functional-group perception now lives in one place: the SMARTS rule pack, and
// its full flag-by-flag coverage is in test_chem_groups.cpp. Keeping a second,
// smaller copy of those expectations here would give two answers to the same
// question the first time a rule is corrected.

// Wildman-Crippen logP against experimental/consensus reference logP.
//
// WHAT THIS CASE CAN AND CANNOT ASSERT. crippenLogP() is now a faithful
// implementation of the published additive method (chem/Crippen.h), whose OWN
// reported performance is an RMS error of about 0.67 log units against
// experiment. So a tight per-compound band here would not be testing this
// codebase, it would be asserting that a 1999 additive model is better than its
// authors measured it to be - and the previous +/-0.7 band only held because the
// estimator it replaced was an ad-hoc invention tuned to look good on these ten
// rows. The band below is the method's published error doubled, and the
// aggregate check is stated in the method's own terms (RMS), so a real
// regression - a mis-parsed pack, a lost hydrogen class, a broken aromaticity
// perception - still fails it loudly. Exact-value oracles verified against the
// reference implementation live in test_chem_crippen.cpp.
TEST_CASE("Wildman-Crippen logP tracks reference for the library", "[chem][logp]") {
    struct Case { const char* smiles; const char* name; double ref; };
    const Case cases[] = {
        {"CC(N)Cc1ccccc1", "amphetamine", 1.76},
        {"CNC(C)Cc1ccccc1", "methamphetamine", 2.07},
        {"CNC(C)Cc1ccc2OCOc2c1", "mdma", 2.15},
        {"CC(N)Cc1ccc2OCOc2c1", "mda", 1.64},
        {"NCCc1ccccc1", "phenethylamine", 1.41},
        {"NCCc1ccc(O)cc1", "tyramine", 0.86},
        {"CNC(C)C(O)c1ccccc1", "ephedrine", 1.13},
        {"CC(N)C(=O)c1ccccc1", "cathinone", 0.59},
        {"CC(=O)Nc1ccc(O)cc1", "acetaminophen", 0.46},
        {"CC(N)Cc1ccc(F)cc1", "4-fluoroamphetamine", 1.90},
    };
    double sumSqErr = 0.0;
    for (const auto& c : cases) {
        const auto m = mustParse(c.smiles);
        const double got = crippenLogP(m);
        INFO(c.name << "  ref=" << c.ref << "  got=" << got);
        REQUIRE_THAT(got, WithinAbs(c.ref, 1.35));  // 2x the method's published RMS
        sumSqErr += (got - c.ref) * (got - c.ref);
    }
    const double rms = std::sqrt(sumSqErr / static_cast<double>(std::size(cases)));
    INFO("RMS deviation from experiment = " << rms << " log units");
    REQUIRE(rms < 0.67);  // this set is no worse than the method's own average
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
