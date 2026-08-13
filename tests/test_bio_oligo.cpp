// Nearest-neighbour oligo thermodynamics.
//
// The three things that actually go wrong in a Tm implementation are all tested
// here: the x factor in the Tm equation (1 for a self-complementary strand, 4
// otherwise), the salt correction being applied to entropy rather than to the
// temperature, and a degenerate base being averaged instead of refused. The
// numeric anchor is a published worked example - Biopython's MeltingTemp module
// documents Tm_NN("CGTTCCAAAGATGTGGGCATGAGCTTAC") = 60.32 degC at 50 mM Na+ with
// 25 nM of each strand, using the same Allawi/SantaLucia unified parameters - so
// this file checks against a number computed elsewhere, not against itself.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <filesystem>
#include <string>

#include "bio/NucSeq.h"
#include "bio/OligoThermo.h"

using namespace biocad;
using namespace biocad::bio;

namespace {

const NnParameters& params() {
    static const NnParameters p = loadNnParameters(
        std::filesystem::path(BIOCAD_ASSETS_DIR) / "packs" / "nucleic" / "nn-thermodynamics.json");
    return p;
}

ThermoOptions standardConditions() {
    ThermoOptions o;
    o.naMolar = 0.05;      // 50 mM monovalent
    o.oligoMolar = 5e-8;   // TOTAL strand concentration: 25 nM + 25 nM
    return o;
}

}  // namespace

TEST_CASE("the NN pack carries the ten unique parameters and its own constants", "[bio][oligo]") {
    const NnParameters& p = params();
    REQUIRE(p.step('A', 'A') != nullptr);
    REQUIRE(p.step('A', 'A')->dH == -7.9);
    REQUIRE(p.step('A', 'A')->dS == -22.2);
    // The complement identity: reading a duplex 5'->3' on either strand must give
    // the same numbers, which is what makes ten parameters cover sixteen steps.
    REQUIRE(p.step('T', 'T')->dH == p.step('A', 'A')->dH);
    REQUIRE(p.step('C', 'G')->dH == -10.6);
    REQUIRE(p.step('G', 'C')->dH == -9.8);
    REQUIRE(p.step('T', 'G')->dH == p.step('C', 'A')->dH);
    REQUIRE(p.initiation('G').dH == 0.1);
    REQUIRE(p.initiation('A').dH == 2.3);   // the terminal-AT penalty
    REQUIRE(p.symmetry().dS == -1.4);
    REQUIRE(p.gasConstant() > 1.98);
    REQUIRE(p.saltCoefficient() == 0.368);
}

TEST_CASE("Tm reproduces a published worked example", "[bio][oligo]") {
    const OligoThermo t = oligoThermo("CGTTCCAAAGATGTGGGCATGAGCTTAC", params(),
                                      standardConditions());
    // The dH/dS decomposition is asserted alongside the Tm, so a future failure
    // shows whether the parameters or the equation moved.
    REQUIRE_THAT(t.deltaH.value, Catch::Matchers::WithinAbs(-222.9, 1e-9));
    REQUIRE_THAT(t.deltaS.value, Catch::Matchers::WithinAbs(-602.5, 1e-9));
    REQUIRE_THAT(t.tm.value, Catch::Matchers::WithinAbs(60.32, 0.5));
    REQUIRE(t.deltaH.unit == "kcal/mol");
    REQUIRE(t.deltaS.unit == "cal/(mol*K)");
    REQUIRE(t.tm.unit == "degC");
    REQUIRE(t.tm.provenance == Provenance::Predicted);
    REQUIRE(t.tm.source.find("SantaLucia") != std::string::npos);
    // dG37 is the standard-state value derived from the same dH and dS.
    REQUIRE_THAT(t.deltaG37.value,
                 Catch::Matchers::WithinAbs(t.deltaH.value - 310.15 * t.deltaS.value / 1000.0,
                                            1e-12));
}

TEST_CASE("the x factor follows self-complementarity", "[bio][oligo]") {
    const NnParameters& p = params();
    const ThermoOptions o = standardConditions();
    REQUIRE(isSelfComplementary("CGCGAATTCGCG"));
    REQUIRE_FALSE(isSelfComplementary("CGCGAATTCGCA"));

    const OligoThermo selfComp = oligoThermo("CGCGAATTCGCG", p, o);
    const double dSsalt = selfComp.deltaS.value +
                          p.saltCoefficient() * 11.0 * std::log(o.naMolar);
    const double tmX1 = selfComp.deltaH.value * 1000.0 /
                            (dSsalt + p.gasConstant() * std::log(o.oligoMolar / 1.0)) - 273.15;
    const double tmX4 = selfComp.deltaH.value * 1000.0 /
                            (dSsalt + p.gasConstant() * std::log(o.oligoMolar / 4.0)) - 273.15;
    // The x = 1 branch must be the one used, and the two branches must differ by
    // enough that getting it wrong is not a rounding matter.
    REQUIRE_THAT(selfComp.tm.value, Catch::Matchers::WithinAbs(tmX1, 1e-9));
    REQUIRE(std::abs(tmX1 - tmX4) > 1.0);
    // The symmetry correction is applied to entropy only for the palindrome.
    const OligoThermo nonSelf = oligoThermo("CGCGAATTCGCA", p, o);
    REQUIRE(nonSelf.tm.value < selfComp.tm.value);
}

TEST_CASE("salt, oligo concentration and unused ions are all reported", "[bio][oligo]") {
    const NnParameters& p = params();
    ThermoOptions low = standardConditions();
    ThermoOptions high = standardConditions();
    high.naMolar = 1.0;
    const double tmLow = oligoThermo("CGTTCCAAAGATGTGGGCATGAGCTTAC", p, low).tm.value;
    const double tmHigh = oligoThermo("CGTTCCAAAGATGTGGGCATGAGCTTAC", p, high).tm.value;
    REQUIRE(tmHigh > tmLow);   // more monovalent salt stabilises the duplex

    ThermoOptions dilute = standardConditions();
    dilute.oligoMolar = 5e-10;
    REQUIRE(oligoThermo("CGTTCCAAAGATGTGGGCATGAGCTTAC", p, dilute).tm.value < tmLow);

    // Mg2+ and dNTPs are echoed and explicitly unused: no equivalence is invented.
    ThermoOptions withMg = standardConditions();
    withMg.mgMolar = 0.0015;
    withMg.dntpMolar = 0.0006;
    const OligoThermo t = oligoThermo("CGTTCCAAAGATGTGGGCATGAGCTTAC", p, withMg);
    REQUIRE_THAT(t.tm.value, Catch::Matchers::WithinAbs(tmLow, 1e-12));
    REQUIRE(t.mgMolar == 0.0015);
    bool saidUnused = false;
    for (const auto& a : t.assumptions) {
        if (a.find("UNUSED") != std::string::npos) saidUnused = true;
    }
    REQUIRE(saidUnused);
}

TEST_CASE("a degenerate or impossible input is notComputed, not averaged", "[bio][oligo]") {
    const NnParameters& p = params();
    const OligoThermo degenerate = oligoThermo("ACGTNACGT", p, standardConditions());
    REQUIRE(degenerate.tm.provenance == Provenance::NotComputed);
    REQUIRE(degenerate.deltaH.provenance == Provenance::NotComputed);
    REQUIRE(degenerate.tm.source.find('N') != std::string::npos);

    REQUIRE(oligoThermo("A", p, standardConditions()).tm.provenance == Provenance::NotComputed);
    ThermoOptions zero = standardConditions();
    zero.oligoMolar = 0.0;
    REQUIRE(oligoThermo("ACGTACGT", p, zero).tm.provenance == Provenance::NotComputed);
}

TEST_CASE("structure scans find the obvious structures and label their convention",
          "[bio][oligo]") {
    const NnParameters& p = params();
    const std::string hairpinSeq = "GCGCGGTTTTCCGCGC";
    const auto hp = hairpins(hairpinSeq, p, -1.0);
    REQUIRE_FALSE(hp.empty());
    REQUIRE(hp[0].kind == "hairpin");
    REQUIRE(hp[0].deltaG37.value < -1.0);
    REQUIRE(hp[0].deltaG37.unit == "kcal/mol");
    REQUIRE(hp[0].deltaG37.source.find("loop-closure") != std::string::npos);
    REQUIRE_FALSE(hp[0].alignment.empty());
    // Sorted most stable first.
    for (std::size_t i = 1; i < hp.size(); ++i) {
        REQUIRE(hp[i - 1].deltaG37.value <= hp[i].deltaG37.value);
    }

    const auto self = selfDimers(hairpinSeq, p, -1.0);
    REQUIRE_FALSE(self.empty());
    REQUIRE(self[0].kind == "self-dimer");

    const auto hetero = heteroDimers(hairpinSeq, "GCGCGGAAAACCGCGC", p, -1.0);
    REQUIRE_FALSE(hetero.empty());
    REQUIRE(hetero[0].kind == "hetero-dimer");
    // The perfect register against the exact reverse complement IS the whole duplex,
    // so the scan's helix dG37 must equal the duplex dG37 to the bit. That ties the
    // scan and the Tm path to one parameter table instead of asserting a vague
    // "more stable than" relation, which is not even true in general: a GC-rich
    // partial duplex beats an AT-rich perfect one.
    const std::string probe = "ATGGCCATTGTAATGG";
    const auto perfect = heteroDimers(probe, reverseComplement(probe), p, -1.0);
    REQUIRE_FALSE(perfect.empty());
    REQUIRE_THAT(perfect[0].deltaG37.value,
                 Catch::Matchers::WithinAbs(
                     oligoThermo(probe, p, standardConditions()).deltaG37.value, 1e-9));

    // A sequence with nothing to pair has no structures at a stability cutoff.
    REQUIRE(hairpins("AAAAAAAAAAAAAAAA", p, -1.0).empty());
    // Ambiguity codes are refused here too rather than scored.
    REQUIRE(selfDimers("ACGTNACGTNACGT", p, -1.0).empty());
}
