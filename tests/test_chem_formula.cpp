// tests/test_chem_formula.cpp - exact formula arithmetic.
//
// What these cases defend, and why each one is here:
//  * the masses are the NIST table's masses - water and caffeine are checked
//    against published values, so a corrupted digit in isotopes.json fails here
//    rather than shifting every mass in the app by a quiet amount;
//  * the parser's four awkward cases: nesting, hydrate dots, isotope labels and
//    a trailing charge, each of which is silently wrong in the obvious
//    implementation;
//  * the envelope is a real convolution - bromine's 1:1 doublet and glucose's
//    M+1 both come out of the shipped abundances, not out of a fitted constant;
//  * an ion's m/z accounts for the electron mass, and a neutral's m/z is
//    NotComputed rather than a number nobody should read;
//  * balancing is exact: a right equation balances, and two kinds of wrong
//    equation report why instead of being approximated.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstdlib>
#include <map>
#include <string>
#include <utility>

#include "chem/Formula.h"

using namespace biocad;
using Catch::Matchers::WithinAbs;

namespace {

// The isotope pack is data on disk. ctest runs from the build directory, so
// point the loader at the in-tree copy; only a successful load is cached.
void useInTreePack() {
    const std::string dir = std::string(BIOCAD_ASSETS_DIR) + "/packs/descriptors";
#if defined(_WIN32)
    _putenv_s("BIOCAD_DESCRIPTOR_DIR", dir.c_str());
#else
    setenv("BIOCAD_DESCRIPTOR_DIR", dir.c_str(), 1);
#endif
}

chem::ParsedFormula parse(const std::string& text) {
    useInTreePack();
    auto f = chem::parseFormula(text);
    REQUIRE(f.has_value());
    return *f;
}

int countOf(const chem::ParsedFormula& f, int z, int massNumber = 0) {
    for (const auto& t : f.terms)
        if (t.z == z && t.massNumber == massNumber) return t.count;
    return 0;
}

}  // namespace

TEST_CASE("The NIST isotope pack loads", "[chem][formula]") {
    useInTreePack();
    INFO(chem::isotopeTableNote());
    REQUIRE(chem::isotopeTableOk());
    CHECK(std::string(chem::formulaCitation()).find("NIST") != std::string::npos);
}

TEST_CASE("Monoisotopic and average masses are the published values", "[chem][formula]") {
    const auto water = parse("H2O");
    // 15.9949146196 + 2 x 1.00782503223, and 15.999 + 2 x 1.008.
    CHECK_THAT(chem::monoisotopicMass(water), WithinAbs(18.0105646863, 1e-7));
    CHECK_THAT(chem::averageMass(water), WithinAbs(18.015, 1e-9));

    const auto caffeine = parse("C8H10N4O2");
    CHECK_THAT(chem::monoisotopicMass(caffeine), WithinAbs(194.080376, 1e-5));
    CHECK_THAT(chem::averageMass(caffeine), WithinAbs(194.19, 0.01));
    CHECK(chem::hillFormula(caffeine) == "C8H10N4O2");
    // 1 + 8 - 10/2 + 4/2: the three rings and three double bonds of caffeine.
    CHECK(chem::ringPlusDoubleBondEquivalents(caffeine, nullptr) == 6.0);
    CHECK(chem::electronCount(caffeine) == 102);
}

TEST_CASE("The parser handles nesting, hydrates, labels and charge", "[chem][formula]") {
    CHECK(chem::hillFormula(parse("C6H5(CH3)")) == "C7H8");
    CHECK(chem::hillFormula(parse("K[Al(SO4)2].12H2O")) == "H24AlKO20S2");

    const auto blue = parse("CuSO4.5H2O");
    CHECK(countOf(blue, 29) == 1);
    CHECK(countOf(blue, 16) == 1);
    CHECK(countOf(blue, 8) == 9);
    CHECK(countOf(blue, 1) == 10);
    CHECK(chem::hillFormula(blue) == "H10CuO9S");

    // A leading digit run multiplies the segment when it cannot be an isotope
    // (there is no 2C), and reads as a label when it can (2H is deuterium).
    CHECK(chem::hillFormula(parse("2C6H6")) == "C12H12");
    CHECK(chem::hillFormula(parse("2H2O")) == "[2H]2O");

    const auto labelled = parse("[13C]");
    CHECK_THAT(chem::monoisotopicMass(labelled) - chem::monoisotopicMass(parse("C")),
               WithinAbs(1.00336, 1e-5));
    CHECK(chem::hillFormula(labelled) == "[13C]");
    // A labelled position is one isotope, so it is not averaged.
    const auto heavyGlucose = parse("13C6H12O6");
    CHECK(countOf(heavyGlucose, 6, 13) == 6);
    CHECK_THAT(chem::averageMass(heavyGlucose) - chem::averageMass(parse("C6H12O6")),
               WithinAbs(6 * (13.00335483507 - 12.011), 1e-9));

    CHECK(parse("NH4+").charge == 1);
    CHECK(countOf(parse("NH4+"), 1) == 4);   // the 4 is an atom count, not a charge
    CHECK(parse("SO4 2-").charge == -2);
    CHECK(parse("SO4^2-").charge == -2);
    CHECK(parse("Ca+2").charge == 2);
    CHECK(parse("Mg++").charge == 2);
    CHECK(chem::hillFormula(parse("SO4 2-")) == "O4S2-");
}

TEST_CASE("FormulaMass carries measured masses and an honest m/z", "[chem][formula]") {
    useInTreePack();
    const auto neutral = chem::toFormulaMass("H2O");
    CHECK(neutral.monoisotopic.provenance == Provenance::Measured);
    CHECK(neutral.monoisotopic.unit == "Da");
    CHECK(neutral.average.provenance == Provenance::Measured);
    CHECK(neutral.mz.provenance == Provenance::NotComputed);
    CHECK(neutral.mz.source == "charge is zero");

    const auto ion = chem::toFormulaMass("NH4+");
    REQUIRE(ion.mz.provenance == Provenance::Measured);
    // An ion is the neutral minus one electron, so m/z is not simply M.
    CHECK_THAT(ion.mz.value, WithinAbs(ion.monoisotopic.value - 0.00054857990907, 1e-12));
    CHECK(ion.electrons == 10);

    const auto nonsense = chem::toFormulaMass("!!");
    CHECK(nonsense.monoisotopic.provenance == Provenance::NotComputed);
    CHECK_FALSE(nonsense.warnings.empty());
}

TEST_CASE("Isotope envelopes come out of the shipped abundances", "[chem][formula]") {
    const auto glucose = parse("C6H12O6");
    const auto env = chem::isotopeEnvelope(glucose, 1e-8);
    REQUIRE(env.peaks.size() >= 3);
    CHECK(env.peaks[0].intensity == 1.0);
    CHECK(env.peaks[0].nominalShift == 0);
    CHECK_THAT(env.peaks[0].mass, WithinAbs(chem::monoisotopicMass(glucose), 1e-9));
    for (std::size_t i = 1; i < env.peaks.size(); ++i)
        CHECK(env.peaks[i].mass > env.peaks[i - 1].mass);
    // 6 x 13C/12C + 12 x 2H/1H + 6 x 17O/16O from isotopes.json = 0.068560.
    const double expected = 6 * 0.0107 / 0.9893 + 12 * 0.000115 / 0.999885
                            + 6 * 0.00038 / 0.99757;
    CHECK_THAT(env.peaks[1].intensity, WithinAbs(expected, 2e-4));
    CHECK_THAT(env.peaks[1].intensity, WithinAbs(0.068560, 1e-5));

    // Bromine's near-1:1 doublet is the classic check that the convolution is a
    // convolution: M and M+2 only, with no M+1.
    const auto br2 = chem::isotopeEnvelope(parse("Br2"), 1e-6);
    REQUIRE(br2.peaks.size() == 3);
    CHECK(br2.peaks[0].nominalShift == 0);
    CHECK(br2.peaks[1].nominalShift == 2);
    CHECK(br2.peaks[2].nominalShift == 4);
    CHECK_THAT(br2.peaks[1].intensity, WithinAbs(1.0, 1e-12));
    CHECK_THAT(br2.peaks[0].intensity, WithinAbs(0.5140, 1e-3));
}

TEST_CASE("Formula finding uses only the bounds it was given", "[chem][formula]") {
    useInTreePack();
    const std::map<std::string, std::pair<int, int>> bounds{
        {"C", {1, 10}}, {"H", {0, 20}}, {"N", {0, 5}}, {"O", {0, 5}}};
    const auto hits = chem::findFormulas(194.080376, 0.003, bounds);
    REQUIRE_FALSE(hits.empty());
    // Best match first, and every hit inside the tolerance.
    CHECK(chem::hillFormula(hits.front()) == "C8H10N4O2");
    for (const auto& h : hits)
        CHECK(std::abs(chem::monoisotopicMass(h) - 194.080376) <= 0.003);
    // No element is guessed: without bounds there are no candidates at all.
    CHECK(chem::findFormulas(194.080376, 0.003, {}).empty());
}

TEST_CASE("Equations balance exactly or not at all", "[chem][formula]") {
    useInTreePack();
    const auto combustion = chem::balanceEquation({"C3H8", "O2"}, {"CO2", "H2O"});
    REQUIRE(combustion.balanced);
    CHECK(combustion.reactantCoefficients == std::vector<int>{1, 5});
    CHECK(combustion.productCoefficients == std::vector<int>{3, 4});
    CHECK(combustion.atomEconomy.provenance == Provenance::Measured);
    CHECK_THAT(combustion.atomEconomy.value, WithinAbs(64.6915, 1e-3));
    CHECK(combustion.theoreticalYield.provenance == Provenance::NotComputed);

    const auto rust = chem::balanceEquation({"Fe", "O2"}, {"Fe2O3"});
    REQUIRE(rust.balanced);
    CHECK(rust.reactantCoefficients == std::vector<int>{4, 3});
    CHECK(rust.productCoefficients == std::vector<int>{2});

    // Amounts give a limiting reagent and a theoretical yield of the first
    // product - arithmetic on the user's own numbers, nothing more.
    const auto withAmounts =
        chem::balanceEquation({"C3H8", "O2"}, {"CO2", "H2O"}, {44.097, 100.0});
    REQUIRE(withAmounts.balanced);
    CHECK(withAmounts.limitingReagent == "O2");
    CHECK(withAmounts.theoreticalYield.provenance == Provenance::Measured);
    CHECK(withAmounts.theoreticalYield.unit == "g");
    CHECK_THAT(withAmounts.theoreticalYield.value, WithinAbs(82.522, 0.01));

    const auto impossible = chem::balanceEquation({"H2O"}, {"CO2"});
    CHECK_FALSE(impossible.balanced);
    CHECK_FALSE(impossible.warnings.empty());
    CHECK(impossible.atomEconomy.provenance == Provenance::NotComputed);
    CHECK(impossible.reactantCoefficients.empty());

    // Two independent reactions written as one: the null space is 2-dimensional,
    // so there is no unique answer and none is invented.
    const auto ambiguous = chem::balanceEquation({"C", "O2", "H2"}, {"CO2", "H2O"});
    CHECK_FALSE(ambiguous.balanced);
    REQUIRE_FALSE(ambiguous.warnings.empty());
    CHECK(ambiguous.warnings.front().find("null space") != std::string::npos);
}
