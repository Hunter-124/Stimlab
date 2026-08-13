// tests/test_chem_speciation.cpp - the component-tableau equilibrium solver and
// the independent-group microspecies ladder.
//
// What these cases defend, and why each one is here:
//  * published fixtures the solver must land on with no per-case tuning: pure
//    water pH 7.000, 0.1 M acetic acid pH 2.881, 0.1 M NH4Cl pH 5.13, glycine
//    pI 5.97;
//  * 1e-8 M strong acid, which is pH 6.978 and not 8 - the case that proves
//    water is in the tableau as chemistry rather than as a post-hoc correction;
//  * the same acid solved two ways (totals route and charge-balance route) must
//    agree, because they are the same equilibrium seen from two sides;
//  * the residual is the contract: max relative mass balance stays far below the
//    1e-10 acceptance bound at every pH of a diprotic titration;
//  * fractions sum to 1 - a distribution diagram whose bars do not is lying;
//  * Davies is refused, not extrapolated, above I = 0.5 M;
//  * a NotComputed pKa yields NotComputed dependents naming the group, never a
//    curve computed from an assumed pKa.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

#include "chem/Speciation.h"

using namespace biocad;
using biocad::chem::solveSpeciation;
using biocad::chem::solveSpeciationPh;
using biocad::chem::titrateGroups;
using Catch::Matchers::WithinAbs;

namespace {

SpeciationProblem monoproticAcid(double total, double pKa) {
    SpeciationProblem p;
    p.components = {"H", "A"};
    p.totals = {total, total};   // the acid was added protonated
    p.species = {"HA", "A-"};
    p.stoichiometry = {{1, 1}, {0, 1}};
    p.logK = {pKa, 0.0};
    p.charges = {0.0, -1.0};
    return p;
}

SpeciationProblem ammoniumChloride(double total, double pKa) {
    SpeciationProblem p;
    p.components = {"H", "N", "Cl"};
    p.totals = {total, total, total};
    p.species = {"NH4+", "NH3", "Cl-"};
    p.stoichiometry = {{1, 1, 0}, {0, 1, 0}, {0, 0, 1}};
    p.logK = {pKa, 0.0, 0.0};
    p.charges = {1.0, 0.0, -1.0};
    return p;
}

SpeciationProblem diproticAcid(double total, double pKa1, double pKa2) {
    SpeciationProblem p;
    p.components = {"H", "L"};
    p.totals = {0.0, total};
    p.species = {"H2L", "HL-", "L2-"};
    p.stoichiometry = {{2, 1}, {1, 1}, {0, 1}};
    p.logK = {pKa1 + pKa2, pKa2, 0.0};
    p.charges = {0.0, -1.0, -2.0};
    return p;
}

Quantity measuredPKa(double v) {
    return makeQuantity(v, "pKa", 0.0, Provenance::Measured, "cited experimental value");
}

}  // namespace

TEST_CASE("pure water is pH 7.000 at 25 C", "[chem][speciation]") {
    SpeciationProblem p;
    p.components = {"H"};
    p.totals = {0.0};   // water only: H+ and OH- are added by the solver

    const SpeciationResult r = solveSpeciation(p);
    REQUIRE(r.converged);
    REQUIRE_THAT(r.pH, WithinAbs(7.0, 1e-6));
    REQUIRE(r.massBalanceResidual < 1e-10);
    REQUIRE(r.concentrations.size() == 2);   // H+ and OH- were appended
}

TEST_CASE("1e-8 M strong acid is pH 6.978, not 8", "[chem][speciation]") {
    SpeciationProblem p;
    p.components = {"H", "Cl"};
    p.totals = {1e-8, 1e-8};
    p.species = {"Cl-"};
    p.stoichiometry = {{0, 1}};
    p.logK = {0.0};
    p.charges = {-1.0};

    const SpeciationResult r = solveSpeciation(p);
    const double exactH = 0.5 * (1e-8 + std::sqrt(1e-16 + 4e-14));
    REQUIRE(r.converged);
    REQUIRE(r.iterations > 0);
    REQUIRE_THAT(r.pH, WithinAbs(-std::log10(exactH), 1e-9));
    REQUIRE(r.massBalanceResidual < 1e-10);
}

TEST_CASE("0.1 M acetic acid is pH 2.881 by both routes", "[chem][speciation]") {
    const SpeciationResult byTotals = solveSpeciation(monoproticAcid(0.1, 4.756));
    const SpeciationResult byCharge = solveSpeciationPh(monoproticAcid(0.1, 4.756));

    REQUIRE(byTotals.converged);
    REQUIRE(byCharge.converged);
    REQUIRE_THAT(byTotals.pH, WithinAbs(2.881, 0.002));
    REQUIRE_THAT(byCharge.pH, WithinAbs(byTotals.pH, 1e-6));
    REQUIRE(byTotals.massBalanceResidual < 1e-10);
    REQUIRE_THAT(byTotals.fractions[0] + byTotals.fractions[1], WithinAbs(1.0, 1e-12));
}

TEST_CASE("0.1 M NH4Cl is pH 5.13 with chloride as a component", "[chem][speciation]") {
    const SpeciationResult r = solveSpeciationPh(ammoniumChloride(0.1, 9.25));
    REQUIRE(r.converged);
    REQUIRE_THAT(r.pH, WithinAbs(5.13, 0.02));
    REQUIRE(r.massBalanceResidual < 1e-10);
    REQUIRE(r.chargeBalanceResidual < 1e-10);
}

TEST_CASE("solvePh refuses a tableau with no proton component", "[chem][speciation]") {
    // Without a proton component there is nothing to solve pH on, and water's
    // self-ionization cannot be attached either. Saying so beats reporting 0.
    SpeciationProblem p;
    p.components = {"L"};
    p.totals = {0.1};
    p.species = {"L-"};
    p.stoichiometry = {{1}};
    p.logK = {0.0};
    p.charges = {-1.0};

    const SpeciationResult r = solveSpeciationPh(p);
    REQUIRE_FALSE(r.converged);
    REQUIRE_FALSE(r.warnings.empty());

    const SpeciationResult direct = solveSpeciation(p);
    REQUIRE_FALSE(direct.warnings.empty());   // water was not silently invented
}

TEST_CASE("a diprotic titration keeps its mass balance at every pH", "[chem][speciation]") {
    double worstMass = 0.0;
    double worstFractionSum = 0.0;
    for (int i = 0; i < 40; ++i) {
        SpeciationProblem p = diproticAcid(0.05, 2.34, 9.60);
        p.fixedComponent = 0;
        p.fixedLog10Activity = -(0.5 + 0.35 * i);
        const SpeciationResult r = solveSpeciation(p);
        REQUIRE(r.converged);
        worstMass = std::max(worstMass, r.massBalanceResidual);
        const double sum = r.fractions[0] + r.fractions[1] + r.fractions[2];
        worstFractionSum = std::max(worstFractionSum, std::fabs(sum - 1.0));
    }
    REQUIRE(worstMass < 1e-10);
    REQUIRE(worstFractionSum < 1e-12);
}

TEST_CASE("Davies is used inside its domain and refused outside it", "[chem][speciation]") {
    SpeciationProblem inside = monoproticAcid(0.1, 4.756);
    inside.daviesActivities = true;
    inside.ionicStrength = 0.1;
    const SpeciationResult ok = solveSpeciation(inside);
    const double s = std::sqrt(0.1);
    const double expectedGamma =
        std::pow(10.0, -0.5085 * (s / (1.0 + s) - 0.3 * 0.1));
    REQUIRE(ok.converged);
    REQUIRE_THAT(ok.ionicStrength, WithinAbs(0.1, 1e-15));
    REQUIRE_THAT(ok.activityCoefficients[1], WithinAbs(expectedGamma, 1e-12));
    REQUIRE(ok.massBalanceResidual < 1e-10);

    SpeciationProblem outside = inside;
    outside.ionicStrength = 0.6;
    const SpeciationResult refused = solveSpeciation(outside);
    REQUIRE_FALSE(refused.converged);
    REQUIRE(refused.concentrations.empty());
    REQUIRE_FALSE(refused.warnings.empty());
    REQUIRE(refused.warnings[0].find("0.5") != std::string::npos);
}

TEST_CASE("glycine's isoelectric point is 5.97", "[chem][speciation]") {
    const std::vector<IonizableGroup> groups = {
        {"carboxyl", measuredPKa(2.34), true},
        {"amine", measuredPKa(9.60), false},
    };
    const Quantity logP =
        makeQuantity(-3.21, "log10 P", 0.0, Provenance::Measured, "experimental logP");
    const SpeciationCurve c = titrateGroups(groups, logP, 0.0, 14.0, 0.1);

    REQUIRE(c.labels.size() == 4);   // two groups, four microstates
    REQUIRE(c.points.size() == 141);
    REQUIRE(c.isoelectricPoint.provenance == Provenance::Measured);
    REQUIRE_THAT(c.isoelectricPoint.value, WithinAbs(5.97, 0.01));
    for (const SpeciationPoint& pt : c.points) {
        double sum = 0.0;
        for (double f : pt.microspeciesFractions) sum += f;
        REQUIRE_THAT(sum, WithinAbs(1.0, 1e-12));
    }
}

TEST_CASE("logD of a base is logP plus log10 of its neutral fraction", "[chem][speciation]") {
    const Quantity logP =
        makeQuantity(2.0, "log10 P", 0.0, Provenance::Measured, "experimental logP");
    const std::vector<IonizableGroup> base = {{"amine", measuredPKa(9.0), false}};
    const SpeciationCurve c = titrateGroups(base, logP, 0.0, 14.0, 0.1);

    const double fNeutral = 1.0 / (1.0 + std::pow(10.0, 9.0 - 7.4));
    REQUIRE_THAT(c.logDAtPh74.value, WithinAbs(2.0 + std::log10(fNeutral), 1e-12));
    REQUIRE_THAT(c.logDAtPh74.value, WithinAbs(0.389226, 1e-6));

    // The mirror case: an acid of the same pKa is barely ionized at pH 7.4, so
    // its logD sits 0.0107 log unit below logP.
    const std::vector<IonizableGroup> acid = {{"phenol", measuredPKa(9.0), true}};
    const SpeciationCurve ac = titrateGroups(acid, logP, 0.0, 14.0, 0.1);
    REQUIRE_THAT(logP.value - ac.logDAtPh74.value, WithinAbs(0.010774, 1e-6));
}

TEST_CASE("a missing pKa yields NotComputed naming the group", "[chem][speciation]") {
    const std::vector<IonizableGroup> groups = {
        {"carboxyl", measuredPKa(4.2), true},
        {"secondary amine", notComputed("pKa"), false},
    };
    const Quantity logP =
        makeQuantity(1.5, "log10 P", 0.0, Provenance::Measured, "experimental logP");
    const SpeciationCurve c = titrateGroups(groups, logP, 0.0, 14.0, 0.5);

    REQUIRE(c.points.empty());
    REQUIRE(c.isoelectricPoint.provenance == Provenance::NotComputed);
    REQUIRE(c.logDAtPh74.provenance == Provenance::NotComputed);
    REQUIRE(c.isoelectricPoint.source.find("secondary amine") != std::string::npos);
}

TEST_CASE("logD is NotComputed without a logP", "[chem][speciation]") {
    const std::vector<IonizableGroup> groups = {{"carboxyl", measuredPKa(4.2), true}};
    const SpeciationCurve c = titrateGroups(groups, notComputed("logP"), 0.0, 14.0, 1.0);
    REQUIRE_FALSE(c.points.empty());
    REQUIRE(c.logDAtPh74.provenance == Provenance::NotComputed);
}
