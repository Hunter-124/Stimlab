// Phase 13.3 - mechanistic drug interactions, and 13.1 - the population layer.
//
// Two properties matter more than any single number here: a missing fm must stop
// the AUC ratio rather than default it, and the dynamic enzyme model's steady state
// must EQUAL the static model, because two implementations of the same mechanism
// that disagree mean at least one of them is wrong.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

#include "core/Physiology.h"
#include "sim/Ddi.h"
#include "sim/Population.h"

using namespace biocad;
using Catch::Matchers::WithinAbs;

namespace {

PerpetratorSpec inhibitor() {
    PerpetratorSpec p;
    p.label = "reversible + time-dependent CYP3A4 inhibitor";
    p.enzyme = "CYP3A4";
    p.ki = 0.5;
    p.kinact = 0.6;
    p.kI = 1.2;
    p.unboundHepaticInletUM = 2.0;
    p.unboundSystemicUM = 1.0;
    p.enterocyteUM = 20.0;
    p.source = "test fixture, not a real compound";
    return p;
}

VictimSpec victim() {
    VictimSpec v;
    v.label = "victim";
    v.fractionMetabolizedByEnzyme = 0.9;
    v.intestinalAvailability = 0.5;
    v.fractionExcretedUnchanged = 0.2;
    return v;
}

PkModelSpec ivSpec() {
    PkModelSpec m;
    m.model = PkModel::IvBolus;
    m.clearance = makeQuantity(5.0, "L/h", 0, Provenance::Measured, "fixture");
    m.volume = makeQuantity(50.0, "L", 0, Provenance::Measured, "fixture");
    m.bioavailability = makeQuantity(1.0, "", 0, Provenance::Measured, "IV");
    m.unboundFraction = makeQuantity(1.0, "", 0, Provenance::Measured, "fixture");
    m.stepH = 0.05;
    m.horizonH = 48.0;
    return m;
}

}  // namespace

TEST_CASE("The physiology pack is the only source of Qh, Qen and kdeg", "[sim][ddi]") {
    const auto& p = core::physiology();
    REQUIRE(p.loaded);
    REQUIRE_THAT(p.hepaticBloodFlowLPerH, WithinAbs(97.0, 1e-12));
    REQUIRE_THAT(p.enterocyteBloodFlowLPerH, WithinAbs(18.0, 1e-12));
    REQUIRE_THAT(core::enzymeDegradationRate("CYP3A4_hepatic"), WithinAbs(0.0193, 1e-12));
    // An enzyme the pack does not carry yields -1, never a plausible-looking default.
    REQUIRE(core::enzymeDegradationRate("CYP9Z9_hepatic") < 0.0);
}

TEST_CASE("FDA basic-model R-values are the published expressions", "[sim][ddi]") {
    const InteractionReport r = sim::interaction(inhibitor(), victim());
    REQUIRE_THAT(r.r1.value, WithinAbs(1.0 + 2.0 / 0.5, 1e-12));
    REQUIRE_THAT(r.r1Gut.value, WithinAbs(1.0 + 20.0 / 0.5, 1e-12));
    const double kdeg = core::enzymeDegradationRate("CYP3A4_hepatic");
    const double kobs = 0.6 * 1.0 / (1.2 + 1.0);   // at the systemic concentration
    REQUIRE_THAT(r.r2.value, WithinAbs((kdeg + kobs) / kdeg, 1e-12));
    REQUIRE(r.rInduction.provenance == Provenance::NotComputed);   // no Emax supplied
}

TEST_CASE("A missing fm makes the AUC ratio NotComputed", "[sim][ddi]") {
    VictimSpec v = victim();
    v.fractionMetabolizedByEnzyme = -1.0;
    const InteractionReport r = sim::interaction(inhibitor(), v);
    REQUIRE(r.aucRatio.provenance == Provenance::NotComputed);
    REQUIRE(r.aucRatio.source == "fm");
    REQUIRE(r.aucRatioHepaticOnly.provenance == Provenance::NotComputed);
    REQUIRE(r.theoreticalCeiling.provenance == Provenance::NotComputed);
}

TEST_CASE("A missing Fg reports the hepatic term only", "[sim][ddi]") {
    VictimSpec v = victim();
    v.intestinalAvailability = -1.0;
    const InteractionReport r = sim::interaction(inhibitor(), v);
    REQUIRE_FALSE(r.gutIncluded);
    REQUIRE(r.aucRatio.provenance == Provenance::NotComputed);
    REQUIRE(r.aucRatioHepaticOnly.provenance == Provenance::Predicted);
    REQUIRE_THAT(r.theoreticalCeiling.value, WithinAbs(10.0, 1e-12));
    // The hepatic ratio can approach the 1/(1-fm) ceiling but never pass it.
    REQUIRE(r.aucRatioHepaticOnly.value < r.theoreticalCeiling.value);
}

TEST_CASE("The dynamic enzyme model's steady state equals the static model",
          "[sim][ddi]") {
    const EnzymeTimeCourse tc = sim::enzymeTimeCourse(inhibitor(), 72.0);
    REQUIRE(tc.kdegUsed > 0.0);
    REQUIRE_THAT(tc.relativeActivity.front(), WithinAbs(1.0, 1e-12));
    REQUIRE(tc.agreement < 1e-9);
    REQUIRE_THAT(tc.steadyStateActivity, WithinAbs(tc.staticModelActivity, 1e-9));

    // ... and the static hepatic factor A*B*C implied by the AUC ratio is the same
    // number again, which is what stops the two code paths from drifting apart.
    const InteractionReport r = sim::interaction(inhibitor(), victim());
    const double abcFromTimeCourse = tc.staticModelActivity / (1.0 + 2.0 / 0.5);
    const double fm = 0.9;
    const double abcFromRatio = (1.0 / r.aucRatioHepaticOnly.value - (1.0 - fm)) / fm;
    REQUIRE_THAT(abcFromTimeCourse, WithinAbs(abcFromRatio, 1e-9));
}

TEST_CASE("Impairment is an exposure ratio, never a dose", "[sim][ddi]") {
    const ImpairmentScenario renal = sim::impairment(victim(), 0.3, 1.0);
    REQUIRE_THAT(renal.exposureRatio.value, WithinAbs(1.0 / (0.2 * 0.3 + 0.8), 1e-12));
    REQUIRE(renal.boundaryStatement.find("not a dose") != std::string::npos);

    const ImpairmentScenario hepatic = sim::impairment(victim(), 1.0, 0.5);
    REQUIRE_THAT(hepatic.exposureRatio.value, WithinAbs(1.0 / (0.2 + 0.8 * 0.5), 1e-12));

    VictimSpec noFe = victim();
    noFe.fractionExcretedUnchanged = -1.0;
    REQUIRE(sim::impairment(noFe, 0.3, 1.0).exposureRatio.provenance ==
            Provenance::NotComputed);
}

TEST_CASE("A population band is reproducible from its seed", "[sim][population]") {
    VariabilitySpec v;
    v.betweenSubject = true;
    v.parameterUncertainty = true;
    v.residualError = true;
    v.parameters = {"CL", "V"};
    v.omega = {0.09, 0.02, 0.02, 0.04};
    v.parameterCovariance = {0.01, 0.004, 0.004, 0.01};
    v.proportionalResidualCv = 0.1;
    v.seed = 20260813;
    v.subjects = 120;
    v.sampler = "latin-hypercube";

    DoseRegimen regimen;
    regimen.doses.push_back(DoseEvent{0.0, 100.0, 0.0});

    const PopulationProfile a = sim::simulatePopulation(ivSpec(), regimen, v);
    const PopulationProfile b = sim::simulatePopulation(ivSpec(), regimen, v);
    VariabilitySpec other = v;
    other.seed = 20260814;
    const PopulationProfile c = sim::simulatePopulation(ivSpec(), regimen, other);

    REQUIRE(a.bands.size() == b.bands.size());
    for (std::size_t i = 0; i < a.bands.size(); ++i) {
        REQUIRE(a.bands[i].p5 == b.bands[i].p5);     // bit equality
        REQUIRE(a.bands[i].p50 == b.bands[i].p50);
        REQUIRE(a.bands[i].p95 == b.bands[i].p95);
        REQUIRE(a.bands[i].p5 <= a.bands[i].p50);
        REQUIRE(a.bands[i].p50 <= a.bands[i].p95);
    }
    bool differs = false;
    for (std::size_t i = 0; i < a.bands.size(); ++i)
        if (a.bands[i].p50 != c.bands[i].p50) differs = true;
    REQUIRE(differs);

    REQUIRE(static_cast<int>(a.sampleTrajectories.size()) == sim::kMaxStoredTrajectories);
    REQUIRE(a.provenanceStatement.find("not a prediction about any individual") !=
            std::string::npos);
}

TEST_CASE("With every layer off the band collapses onto the typical profile",
          "[sim][population]") {
    VariabilitySpec v;
    v.parameters = {"CL"};
    v.seed = 1;
    v.subjects = 8;
    DoseRegimen regimen;
    regimen.doses.push_back(DoseEvent{0.0, 100.0, 0.0});
    const PopulationProfile p = sim::simulatePopulation(ivSpec(), regimen, v);
    for (const auto& b : p.bands) REQUIRE_THAT(b.p95 - b.p5, WithinAbs(0.0, 1e-12));
    // 100 mg into 50 L cleared at 5 L/h: AUC over 48 h is (D/CL)(1 - exp(-ke*T)).
    REQUIRE_THAT(p.medianAuc.value, WithinAbs(20.0 * (1.0 - std::exp(-4.8)), 1e-3));
    REQUIRE(p.medianAuc.provenance == Provenance::Model);
}
