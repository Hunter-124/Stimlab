// PK engine: the integrator is checked against the closed forms it exists to replace.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

#include "pkpd/PkEngine.h"

using namespace biocad;

namespace {

Quantity measured(double v, const char* unit) {
    return makeQuantity(v, unit, 0.0, Provenance::Measured, "test fixture");
}

}  // namespace

TEST_CASE("IV bolus AUC and half-life match the closed form", "[pkpd]") {
    PkModelSpec spec;
    spec.model = PkModel::IvBolus;
    spec.clearance = measured(5.0, "L/h");
    spec.volume = measured(50.0, "L");
    spec.unboundFraction = measured(1.0, "");
    spec.stepH = 0.005;
    spec.horizonH = 240.0;   // ~35 half-lives: the trapezoidal AUC has converged
    DoseRegimen regimen;
    regimen.doses.push_back(DoseEvent{0.0, 100.0, 0.0});

    const PkProfile p = pkpd::simulate(spec, regimen);
    REQUIRE(p.auc.provenance == Provenance::Model);   // a simulation is a constructed artefact
    REQUIRE(p.auc.unit == "mg*h/L");
    REQUIRE_THAT(p.auc.value, Catch::Matchers::WithinAbs(20.0, 1e-6));
    REQUIRE_THAT(p.halfLife.value, Catch::Matchers::WithinAbs(0.6931471805599453 * 10.0, 1e-9));
    REQUIRE_THAT(p.cmax.value, Catch::Matchers::WithinAbs(2.0, 1e-9));
    REQUIRE(p.tmax.value == 0.0);
    REQUIRE(p.accumulation.provenance == Provenance::NotComputed);
    REQUIRE(p.accumulation.source.find("repeated regimen") != std::string::npos);
}

TEST_CASE("Oral one-compartment RK4 matches the Bateman function", "[pkpd]") {
    PkModelSpec spec;
    spec.model = PkModel::OralOneCompartment;
    spec.bioavailability = measured(0.8, "");
    spec.absorptionRate = measured(1.2, "1/h");
    spec.clearance = measured(7.5, "L/h");   // ke = CL/V = 0.15 /h
    spec.volume = measured(50.0, "L");
    spec.unboundFraction = measured(1.0, "");
    spec.stepH = 0.001;
    spec.horizonH = 24.0;
    DoseRegimen regimen;
    regimen.doses.push_back(DoseEvent{0.0, 100.0, 0.0});

    const PkProfile p = pkpd::simulate(spec, regimen);
    for (double t : {0.5, 1.0, 2.0, 6.0, 24.0}) {
        const double expected = pkpd::batemanConcentration(100.0, 0.8, 1.2, 0.15, 50.0, t);
        std::size_t best = 0;
        for (std::size_t i = 1; i < p.timeH.size(); ++i) {
            if (std::fabs(p.timeH[i] - t) < std::fabs(p.timeH[best] - t)) best = i;
        }
        REQUIRE_THAT(p.timeH[best], Catch::Matchers::WithinAbs(t, 1e-9));
        REQUIRE_THAT(p.concentrationMgPerL[best], Catch::Matchers::WithinAbs(expected, 1e-6));
    }
    REQUIRE_FALSE(p.flipFlop);
}

TEST_CASE("The ka == ke Bateman limit is finite and correct", "[pkpd]") {
    const double ke = 0.5, t = 3.0;
    const double limit = pkpd::batemanConcentration(100.0, 1.0, ke, ke, 20.0, t);
    REQUIRE(std::isfinite(limit));
    REQUIRE_THAT(limit, Catch::Matchers::WithinAbs(100.0 * ke * t * std::exp(-ke * t) / 20.0,
                                                   1e-12));
    // The limit is approached continuously from either side.
    const double near = pkpd::batemanConcentration(100.0, 1.0, ke + 1e-7, ke, 20.0, t);
    REQUIRE_THAT(near, Catch::Matchers::WithinAbs(limit, 1e-6));
}

TEST_CASE("Accumulation ratio over 10 q12h doses matches the closed form", "[pkpd]") {
    PkModelSpec spec;
    spec.model = PkModel::IvBolus;
    spec.clearance = measured(5.0, "L/h");
    spec.volume = measured(50.0, "L");
    spec.unboundFraction = measured(1.0, "");
    spec.stepH = 0.01;
    spec.horizonH = 120.0;
    DoseRegimen regimen;
    for (int i = 0; i < 10; ++i) regimen.doses.push_back(DoseEvent{12.0 * i, 100.0, 0.0});

    const PkProfile p = pkpd::simulate(spec, regimen);
    const double ke = 0.1;
    REQUIRE_THAT(p.accumulation.value,
                 Catch::Matchers::WithinAbs(1.0 / (1.0 - std::exp(-ke * 12.0)), 1e-12));
    REQUIRE(p.accumulation.unit.empty());
    REQUIRE_THAT(pkpd::accumulationRatio(ke, 12.0),
                 Catch::Matchers::WithinAbs(p.accumulation.value, 1e-15));
    REQUIRE_THAT(pkpd::steadyStateAverage(1.0, 100.0, 5.0, 12.0),
                 Catch::Matchers::WithinAbs(100.0 / 60.0, 1e-12));
}

TEST_CASE("Michaelis-Menten elimination has no single half-life", "[pkpd]") {
    PkModelSpec spec;
    spec.model = PkModel::IvBolus;
    spec.clearance = measured(5.0, "L/h");
    spec.volume = measured(50.0, "L");
    spec.unboundFraction = measured(1.0, "");
    spec.vmax = measured(10.0, "mg/h");
    spec.km = measured(1.0, "mg/L");
    spec.stepH = 0.01;
    spec.horizonH = 48.0;
    DoseRegimen regimen;
    regimen.doses.push_back(DoseEvent{0.0, 100.0, 0.0});

    const PkProfile p = pkpd::simulate(spec, regimen);
    REQUIRE(p.halfLife.provenance == Provenance::NotComputed);
    REQUIRE(p.halfLife.source == "linear elimination");
    REQUIRE(p.auc.value > 0.0);
}

TEST_CASE("Occupancy refuses a missing Kd and computes one from an unbound series",
          "[pkpd]") {
    PkModelSpec spec;
    spec.model = PkModel::IvBolus;
    spec.clearance = measured(5.0, "L/h");
    spec.volume = measured(50.0, "L");
    spec.unboundFraction = measured(0.5, "");
    spec.stepH = 0.01;
    spec.horizonH = 24.0;
    DoseRegimen regimen;
    regimen.doses.push_back(DoseEvent{0.0, 100.0, 0.0});
    const PkProfile p = pkpd::simulate(spec, regimen);

    const OccupancyCurve none = pkpd::occupancy(p, notComputed("a measured Kd"));
    REQUIRE(none.occupancy.empty());
    REQUIRE(none.peakOccupancy.provenance == Provenance::NotComputed);

    // Cu(0) = 0.5 * 2.0 = 1.0 mg/L; with Kd = 1.0 mg/L the peak occupancy is exactly 0.5.
    const OccupancyCurve curve = pkpd::occupancy(p, measured(1.0, "mg/L"));
    REQUIRE(curve.occupancy.size() == p.timeH.size());
    REQUIRE_THAT(curve.peakOccupancy.value, Catch::Matchers::WithinAbs(0.5, 1e-12));
    REQUIRE(curve.timeAbove50Pct.value >= 0.0);
    REQUIRE(curve.timeAbove50Pct.value < 1e-6);   // it starts at 0.5 and only falls
}

TEST_CASE("Assumptions name every non-measured parameter", "[pkpd]") {
    PkModelSpec spec;
    spec.model = PkModel::OralOneCompartment;
    spec.bioavailability = makeQuantity(0.8, "", 0.0, Provenance::Heuristic, "stated default");
    spec.absorptionRate = makeQuantity(1.2, "1/h", 0.0, Provenance::Model, "stated default");
    spec.clearance = makeQuantity(7.5, "L/h", 0.0, Provenance::Predicted, "PKSmart");
    spec.volume = measured(50.0, "L");
    spec.unboundFraction = makeQuantity(0.3, "", 0.0, Provenance::Heuristic, "stated default");
    spec.stepH = 0.01;
    spec.horizonH = 24.0;
    DoseRegimen regimen;
    regimen.doses.push_back(DoseEvent{0.0, 100.0, 0.0});

    const PkProfile p = pkpd::simulate(spec, regimen);
    auto mentions = [&](const char* needle) {
        for (const auto& a : p.assumptions) {
            if (a.find(needle) != std::string::npos) return true;
        }
        return false;
    };
    REQUIRE(mentions("F = 0.80"));
    REQUIRE(mentions("ka = 1.200"));
    REQUIRE(mentions("fu = 0.300"));
    REQUIRE(mentions("10.4%"));            // the PKSmart generalisation figure
    REQUIRE(mentions("not a dose recommendation"));
    REQUIRE_FALSE(mentions("V = 50"));     // measured parameters are not assumptions
}
