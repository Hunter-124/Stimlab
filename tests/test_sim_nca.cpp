// Phase 13.2 - noncompartmental analysis.
//
// The fixtures are monoexponential on purpose: every NCA parameter then has a
// closed form, so the test can assert the ARITHMETIC rather than agreement with a
// previous run. IV bolus D = 100 mg with CL = 5 L/h and V = 50 L gives ke = 0.1 /h,
// C0 = 2 mg/L, AUCinf = D/CL = 20 mg*h/L, MRT = 1/ke = 10 h and Vss = CL*MRT = 50 L.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <string>

#include "sim/Nca.h"

using namespace biocad;
using Catch::Matchers::WithinAbs;

namespace {

ConcentrationSeries monoexponential(double lastH, int points) {
    ConcentrationSeries s;
    s.subjectId = "iv-bolus";
    s.dose = 100.0;
    s.intravenous = true;
    for (int i = 0; i <= points; ++i) {
        const double t = lastH * i / points;
        s.timeH.push_back(t);
        s.concentration.push_back(2.0 * std::exp(-0.1 * t));
    }
    return s;
}

}  // namespace

TEST_CASE("IV bolus NCA recovers the closed-form parameters", "[sim][nca]") {
    const NcaResult r = sim::noncompartmental(monoexponential(48.0, 96));

    REQUIRE_THAT(r.cmax.value, WithinAbs(2.0, 1e-12));
    REQUIRE_THAT(r.tmax.value, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(r.lambdaZ.value, WithinAbs(0.1, 1e-10));
    REQUIRE_THAT(r.halfLife.value, WithinAbs(std::log(2.0) / 0.1, 1e-9));
    REQUIRE_THAT(r.aucInfinity.value, WithinAbs(20.0, 1e-9));
    REQUIRE_THAT(r.clearance.value, WithinAbs(5.0, 1e-9));
    REQUIRE_THAT(r.volumeZ.value, WithinAbs(50.0, 1e-7));
    REQUIRE_THAT(r.meanResidenceTime.value, WithinAbs(10.0, 1e-7));
    REQUIRE_THAT(r.volumeSteadyState.value, WithinAbs(50.0, 1e-7));
    REQUIRE(r.lambdaZPointCount >= 3);
    REQUIRE_FALSE(r.extrapolationUnreliable);
}

TEST_CASE("The log-down trapezoid is exact on monoexponential data", "[sim][nca]") {
    // A LINEAR trapezoid over a decaying exponential always overestimates; the
    // log-down rule is the exact integral of the segment, so AUClast must equal
    // (C0/ke)(1 - exp(-ke*T)) to machine precision at any sampling density.
    for (int points : {12, 48, 96}) {
        const NcaResult r = sim::noncompartmental(monoexponential(48.0, points));
        const double exact = (2.0 / 0.1) * (1.0 - std::exp(-0.1 * 48.0));
        REQUIRE_THAT(r.aucLast.value, WithinAbs(exact, 1e-11));
    }
}

TEST_CASE("Truncating at one half-life extrapolates exactly 50%", "[sim][nca]") {
    const double halfLife = std::log(2.0) / 0.1;
    const NcaResult r = sim::noncompartmental(monoexponential(halfLife, 20));

    REQUIRE_THAT(r.percentExtrapolated.value, WithinAbs(50.0, 1e-9));
    REQUIRE(r.extrapolationUnreliable);
    // Every quantity that inherits the extrapolation must be NAMED, not merely
    // implied by a single global flag.
    for (const char* derived : {"AUCinf", "CL (or CL/F)", "Vz (or Vz/F)", "Vss", "MRT", "AUMC"}) {
        bool named = false;
        for (const auto& w : r.warnings)
            if (w.find(derived) != std::string::npos && w.find("unreliable") != std::string::npos)
                named = true;
        REQUIRE(named);
    }
}

TEST_CASE("Vss is IV-only and lambda_z needs three post-Tmax points", "[sim][nca]") {
    ConcentrationSeries oral = monoexponential(48.0, 96);
    oral.intravenous = false;
    const NcaResult r = sim::noncompartmental(oral);
    REQUIRE(r.volumeSteadyState.provenance == Provenance::NotComputed);
    REQUIRE(r.clearance.provenance == Provenance::Measured);   // CL/F is still computable

    ConcentrationSeries few;
    few.dose = 100.0;
    few.intravenous = true;
    few.timeH = {0.0, 1.0, 2.0};
    few.concentration = {2.0, 1.8, 1.6};   // Tmax at t = 0, so only two points follow it
    const NcaResult r2 = sim::noncompartmental(few);
    REQUIRE(r2.lambdaZ.provenance == Provenance::NotComputed);
    REQUIRE(r2.aucInfinity.provenance == Provenance::NotComputed);
    REQUIRE(r2.aucLast.provenance == Provenance::Measured);
}

TEST_CASE("Steady-state metrics come from the interval, not the tail", "[sim][nca]") {
    ConcentrationSeries ss;
    ss.dose = 100.0;
    ss.intravenous = true;
    ss.tauH = 12.0;
    for (int i = 0; i <= 48; ++i) {
        const double t = 12.0 * i / 48.0;
        ss.timeH.push_back(t);
        ss.concentration.push_back(2.0 * std::exp(-0.1 * t) + 0.5);
    }
    const NcaResult r = sim::noncompartmental(ss);
    REQUIRE(r.aucTau.provenance == Provenance::Measured);
    REQUIRE_THAT(r.cAverage.value, WithinAbs(r.aucTau.value / 12.0, 1e-12));
    const double cMin = 2.0 * std::exp(-1.2) + 0.5;
    REQUIRE_THAT(r.swing.value, WithinAbs((2.5 - cMin) / cMin, 1e-9));
    // Clearance at steady state uses AUCtau, so it is dose/AUCtau, not dose/AUCinf.
    REQUIRE_THAT(r.clearance.value, WithinAbs(100.0 / r.aucTau.value, 1e-9));
}
