// tests/test_chem_solubility.cpp - buffer capacity, pH-solubility, BCS numbers and
// the dissolution/precipitation time course.
//
// What these cases defend, and why each one is here:
//  * the Van Slyke headline fixture - 0.1 M buffer at pH == pKa is 0.0576 mol/L/pH -
//    and the fact that pure water's buffer value is the nonzero water-only term, not
//    a zero somebody optimized away;
//  * the General Solubility Equation is the published expression exactly, carries an
//    error bar (a Predicted number without one is not acceptable), and is NEVER run
//    without a melting point;
//  * the salt-limited branch has exactly one kink, at the analytic intersection - a
//    kink drawn in the wrong place misstates the most important feature of the plot;
//  * without a Ksp there is no salt branch at all and pHmax names what is missing;
//  * dissolution conserves solid + dissolved + precipitated mass to well under
//    1e-9 mg, because a model that loses mass is wrong in a way one curve hides;
//  * kppt is supplied or fitted and never predicted: asking for precipitation
//    without one is a NotComputed, not an assumed rate.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <string>
#include <vector>

#include "chem/Solubility.h"

using namespace biocad;
using namespace biocad::chem;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

SolubilityInput measuredAcid() {
    SolubilityInput in;
    in.moleculeId = "acid";
    in.kind = IonizationKind::MonoproticAcid;
    in.pKa = 4.5;
    in.hasPKa = true;
    in.measuredS0Molar = 1.0e-4;
    in.hasMeasuredS0 = true;
    in.measuredS0Source = "measured intrinsic solubility, user input";
    in.pHMin = 1.0;
    in.pHMax = 10.0;
    in.pHStep = 0.001;
    return in;
}

DissolutionInput baseDissolution() {
    DissolutionInput in;
    in.doseMg = 100.0;
    in.molWeight = 180.16;
    in.volumeL = 0.250;
    in.initialRadiusUm = 25.0;
    in.densityGPerCm3 = 1.3;
    in.diffusivityCm2PerS = 5.0e-6;
    in.diffusionLayerUm = 25.0;
    in.solubilityMolar = 5.0e-3;
    in.horizonS = 3600.0;
    in.stepS = 0.5;
    return in;
}

}  // namespace

TEST_CASE("Van Slyke buffer capacity", "[chem][solubility]") {
    SECTION("a 0.1 M buffer at pH == pKa is 0.0576 mol/L per pH unit") {
        const std::vector<BufferComponent> comp{{"acetate", 4.75, 0.1}};
        REQUIRE_THAT(vanSlykeBeta(comp, 4.75), WithinAbs(0.0576, 5e-5));
    }
    SECTION("pure water at pH 7 is exactly the water-only term") {
        const double expected = 2.302585092994046 * (1.0e-14 / 1.0e-7 + 1.0e-7);
        REQUIRE_THAT(vanSlykeBeta({}, 7.0), WithinRel(expected, 1e-15));
        REQUIRE(vanSlykeBeta({}, 7.0) > 0.0);
    }
    SECTION("the maximum sits at the pKa and is grid-resolved") {
        BufferSpec spec;
        spec.components = {{"acetate", 4.75, 0.1}};
        spec.pHMin = 2.0;
        spec.pHMax = 8.0;
        spec.pHStep = 0.01;
        const BufferReport r = bufferCapacity(spec);
        REQUIRE(r.curve.size() == 601);
        REQUIRE_THAT(r.maxCapacityPh.value, WithinAbs(4.75, 0.011));
        REQUIRE_THAT(r.maxCapacity.value, WithinAbs(0.0576, 5e-5));
        REQUIRE(r.maxCapacityPh.error == 0.5 * spec.pHStep);
        REQUIRE(r.betaAtPh74.provenance == Provenance::Model);
        REQUIRE(r.betaAtPh74.unit == "mol/L/pH");
        REQUIRE_FALSE(r.assumptions.empty());
    }
}

TEST_CASE("General Solubility Equation", "[chem][solubility]") {
    const Quantity q = gseIntrinsicSolubility(2.0, 150.0);
    REQUIRE_THAT(std::log10(q.value), WithinAbs(-2.75, 1e-12));
    REQUIRE(q.provenance == Provenance::Predicted);
    REQUIRE(q.unit == "mol/L");
    // A Predicted quantity without an error bar is not acceptable, and the source
    // must state the benchmark the bar came from.
    REQUIRE(q.error > 0.0);
    REQUIRE(q.source.find("log10 units") != std::string::npos);

    SECTION("no melting point means no GSE") {
        SolubilityInput in;
        in.kind = IonizationKind::MonoproticAcid;
        in.pKa = 4.5;
        in.hasPKa = true;
        in.logP = 2.0;
        in.hasLogP = true;
        const SolubilityReport r = phSolubility(in);
        REQUIRE(r.intrinsic.provenance == Provenance::NotComputed);
        REQUIRE(r.intrinsic.source == "melting point");
        REQUIRE(r.curve.empty());
        REQUIRE(r.solubilityAtPh74.provenance == Provenance::NotComputed);
    }

    SECTION("the GSE refuses to extrapolate outside its fitted domain") {
        // beta-alanine: logP -3.05, MP 200 C. The unguarded equation gives 63 mol/L
        // (5.6 kg/L), which is not a solubility.
        const Quantity out = gseIntrinsicSolubility(-3.05, 200.0);
        REQUIRE(out.provenance == Provenance::NotComputed);
        REQUIRE(out.source.find("fitted domain") != std::string::npos);

        SolubilityInput in;
        in.kind = IonizationKind::MonoproticAcid;
        in.pKa = 3.55;                 in.hasPKa = true;
        in.logP = -3.05;               in.hasLogP = true;
        in.meltingPointC = 200.0;      in.hasMeltingPoint = true;
        const SolubilityReport r = phSolubility(in);
        REQUIRE(r.intrinsic.provenance == Provenance::NotComputed);
        REQUIRE(r.curve.empty());
        REQUIRE(r.solubilityAtPh74.provenance == Provenance::NotComputed);
        REQUIRE_FALSE(r.warnings.empty());
    }
}

TEST_CASE("pH-solubility of a monoprotic acid", "[chem][solubility]") {
    const SolubilityInput in = measuredAcid();
    const SolubilityReport r = phSolubility(in);
    REQUIRE(r.intrinsic.provenance == Provenance::Measured);
    // S(pH)/S0 = 1 + Ka/[H+] = 1 + 10^(pH - pKa).
    REQUIRE_THAT(r.solubilityAtPh74.value / 1.0e-4,
                 WithinAbs(1.0 + std::pow(10.0, 2.9), 1e-9));

    SECTION("without a Ksp there is no salt branch and pHmax says why") {
        REQUIRE(r.pHmax.provenance == Provenance::NotComputed);
        REQUIRE(r.pHmax.source == "salt solubility product");
        for (const auto& p : r.curve) REQUIRE_FALSE(p.saltLimited);
    }
}

TEST_CASE("the salt-limited branch has one kink at the analytic pHmax", "[chem][solubility]") {
    SolubilityInput in = measuredAcid();
    in.ksp = 1.0e-5;
    in.hasKsp = true;
    in.counterionMolar = 0.15;   // common-ion suppression is an explicit input
    in.hasCounterion = true;
    const SolubilityReport r = phSolubility(in);

    const double ka = std::pow(10.0, -in.pKa);
    const double analytic = -std::log10(in.measuredS0Molar * ka * in.counterionMolar / in.ksp);
    REQUIRE(r.pHmax.provenance == Provenance::Model);
    REQUIRE_THAT(r.pHmax.value, WithinAbs(analytic, 1e-9));

    int transitions = 0;
    for (std::size_t i = 1; i < r.curve.size(); ++i)
        if (r.curve[i].saltLimited != r.curve[i - 1].saltLimited) ++transitions;
    REQUIRE(transitions == 1);

    const double plateau = in.measuredS0Molar + in.ksp / in.counterionMolar;
    REQUIRE_THAT(std::pow(10.0, r.curve.back().logS), WithinRel(plateau, 1e-12));

    SECTION("a Ksp without a counterion concentration cannot place the plateau") {
        SolubilityInput noCn = in;
        noCn.hasCounterion = false;
        const SolubilityReport nr = phSolubility(noCn);
        REQUIRE(nr.pHmax.provenance == Provenance::NotComputed);
        REQUIRE(nr.pHmax.source == "counterion concentration");
    }
}

TEST_CASE("BCS dose, dissolution and absorption numbers", "[chem][solubility]") {
    SolubilityInput in = measuredAcid();
    in.doseMg = 100.0;              in.hasDose = true;
    in.molWeight = 180.16;          in.hasMolWeight = true;
    in.particleRadiusUm = 25.0;     in.hasParticleRadius = true;
    in.diffusivityCm2PerS = 5.0e-6; in.hasDiffusivity = true;
    in.densityGPerCm3 = 1.3;        in.hasDensity = true;
    in.peffCmPerS = 2.0e-4;         in.hasPeff = true;
    const SolubilityReport r = phSolubility(in);

    const double s = r.solubilityAtPh74.value;
    const double rCm = 25.0e-4;
    REQUIRE_THAT(r.doseNumber.value,
                 WithinRel(((100.0 / 1000.0 / 180.16) / 0.250) / s, 1e-12));
    REQUIRE_THAT(r.dissolutionNumber.value,
                 WithinRel((3.0 * 5.0e-6 / (rCm * rCm * 1.3)) * (s * 180.16 / 1000.0) * 3600.0,
                           1e-12));
    REQUIRE_THAT(r.absorptionNumber.value, WithinRel(2.0e-4 * 3600.0 / 1.0, 1e-12));
    // Dimensionless, so no unit string may be attached.
    REQUIRE(r.doseNumber.unit.empty());
    REQUIRE(r.dissolutionNumber.unit.empty());
    REQUIRE(r.absorptionNumber.unit.empty());
    // The residence time and radius are assumptions and must be stated as numbers.
    bool stated = false;
    for (const auto& a : r.assumptions)
        if (a.find("3600.0 s") != std::string::npos && a.find("1.00 cm") != std::string::npos)
            stated = true;
    REQUIRE(stated);

    SECTION("each number is NotComputed when its own input is missing") {
        SolubilityInput bare = measuredAcid();
        const SolubilityReport br = phSolubility(bare);
        REQUIRE(br.doseNumber.source == "dose");
        REQUIRE(br.dissolutionNumber.source == "particle radius");
        REQUIRE(br.absorptionNumber.source == "effective permeability Peff");
    }
}

TEST_CASE("dissolution conserves mass and reports when 85% is unreachable",
          "[chem][solubility]") {
    const DissolutionReport r = dissolutionTimeCourse(baseDissolution());
    REQUIRE(r.points.size() == 7201);
    REQUIRE(r.maxMassImbalance < 1e-9);
    REQUIRE(r.timeTo85Pct.provenance == Provenance::Model);
    REQUIRE(r.timeTo85Pct.value > 0.0);
    REQUIRE(r.timeTo85Pct.value < 3600.0);
    // Hixson-Crowell: the particle shrinks, it does not grow.
    REQUIRE(r.points.back().particleRadiusUm < r.points.front().particleRadiusUm);

    SECTION("a solubility-limited dose never reaches 85% and says so") {
        DissolutionInput in = baseDissolution();
        in.solubilityMolar = 1.0e-5;
        const DissolutionReport lr = dissolutionTimeCourse(in);
        REQUIRE(lr.timeTo85Pct.provenance == Provenance::NotComputed);
        REQUIRE(lr.timeTo85Pct.source.find("85%") != std::string::npos);
        REQUIRE(lr.maxMassImbalance < 1e-9);
    }
}

TEST_CASE("precipitation needs a rate constant it is never allowed to invent",
          "[chem][solubility]") {
    DissolutionInput in = baseDissolution();
    in.precipitation = true;
    in.precipSolubilityMolar = 5.0e-4;

    const DissolutionReport missing = dissolutionTimeCourse(in);
    REQUIRE(missing.points.empty());
    REQUIRE(missing.timeTo85Pct.provenance == Provenance::NotComputed);
    REQUIRE(missing.timeTo85Pct.source.find("kppt") != std::string::npos);

    in.kpptPerS = 2.0e-3;
    in.hasKppt = true;
    const DissolutionReport r = dissolutionTimeCourse(in);
    REQUIRE(r.maxMassImbalance < 1e-9);
    REQUIRE(r.points.back().precipitatedMg > 0.0);
    // The solution relaxes toward S from above; it must not undershoot it.
    REQUIRE(r.points.back().dissolvedMolar >= 5.0e-4 - 1e-9);

    SECTION("kppt is recovered exactly from a clean decay") {
        std::vector<double> t, c;
        const double k = 3.5e-3, s = 5.0e-4, c0 = 4.0e-3;
        for (int i = 0; i <= 20; ++i) {
            t.push_back(i * 30.0);
            c.push_back(s + (c0 - s) * std::exp(-k * t.back()));
        }
        const Quantity q = fitPrecipitationRate(t, c, s);
        REQUIRE(q.provenance == Provenance::Model);
        REQUIRE(q.unit == "1/s");
        REQUIRE_THAT(q.value, WithinRel(k, 1e-12));
        REQUIRE(fitPrecipitationRate({}, {}, s).provenance == Provenance::NotComputed);
    }
}
