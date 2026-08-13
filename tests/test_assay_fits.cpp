// Tests for assay::Fits - the curve, kinetics and modality fitters.
//
// Every case here is generated from a known truth model, so the assertion is
// "the fitter recovers the parameters that produced the data", not "the fitter
// returns what it returned last time".
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

#include "assay/Fits.h"

using namespace biocad;
using namespace biocad::assay;
using Catch::Matchers::WithinRel;

namespace {

double fourPl(double conc, double a, double d, double c, double b) {
    return d + (a - d) / (1.0 + std::pow(10.0, b * (std::log10(c) - std::log10(conc))));
}

double fivePl(double x, double a, double d, double c, double b, double g) {
    return d + (a - d) / std::pow(1.0 + std::pow(x / c, b), g);
}

}  // namespace

TEST_CASE("4PL recovers its truth parameters from noiseless data", "[assay][fits]") {
    const double a = 100.0, d = 5.0, c = 1e-7, b = 1.2;
    std::vector<DosePoint> pts;
    for (int i = 0; i < 12; ++i) {
        const double conc = 1e-10 * std::pow(10.0, i * 0.5);
        pts.push_back({conc, fourPl(conc, a, d, c, b), -1.0});
    }

    const FitResult f = fitFourParameterLogistic(pts);
    REQUIRE(f.converged);
    CHECK_THAT(f.parameters[0].value.value, WithinRel(a, 1e-9));
    CHECK_THAT(f.parameters[1].value.value, WithinRel(d, 1e-9));
    CHECK_THAT(f.derivedEc50.value, WithinRel(c, 1e-9));
    CHECK_THAT(f.parameters[3].value.value, WithinRel(b, 1e-9));
    CHECK(f.parameters[2].name == "log10EC50");
    CHECK(f.rank == 4);
    CHECK(f.conditionNumber > 0.0);
    CHECK_FALSE(f.extrapolated);
    // A fitted parameter is a constructed artefact, and the response scale is
    // arbitrary, so the asymptotes must carry no unit.
    CHECK(f.parameters[0].value.provenance == Provenance::Model);
    CHECK(f.parameters[0].value.unit.empty());
    CHECK(f.derivedEc50.unit == "mol/L");
}

TEST_CASE("an EC50 outside the tested ladder is flagged extrapolated", "[assay][fits]") {
    // Truth EC50 sits two decades above the highest concentration tested.
    std::vector<DosePoint> pts;
    for (int i = 0; i < 8; ++i) {
        const double conc = 1e-10 * std::pow(10.0, i * 0.25);
        pts.push_back({conc, fourPl(conc, 100.0, 0.0, 1e-6, 1.0), -1.0});
    }
    const FitResult f = fitFourParameterLogistic(pts);
    CHECK(f.extrapolated);
    CHECK_FALSE(f.warnings.empty());
}

TEST_CASE("5PL reports C*(2^(1/G)-1)^(1/B) as the EC50, never C", "[assay][fits]") {
    const double a = 2.0, d = 120.0, c = 1e-7, b = 1.5, g = 0.6;
    std::vector<DosePoint> pts;
    for (int i = 0; i < 12; ++i) {
        const double conc = 1e-10 * std::pow(10.0, i * 0.4);
        pts.push_back({conc, fivePl(conc, a, d, c, b, g), -1.0});
    }

    const FitResult f = fitFiveParameterLogistic(pts);
    REQUIRE(f.converged);
    CHECK_THAT(f.parameters[3].value.value, WithinRel(b, 1e-8));
    CHECK_THAT(f.parameters[4].value.value, WithinRel(g, 1e-8));
    const double fittedC = std::pow(10.0, f.parameters[2].value.value);
    CHECK_THAT(fittedC, WithinRel(c, 1e-8));

    const double analytic = c * std::pow(std::pow(2.0, 1.0 / g) - 1.0, 1.0 / b);
    CHECK_THAT(f.derivedEc50.value, WithinRel(analytic, 1e-12));
    // The whole point: with G != 1 the inflection scale is not the EC50.
    CHECK(std::abs(f.derivedEc50.value - fittedC) / fittedC > 0.10);
}

TEST_CASE("5PL refuses fewer than eight distinct concentrations", "[assay][fits]") {
    std::vector<DosePoint> pts;
    for (int i = 0; i < 7; ++i) {
        const double conc = 1e-9 * std::pow(10.0, i * 0.5);
        pts.push_back({conc, fivePl(conc, 0.0, 100.0, 1e-7, 1.5, 0.6), -1.0});
    }
    const FitResult f = fitFiveParameterLogistic(pts);
    CHECK_FALSE(f.converged);
    CHECK(f.note.find("8 distinct concentrations") != std::string::npos);
    CHECK(f.derivedEc50.provenance == Provenance::NotComputed);
}

TEST_CASE("Michaelis-Menten and substrate inhibition recover their parameters",
          "[assay][fits]") {
    SECTION("Michaelis-Menten") {
        const double vmax = 2.5, km = 5e-5;
        std::vector<KineticPoint> pts;
        for (int i = 0; i < 10; ++i) {
            const double s = 1e-6 * std::pow(10.0, i * 0.3);
            pts.push_back({s, vmax * s / (km + s), -1.0});
        }
        const FitResult f = fitMichaelisMenten(pts);
        REQUIRE(f.converged);
        CHECK_THAT(f.parameters[0].value.value, WithinRel(vmax, 1e-9));
        CHECK_THAT(f.parameters[1].value.value, WithinRel(km, 1e-9));
    }

    SECTION("substrate inhibition peaks at sqrt(Km*Ki)") {
        const double vmax = 3.0, km = 1e-5, ki = 1e-3;
        std::vector<KineticPoint> pts;
        for (int i = 0; i < 14; ++i) {
            const double s = 1e-6 * std::pow(10.0, i * 0.3);
            pts.push_back({s, vmax * s / (km + s * (1.0 + s / ki)), -1.0});
        }
        const FitResult f = fitSubstrateInhibition(pts);
        REQUIRE(f.converged);
        CHECK_THAT(f.parameters[1].value.value, WithinRel(km, 1e-7));
        CHECK_THAT(f.parameters[2].value.value, WithinRel(ki, 1e-7));
        CHECK_THAT(f.derivedKd.value, WithinRel(std::sqrt(km * ki), 1e-6));
    }
}

TEST_CASE("Morrison's quadratic solution departs from the classic approximation",
          "[assay][fits]") {
    const double ki = 1e-8, et = 1e-8;
    // At [E]t = [I] = Ki_app the approximation's [I] >> [E]t premise is false.
    const double quadratic = morrisonFraction(et, ki, ki);
    const double classic = classicInhibitionFraction(ki, ki);
    CHECK_THAT(quadratic, WithinRel(0.6180339887498949, 1e-12));
    CHECK_THAT(classic, WithinRel(0.5, 1e-12));
    CHECK(std::abs(quadratic - classic) / classic > 0.20);

    std::vector<DosePoint> pts;
    for (int i = 0; i < 10; ++i) {
        const double conc = 1e-10 * std::pow(10.0, i * 0.35);
        pts.push_back({conc, morrisonFraction(et, conc, ki), -1.0});
    }
    const FitResult f = fitMorrisonTightBinding(pts, et);
    REQUIRE(f.converged);
    CHECK_THAT(f.parameters[1].value.value, WithinRel(ki, 1e-8));
    CHECK_THAT(f.derivedKd.value, WithinRel(ki, 1e-8));

    // [E]t is measured, never fitted: without it there is no quadratic to solve.
    CHECK_FALSE(fitMorrisonTightBinding(pts, 0.0).converged);
}

TEST_CASE("the global [S]x[I] fit recovers the true modality", "[assay][fits]") {
    const double vmax = 1.0, km = 2e-5, ki = 5e-7;
    const std::vector<double> svals{2e-6, 5e-6, 1e-5, 2e-5, 4e-5, 8e-5, 1.6e-4, 3.2e-4};
    const std::vector<double> ivals{0.0, 1e-7, 2.5e-7, 5e-7, 1e-6, 2e-6};

    SECTION("competitive") {
        std::vector<InhibitionPoint> m;
        for (double s : svals)
            for (double i : ivals)
                m.push_back({s, i, vmax * s / (km * (1.0 + i / ki) + s), -1.0});
        const ModelComparison mc = fitInhibitionModality(m);
        REQUIRE(mc.candidates.size() == 4);
        CHECK(mc.decisive);
        CHECK(mc.candidates.front().modality == InhibitionModality::Competitive);
        // Candidates come back ascending in AICc so the runner-up is auditable.
        CHECK(mc.candidates[0].aicc <= mc.candidates[1].aicc);
    }

    SECTION("uncompetitive") {
        std::vector<InhibitionPoint> m;
        for (double s : svals)
            for (double i : ivals)
                m.push_back({s, i, vmax * s / (km + s * (1.0 + i / ki)), -1.0});
        const ModelComparison mc = fitInhibitionModality(m);
        CHECK(mc.decisive);
        CHECK(mc.candidates.front().modality == InhibitionModality::Uncompetitive);
    }

    SECTION("a single [S] = Km cannot distinguish modality") {
        // Competitive and uncompetitive rate laws are algebraically identical at
        // [S] = Km, so the AICc difference is zero and the answer is Unknown.
        std::vector<InhibitionPoint> m;
        for (int rep = 0; rep < 8; ++rep)
            for (double i : ivals)
                m.push_back({km, i, vmax * km / (km * (1.0 + i / ki) + km), -1.0});
        const ModelComparison mc = fitInhibitionModality(m);
        CHECK_FALSE(mc.decisive);
        CHECK(mc.deltaAicc < 2.0);
        CHECK(mc.candidates.front().modality == InhibitionModality::Unknown);
        CHECK(mc.conclusion.find("unknown") != std::string::npos);
    }
}

TEST_CASE("inverse prediction returns the EC50 at the midpoint response", "[assay][fits]") {
    const double a = 100.0, d = 0.0, c = 1e-7, b = 1.0;
    std::vector<DosePoint> pts;
    for (int i = 0; i < 12; ++i) {
        const double conc = 1e-10 * std::pow(10.0, i * 0.5);
        pts.push_back({conc, fourPl(conc, a, d, c, b), -1.0});
    }
    FitOptions opt;
    opt.profileLikelihood = true;
    const FitResult f = fitFourParameterLogistic(pts, opt);
    REQUIRE(f.converged);

    const InversePrediction ip = inversePredict(f, 0.5 * (a + d), pts);
    CHECK_THAT(ip.concentration, WithinRel(c, 1e-8));
    CHECK_FALSE(ip.extrapolated);
    CHECK(ip.lower <= ip.concentration);
    CHECK(ip.upper >= ip.concentration);

    // A response outside the fitted asymptotes has no preimage; saying so beats
    // returning a number.
    CHECK_FALSE(inversePredict(f, 500.0, pts).intervalFound);
}

TEST_CASE("well adapters keep exclusions out of the fit", "[assay][fits]") {
    std::vector<Well> wells;
    for (int i = 0; i < 12; ++i) {
        Well w;
        w.concentration = 1e-10 * std::pow(10.0, i * 0.5);
        w.readout = fourPl(w.concentration, 100.0, 0.0, 1e-7, 1.0);
        w.role = WellRole::Sample;
        wells.push_back(w);
    }
    Well bad = wells.back();
    bad.readout = 1e6;
    bad.excluded = true;
    bad.exclusionRule = "grubbs";
    wells.push_back(bad);

    CHECK(doseSeriesFromWells(wells).size() == 12);
    const FitResult f = fitSeries(wells, AssayModel::FourParameterLogistic);
    REQUIRE(f.converged);
    CHECK_THAT(f.derivedEc50.value, WithinRel(1e-7, 1e-9));

    // Lineweaver-Burk exists for plotting only; it is not reachable as a fit.
    CHECK(lineweaverBurk(kineticSeriesFromWells(wells)).size() == 12);
}
