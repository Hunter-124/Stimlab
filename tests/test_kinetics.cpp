// tests/test_kinetics.cpp - Arrhenius, Eyring, the pH-rate profile and the shelf-life
// extrapolation that replaced BioCAD's functional-group shelf-life bucket.
//
// The load-bearing assertion is the REFUSAL: two temperatures fit a straight line
// exactly and say nothing about how wrong it is, so no extrapolation and no shelf life
// are offered from them. That is the whole reason the old predictThermalWindow /
// predictPhWindow / shelf-life-string trio was deleted rather than tuned.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

#include "sim/Kinetics.h"

using namespace biocad;
using namespace biocad::sim;

namespace {
constexpr double kATrue = 4.7e11;    // 1/day
constexpr double kEaTrue = 92.5;     // kJ/mol

std::vector<double> temperatures() { return {293.15, 303.15, 313.15, 323.15, 333.15}; }

std::vector<double> arrheniusRates() {
    std::vector<double> k;
    for (double t : temperatures())
        k.push_back(kATrue * std::exp(-kEaTrue / (kGasConstantKjPerMolK * t)));
    return k;
}
}  // namespace

TEST_CASE("Arrhenius recovers A and Ea from noiseless data", "[kinetics]") {
    const KineticsFit fit = arrhenius(temperatures(), arrheniusRates());
    REQUIRE(std::abs(fit.preExponential.value - kATrue) / kATrue < 1e-9);
    REQUIRE(std::abs(fit.activationEnergy.value - kEaTrue) / kEaTrue < 1e-9);
    REQUIRE(fit.activationEnergy.unit == "kJ/mol");
    REQUIRE(fit.activationEnergy.provenance == Provenance::Predicted);
    REQUIRE(fit.rSquared > 1.0 - 1e-12);
    REQUIRE(fit.extrapolationSupported);
    REQUIRE_FALSE(fit.confidenceEllipse.empty());

    const double expected = kATrue * std::exp(-kEaTrue / (kGasConstantKjPerMolK * 298.15));
    REQUIRE(std::abs(fit.predictedRateAt25C.value - expected) / expected < 1e-9);
    REQUIRE(fit.predictedRateAt25C.provenance == Provenance::Predicted);
}

TEST_CASE("two temperatures refuse to extrapolate", "[kinetics]") {
    const std::vector<double> k = arrheniusRates();
    const KineticsFit two = arrhenius({303.15, 323.15}, {k[1], k[3]});
    // The LINE is exact - Ea comes back right - but its uncertainty is unknowable.
    REQUIRE(std::abs(two.activationEnergy.value - kEaTrue) / kEaTrue < 1e-9);
    REQUIRE_FALSE(two.extrapolationSupported);
    REQUIRE(two.predictedRateAt25C.provenance == Provenance::NotComputed);
    REQUIRE(two.predictionIntervalLow == 0.0);
    REQUIRE(two.predictionIntervalHigh == 0.0);
    REQUIRE_FALSE(two.warnings.empty());
    // And no shelf life is produced from it.
    const Quantity none = shelfLife(two, 298.15, 0.10);
    REQUIRE(none.provenance == Provenance::NotComputed);
    REQUIRE(none.source.find("three or more") != std::string::npos);
}

TEST_CASE("shelf life is a first-order extrapolation of the fitted rate", "[kinetics]") {
    const KineticsFit fit = arrhenius(temperatures(), arrheniusRates());
    const Quantity sl = shelfLife(fit, 298.15, 0.10);
    REQUIRE(sl.provenance == Provenance::Predicted);
    const double kAt25 = kATrue * std::exp(-kEaTrue / (kGasConstantKjPerMolK * 298.15));
    const double expected = -std::log(0.9) / kAt25;
    REQUIRE(std::abs(sl.value - expected) / expected < 1e-9);
    // A fraction outside (0, 1) is refused rather than producing a signed time.
    REQUIRE(shelfLife(fit, 298.15, 0.0).provenance == Provenance::NotComputed);
    REQUIRE(shelfLife(fit, 298.15, 1.0).provenance == Provenance::NotComputed);
}

TEST_CASE("Eyring recovers dH* and dS* and reports their joint region", "[kinetics]") {
    // Generated from the Eyring form, NOT from the Arrhenius data: the extra factor of
    // T makes the two genuinely different functions, and fitting one to the other's
    // data recovers only an approximation.
    const double dh = 89.4, ds = -25.7;
    std::vector<double> k;
    for (double t : temperatures())
        k.push_back(kBoltzmannJPerK * t / kPlanckJs * std::exp(ds / 8.314462618) *
                    std::exp(-dh / (kGasConstantKjPerMolK * t)));
    const KineticsFit fit = eyring(temperatures(), k);
    REQUIRE(std::abs(fit.enthalpyOfActivation.value - dh) / dh < 1e-9);
    REQUIRE(std::abs(fit.entropyOfActivation.value - ds) / std::abs(ds) < 1e-9);
    REQUIRE(fit.enthalpyOfActivation.unit == "kJ/mol");
    REQUIRE(fit.entropyOfActivation.unit == "J/(mol*K)");
    REQUIRE(fit.confidenceEllipse.size() > 2);
    const double expected = kBoltzmannJPerK * 298.15 / kPlanckJs * std::exp(ds / 8.314462618) *
                            std::exp(-dh / (kGasConstantKjPerMolK * 298.15));
    REQUIRE(std::abs(fit.predictedRateAt25C.value / expected - 1.0) < 1e-9);
}

TEST_CASE("real scatter produces a real interval and a real ellipse", "[kinetics]") {
    const std::vector<double> clean = arrheniusRates();
    const double bump[5] = {1.06, 0.95, 1.03, 0.97, 1.02};
    std::vector<double> noisy;
    for (std::size_t i = 0; i < clean.size(); ++i) noisy.push_back(clean[i] * bump[i]);
    const KineticsFit fit = arrhenius(temperatures(), noisy);
    REQUIRE(fit.activationEnergy.error > 0.0);
    REQUIRE(fit.predictionIntervalHigh > fit.predictionIntervalLow);
    double minH = 1e300, maxH = -1e300, minS = 1e300, maxS = -1e300;
    for (const auto& [h, s] : fit.confidenceEllipse) {
        minH = std::min(minH, h);
        maxH = std::max(maxH, h);
        minS = std::min(minS, s);
        maxS = std::max(maxS, s);
    }
    REQUIRE(maxH - minH > 1e-6);
    REQUIRE(maxS - minS > 1e-6);
}

TEST_CASE("the pH-rate minimum comes from its closed form", "[kinetics]") {
    const double kh = 3.0e-2, k0 = 4.0e-6, koh = 8.0e1, pKw = 14.0;
    const double kw = std::pow(10.0, -pKw);
    std::vector<double> phs, kobs;
    for (int p = 1; p <= 10; ++p) {
        const double h = std::pow(10.0, -p);
        phs.push_back(p);
        kobs.push_back(kh * h + k0 + koh * kw / h);
    }
    const PhRateProfile pr = phRate(phs, kobs, pKw);
    REQUIRE(std::abs(pr.kAcid.value - kh) / kh < 1e-8);
    REQUIRE(std::abs(pr.kBase.value - koh) / koh < 1e-8);
    REQUIRE(std::abs(pr.kNeutral.value - k0) / k0 < 1e-6);

    const double phMin = 0.5 * (pKw + std::log10(kh / koh));
    const double rateMin = k0 + 2.0 * std::sqrt(kh * koh * kw);
    REQUIRE(std::abs(pr.minimumPh.value - phMin) < 1e-9);
    REQUIRE(std::abs(pr.minimumRate.value - rateMin) / rateMin < 1e-8);

    // The closed form really is the minimum: a dense scan of the fitted model agrees.
    double best = 1e300, bestPh = 0;
    for (int i = 0; i <= 200000; ++i) {
        const double p = 0.5 + 13.0 * i / 200000.0;
        const double h = std::pow(10.0, -p);
        const double v = pr.kAcid.value * h + pr.kNeutral.value + pr.kBase.value * kw / h;
        if (v < best) {
            best = v;
            bestPh = p;
        }
    }
    REQUIRE(std::abs(bestPh - pr.minimumPh.value) < 1e-3);
    // pKw is an assumption and is recorded as one.
    bool mentionsPkw = false;
    for (const std::string& w : pr.warnings)
        if (w.find("pKw") != std::string::npos) mentionsPkw = true;
    REQUIRE(mentionsPkw);
}

TEST_CASE("fewer than three pH points is refused", "[kinetics]") {
    const PhRateProfile few = phRate({3.0, 7.0}, {1.0, 2.0});
    REQUIRE(few.kAcid.provenance == Provenance::NotComputed);
    REQUIRE(few.minimumPh.provenance == Provenance::NotComputed);
    REQUIRE_FALSE(few.warnings.empty());
}
