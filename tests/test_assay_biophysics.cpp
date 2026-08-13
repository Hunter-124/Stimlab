// Tests for the SPR/BLI, DSF and ITC fits. Every truth value here is generated from a
// closed form or an independent fine-step integrator, so a passing case means the
// module recovered a known answer rather than reproducing its own arithmetic.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "assay/Biophysics.h"

using namespace biocad;
using namespace biocad::assay;
using Catch::Matchers::ContainsSubstring;

namespace {

const FittedParameter* find(const FitResult& r, const std::string& name) {
    for (const FittedParameter& p : r.parameters) {
        if (p.name == name) {
            return &p;
        }
    }
    return nullptr;
}

double value(const FitResult& r, const std::string& name) {
    const FittedParameter* p = find(r, name);
    REQUIRE(p != nullptr);
    return p->value.value;
}

double relative(double got, double want) { return std::abs(got - want) / std::abs(want); }

// Closed-form 1:1 Langmuir association then dissociation - independent of the
// module's RK4 integrator, which is the point of using it as the oracle.
double langmuirClosed(double t, double conc, double ka, double kd, double rmax, double tDiss) {
    const double lambda = ka * conc + kd;
    const double req = rmax * ka * conc / lambda;
    if (t <= tDiss) {
        return req * (1.0 - std::exp(-lambda * t));
    }
    return req * (1.0 - std::exp(-lambda * tDiss)) * std::exp(-kd * (t - tDiss));
}

// Fine-step (0.01 s) RK4 for the mass-transport-limited model.
std::vector<double> mtlTrace(const std::vector<double>& times, double conc, double ka, double kd,
                             double rmax, double kt, double tDiss) {
    std::vector<double> out;
    out.reserve(times.size());
    double t = 0.0;
    double R = 0.0;
    auto f = [&](double tt, double r) {
        const double c = tt < tDiss ? conc : 0.0;
        const double free = rmax - r;
        return kt * (ka * c * free - kd * r) / (kt + ka * free);
    };
    for (double target : times) {
        while (t < target - 1e-12) {
            const double h = std::min(0.01, target - t);
            const double k1 = f(t, R);
            const double k2 = f(t + h / 2, R + h * k1 / 2);
            const double k3 = f(t + h / 2, R + h * k2 / 2);
            const double k4 = f(t + h, R + h * k3);
            R += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
            t += h;
        }
        out.push_back(R);
    }
    return out;
}

KineticExperiment langmuirSeries(const std::vector<double>& concs, double ka, double kd,
                                 double rmax, double tDiss, double tEnd, double dt) {
    KineticExperiment ex;
    ex.seriesId = "spr";
    for (double c : concs) {
        KineticCurve curve;
        curve.concentrationM = c;
        curve.dissociationStartS = tDiss;
        for (double t = 0.0; t <= tEnd + 1e-9; t += dt) {
            curve.timeS.push_back(t);
            curve.responseRu.push_back(langmuirClosed(t, c, ka, kd, rmax, tDiss));
        }
        ex.curves.push_back(curve);
    }
    return ex;
}

double boltzmann(double T, double fn, double sn, double fu, double su, double tm, double a) {
    const double frac = 1.0 / (1.0 + std::exp((tm - T) / a));
    const double bn = fn + sn * T;
    return bn + frac * ((fu + su * T) - bn);
}

MeltCurve syntheticMelt() {
    MeltCurve melt;
    melt.seriesId = "dsf";
    for (double T = 25.0; T <= 85.0 + 1e-9; T += 0.25) {
        melt.temperatureC.push_back(T);
        melt.signal.push_back(boltzmann(T, 100.0, 0.5, 1000.0, 1.0, 52.0, 1.5));
    }
    return melt;
}

struct ItcTruth {
    double n = 1.2;
    double k = 1e6;
    double dh = -10.0;
    double v0 = 1.4e-3;
    double m0 = 20e-6;
    double x0 = 300e-6;
    double blank = 0.3;
};

ItcExperiment syntheticItc(const ItcTruth& truth, std::size_t injections, double volume) {
    ItcExperiment ex;
    ex.seriesId = "itc";
    ex.cellVolumeL = truth.v0;
    ex.macromoleculeM = truth.m0;
    ex.titrantM = truth.x0;
    ex.temperatureK = 298.15;
    double injected = 0.0;
    double prevQ = 0.0;
    for (std::size_t i = 0; i < injections; ++i) {
        injected += volume;
        const double dilution = std::exp(-injected / truth.v0);
        const double mt = truth.m0 * dilution;
        const double xr = truth.x0 * (1.0 - dilution) / mt;
        const double b = 1.0 + xr / truth.n + 1.0 / (truth.n * truth.k * mt);
        const double q = 0.5 * truth.n * mt * truth.dh * truth.v0
                         * (b - std::sqrt(b * b - 4.0 * xr / truth.n));
        const double w = volume / truth.v0;
        ItcInjection inj;
        inj.volumeL = volume;
        inj.heatUcal = (q - prevQ + w * 0.5 * (q + prevQ)) * 1e9 + truth.blank;
        prevQ = q;
        ex.injections.push_back(inj);
        ex.blankHeatUcal.push_back(truth.blank);
    }
    return ex;
}

}  // namespace

TEST_CASE("Global 1:1 Langmuir recovers ka, kd, Rmax and KD from five concentrations",
          "[assay][biophysics]") {
    const double ka = 1e5;
    const double kd = 1e-3;
    const double rmax = 100.0;
    KineticExperiment ex =
        langmuirSeries({1e-9, 2e-9, 5e-9, 1e-8, 2e-8}, ka, kd, rmax, 3000.0, 6000.0, 20.0);
    ex.theoreticalRmaxRu = 95.0;

    const FitResult r = fitLangmuirKinetics(ex, {});
    REQUIRE(r.converged);
    REQUIRE(relative(value(r, "ka"), ka) < 1e-6);
    REQUIRE(relative(value(r, "kd"), kd) < 1e-6);
    REQUIRE(relative(value(r, "Rmax"), rmax) < 1e-6);
    REQUIRE(relative(r.derivedKd.value, 1e-8) < 1e-9);
    REQUIRE(r.derivedKd.provenance == Provenance::Model);
    REQUIRE(r.derivedKd.unit == "M");

    // Theoretical Rmax is reported beside the fitted one, with the ratio in the note.
    REQUIRE(value(r, "Rmax (theoretical)") == 95.0);
    REQUIRE_THAT(r.note, ContainsSubstring("theoretical Rmax"));

    // The plateau was reached, so a steady-state KD is allowed - and it agrees with the
    // kinetic one to a few percent, which is the independent check worth having.
    const FittedParameter* steady = find(r, "KD (steady state)");
    REQUIRE(steady != nullptr);
    REQUIRE(steady->value.provenance == Provenance::Model);
    REQUIRE(relative(steady->value.value, 1e-8) < 0.05);
}

TEST_CASE("A mass-transport-limited trace fools the 1:1 model and not the two-compartment one",
          "[assay][biophysics]") {
    const double ka = 1e5;
    const double kd = 1e-3;
    const double rmax = 100.0;
    const double kt = 5e6;
    KineticExperiment ex;
    ex.seriesId = "spr-mtl";
    for (double c : {1e-9, 2e-9, 5e-9, 1e-8, 2e-8}) {
        KineticCurve curve;
        curve.concentrationM = c;
        curve.dissociationStartS = 3000.0;
        for (double t = 0.0; t <= 6000.0 + 1e-9; t += 20.0) {
            curve.timeS.push_back(t);
        }
        curve.responseRu = mtlTrace(curve.timeS, c, ka, kd, rmax, kt, 3000.0);
        ex.curves.push_back(curve);
    }

    const FitResult plain = fitLangmuirKinetics(ex, {});
    const FitResult mtl = fitMassTransportKinetics(ex, {});

    // This is the whole reason the second model exists: R2 above 0.99 while ka is out
    // by nearly threefold. A good R2 is not evidence that the mechanism is right.
    REQUIRE(plain.rSquared > 0.99);
    REQUIRE(relative(value(plain, "ka"), ka) > 0.5);
    REQUIRE(mtl.rSquared > plain.rSquared);
    REQUIRE(mtl.converged);
    REQUIRE(relative(value(mtl, "ka"), ka) < 1e-5);
    REQUIRE(relative(value(mtl, "kd"), kd) < 1e-5);
    REQUIRE(relative(value(mtl, "Rmax"), rmax) < 1e-5);
    REQUIRE(relative(value(mtl, "kt"), kt) < 1e-4);
}

TEST_CASE("An association phase short of equilibrium refuses the steady-state KD",
          "[assay][biophysics]") {
    KineticExperiment ex =
        langmuirSeries({1e-9, 2e-9, 5e-9, 1e-8, 2e-8}, 1e5, 1e-3, 100.0, 100.0, 400.0, 5.0);
    const FitResult r = fitLangmuirKinetics(ex, {});

    const FittedParameter* steady = find(r, "KD (steady state)");
    REQUIRE(steady != nullptr);
    REQUIRE(steady->value.provenance == Provenance::NotComputed);
    REQUIRE_THAT(steady->value.source, ContainsSubstring("of the fitted equilibrium response"));
    REQUIRE_THAT(steady->value.source, ContainsSubstring("90%"));
    // The kinetic KD is still valid: it comes from the shape, not from a plateau.
    REQUIRE(r.derivedKd.provenance == Provenance::Model);
}

TEST_CASE("Boltzmann and the Savitzky-Golay derivative agree on a clean melt",
          "[assay][biophysics]") {
    const MeltCurve melt = syntheticMelt();
    const FitResult r = fitBoltzmannMelt(melt, {});
    REQUIRE(r.converged);
    REQUIRE(relative(value(r, "Tm"), 52.0) < 1e-6);
    REQUIRE(relative(value(r, "transition width"), 1.5) < 1e-6);

    const Quantity dtm = derivativeTm(melt);
    REQUIRE(dtm.provenance == Provenance::Measured);
    REQUIRE(std::abs(dtm.value - 52.0) < 0.2);
    REQUIRE(r.warnings.empty());          // agreement is within the stated tolerance
    REQUIRE_THAT(r.note, ContainsSubstring("Savitzky-Golay"));
}

TEST_CASE("A SYPRO Orange trace is truncated at its maximum and still recovers Tm",
          "[assay][biophysics]") {
    MeltCurve sypro = syntheticMelt();
    sypro.syproOrange = true;
    const double peak = boltzmann(65.0, 100.0, 0.5, 1000.0, 1.0, 52.0, 1.5);
    for (std::size_t i = 0; i < sypro.temperatureC.size(); ++i) {
        if (sypro.temperatureC[i] > 65.0) {
            sypro.signal[i] = peak - 20.0 * (sypro.temperatureC[i] - 65.0);
        }
    }

    const FitResult r = fitBoltzmannMelt(sypro, {});
    REQUIRE(r.observations < sypro.temperatureC.size());
    REQUIRE(relative(value(r, "Tm"), 52.0) < 1e-6);
    REQUIRE_THAT(r.warnings.front(), ContainsSubstring("dye behaviour"));
}

TEST_CASE("The two-state thermodynamic model recovers Tm and the van't Hoff enthalpy",
          "[assay][biophysics]") {
    const double dHm = 120.0;
    const double dCp = 1.5;
    const double tmK = 52.0 + 273.15;
    const double gasConstant = 1.987204259e-3;
    MeltCurve melt;
    melt.seriesId = "nanodsf";
    for (double T = 25.0; T <= 85.0 + 1e-9; T += 0.25) {
        const double tk = T + 273.15;
        const double dG = dHm * (1.0 - tk / tmK) - dCp * ((tmK - tk) + tk * std::log(tk / tmK));
        const double fu = 1.0 / (1.0 + std::exp(dG / (gasConstant * tk)));
        melt.temperatureC.push_back(T);
        melt.signal.push_back((100.0 + 0.5 * T) + fu * ((1000.0 + T) - (100.0 + 0.5 * T)));
    }

    MeltFitOptions options;
    options.deltaCpKcalPerMolK = dCp;
    const FitResult r = fitTwoStateMelt(melt, options);
    REQUIRE(r.converged);
    REQUIRE(relative(value(r, "Tm"), 52.0) < 1e-6);
    REQUIRE(relative(value(r, "dHm"), dHm) < 1e-5);
    // dCp is echoed as an assumption, never fitted.
    REQUIRE_THAT(r.assumptions.front(), ContainsSubstring("dCp fixed"));
}

TEST_CASE("Wiseman recovers n, K and dH and reports c before the parameters",
          "[assay][biophysics]") {
    const ItcTruth truth;
    const FitResult r = fitWisemanIsotherm(syntheticItc(truth, 25, 10e-6));
    REQUIRE(r.converged);
    REQUIRE(r.parameters.front().name == "c (Wiseman, n*K*Mt)");
    REQUIRE(relative(r.parameters.front().value.value, truth.n * truth.k * truth.m0) < 1e-6);
    REQUIRE(relative(r.parameters.front().value.value, 24.0) < 1e-6);   // hand calculation
    REQUIRE(relative(value(r, "n"), 1.2) < 1e-6);
    REQUIRE(relative(value(r, "K"), 1e6) < 1e-6);
    REQUIRE(relative(value(r, "dH"), -10.0) < 1e-6);
    // dG = -RT ln K and -T dS = dG - dH, both from the same fit.
    const double dG = -1.987204259e-3 * 298.15 * std::log(1e6);
    REQUIRE(relative(value(r, "dG"), dG) < 1e-6);
    REQUIRE(relative(value(r, "-T dS"), dG - (-10.0)) < 1e-6);
    REQUIRE(relative(r.derivedKd.value, 1e-6) < 1e-6);
    REQUIRE(r.warnings.empty());
    REQUIRE_THAT(r.note, ContainsSubstring("Wiseman c"));
}

TEST_CASE("Wiseman warns outside c ~ 1-1000 and refuses without a blank",
          "[assay][biophysics]") {
    ItcTruth hot;
    hot.k = 1e9;
    hot.m0 = 41.6666666666667e-6;   // n*K*Mt = 5e4
    REQUIRE(relative(hot.n * hot.k * hot.m0, 5e4) < 1e-9);
    const FitResult high = fitWisemanIsotherm(syntheticItc(hot, 25, 10e-6));
    REQUIRE_FALSE(high.warnings.empty());
    REQUIRE_THAT(high.warnings.front(), ContainsSubstring("outside c ~ 1-1000"));

    ItcExperiment noBlank = syntheticItc(ItcTruth{}, 25, 10e-6);
    noBlank.blankHeatUcal.clear();
    const FitResult refused = fitWisemanIsotherm(noBlank);
    REQUIRE_FALSE(refused.converged);
    REQUIRE_THAT(refused.note, ContainsSubstring("blank heat-of-dilution"));
    REQUIRE(refused.derivedKd.provenance == Provenance::NotComputed);
}
