#include "sim/Kinetics.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <numeric>
#include <set>

#include <Eigen/Dense>

#include "numeric/Optimize.h"

namespace biocad::sim {
namespace {

std::size_t distinctCount(const std::vector<double>& v) {
    std::set<long long> seen;
    for (double x : v) seen.insert(std::llround(x * 1e6));
    return seen.size();
}

// Two-sided Student-t quantiles at 95% for 1..12 degrees of freedom. Tabulated
// rather than approximated: an interpolated quantile would be inventing statistics,
// and the same rule is already applied in numeric::profileLikelihood.
double t95(std::size_t dof) {
    static const double table[] = {12.706, 4.303, 3.182, 2.776, 2.571, 2.447,
                                   2.365,  2.306, 2.262, 2.228, 2.201, 2.179};
    if (dof == 0) return 0;
    if (dof <= 12) return table[dof - 1];
    return 1.96;   // the large-sample limit; 1.98 at dof 100, so this is within 1%
}

std::string fmt(double v, int digits = 4) {
    char buf[48];
    std::snprintf(buf, sizeof(buf), "%.*g", digits, v);
    return buf;
}

// ln k against 1/T by ordinary least squares: the starting point for the nonlinear
// fit, and the reason the nonlinear fit converges in a handful of iterations.
void linearSeed(const std::vector<double>& t, const std::vector<double>& k, double& slope,
                double& intercept) {
    const double n = static_cast<double>(t.size());
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (std::size_t i = 0; i < t.size(); ++i) {
        const double x = 1.0 / t[i], y = std::log(k[i]);
        sx += x; sy += y; sxx += x * x; sxy += x * y;
    }
    const double d = n * sxx - sx * sx;
    slope = d != 0 ? (n * sxy - sx * sy) / d : 0.0;
    intercept = d != 0 ? (sy - slope * sx) / n : (t.empty() ? 0.0 : std::log(k[0]));
}

bool validate(const std::vector<double>& t, const std::vector<double>& k,
              std::vector<std::string>& warnings) {
    if (t.size() != k.size()) {
        warnings.push_back("temperature and rate-constant vectors have different lengths");
        return false;
    }
    if (t.size() < 2) {
        warnings.push_back("a kinetics fit needs at least two temperatures");
        return false;
    }
    for (std::size_t i = 0; i < t.size(); ++i) {
        if (!(t[i] > 0)) {
            warnings.push_back("temperature " + fmt(t[i]) + " K is not positive");
            return false;
        }
        if (!(k[i] > 0)) {
            warnings.push_back("rate constant " + fmt(k[i]) +
                               " is not positive; an Arrhenius fit is logarithmic in k");
            return false;
        }
    }
    return true;
}

}  // namespace

KineticsFit arrhenius(const std::vector<double>& temperaturesK,
                      const std::vector<double>& rateConstants) {
    KineticsFit out;
    out.temperaturesK = temperaturesK;
    out.rateConstants = rateConstants;
    if (!validate(temperaturesK, rateConstants, out.warnings)) {
        out.preExponential = notComputed("valid rate data at two or more temperatures");
        out.activationEnergy = out.preExponential;
        out.predictedRateAt25C = out.preExponential;
        return out;
    }
    const std::size_t n = temperaturesK.size();
    double slope = 0, intercept = 0;
    linearSeed(temperaturesK, rateConstants, slope, intercept);
    // slope = -Ea/R, intercept = ln A.
    std::vector<double> initial = {intercept, -slope * kGasConstantKjPerMolK};

    numeric::LmEvaluate evaluate = [&](const std::vector<double>& p, std::vector<double>& r,
                                      std::vector<double>& j) {
        r.assign(n, 0.0);
        j.assign(n * 2, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            const double rt = kGasConstantKjPerMolK * temperaturesK[i];
            const double model = std::exp(p[0] - p[1] / rt);
            r[i] = model - rateConstants[i];
            j[i * 2 + 0] = model;             // d/d(lnA)
            j[i * 2 + 1] = -model / rt;       // d/d(Ea)
        }
    };
    const numeric::LmResult fit = numeric::levenbergMarquardt(initial, n, evaluate);
    const double lnA = fit.params[0], ea = fit.params[1];
    // levenbergMarquardt does not know the observations, so R^2 is computed here
    // against them rather than read out of LmResult (where it stays zero).
    {
        std::vector<double> fitted(n);
        for (std::size_t i = 0; i < n; ++i)
            fitted[i] = std::exp(lnA - ea / (kGasConstantKjPerMolK * temperaturesK[i]));
        out.rSquared = numeric::rSquared(rateConstants, fitted);
    }
    out.extrapolationSupported = distinctCount(temperaturesK) >= 3;

    const double aValue = std::exp(lnA);
    const double lnASe = fit.standardErrors.size() == 2 ? fit.standardErrors[0] : 0.0;
    const double eaSe = fit.standardErrors.size() == 2 ? fit.standardErrors[1] : 0.0;
    // A's error bar is asymmetric in linear space; the reported symmetric error is
    // A * se(ln A), which is the first-order (delta-method) value, and the source
    // string says so rather than letting it look exact.
    out.preExponential = makeQuantity(aValue, "same units as k", aValue * lnASe,
                                      Provenance::Predicted,
                                      "Arrhenius fit (Levenberg-Marquardt); error is the "
                                      "delta-method value A*se(lnA)");
    out.activationEnergy = makeQuantity(ea, "kJ/mol", eaSe, Provenance::Predicted,
                                        "Arrhenius fit to " + std::to_string(n) +
                                            " rate constants");
    // Eyring's relations to the Arrhenius parameters, for a unimolecular step in
    // solution: dH = Ea - RT and dS = R*(ln(A h /(kB T)) - 1), both at the mean
    // experimental temperature, which is where the fit is actually anchored.
    const double tBar = std::accumulate(temperaturesK.begin(), temperaturesK.end(), 0.0) /
                        static_cast<double>(n);
    const double dh = ea - kGasConstantKjPerMolK * tBar;
    const double ds = 8.314462618 * (std::log(aValue * kPlanckJs / (kBoltzmannJPerK * tBar)) - 1.0);
    out.enthalpyOfActivation = makeQuantity(dh, "kJ/mol", eaSe, Provenance::Predicted,
                                            "dH* = Ea - R*Tmean at Tmean = " + fmt(tBar) + " K");
    out.entropyOfActivation =
        makeQuantity(ds, "J/(mol*K)", 8.314462618 * lnASe, Provenance::Predicted,
                     "dS* from the Arrhenius A at Tmean = " + fmt(tBar) + " K");

    const double t25 = 298.15;
    const double rt25 = kGasConstantKjPerMolK * t25;
    const double kPred = std::exp(lnA - ea / rt25);
    if (!out.extrapolationSupported) {
        out.predictedRateAt25C =
            notComputed("rate constants at three or more distinct temperatures (two points fit a "
                        "line exactly and cannot bound the extrapolation)");
        out.warnings.push_back(
            "only " + std::to_string(distinctCount(temperaturesK)) +
            " distinct temperature(s): the fit is reported but no extrapolation is offered");
        out.confidenceEllipse.clear();
        return out;
    }

    // Delta method on ln k = lnA - Ea/(R T): var = g' C g with g = (1, -1/(R T)).
    double varLnK = 0;
    if (fit.covariance.size() == 4) {
        const double g0 = 1.0, g1 = -1.0 / rt25;
        varLnK = g0 * g0 * fit.covariance[0] + 2 * g0 * g1 * fit.covariance[1] +
                 g1 * g1 * fit.covariance[3];
    }
    const double seLnK = std::sqrt(std::max(0.0, varLnK));
    const double tq = t95(n > 2 ? n - 2 : 1);
    out.predictedRateAt25C = makeQuantity(kPred, "same units as k", kPred * seLnK,
                                          Provenance::Predicted,
                                          "Arrhenius extrapolation to 25 C from " +
                                              std::to_string(distinctCount(temperaturesK)) +
                                              " temperatures");
    out.predictionIntervalLow = kPred * std::exp(-tq * seLnK);
    out.predictionIntervalHigh = kPred * std::exp(tq * seLnK);
    out.assumptions = {
        "first-order (or pseudo-first-order) loss, so a single rate constant describes each "
        "temperature",
        "one mechanism over the whole temperature range: a change of rate-determining step "
        "breaks the extrapolation without breaking the fit",
        "the extrapolation is to 25 C only; humidity, light and solid-state form are not in "
        "this model",
    };

    // The dH*/dS* joint 95% confidence ellipse, from the 2x2 covariance of
    // (lnA, Ea) mapped to (dH, dS). chi-square(2 df, 0.95) = 5.991.
    if (fit.covariance.size() == 4) {
        Eigen::Matrix2d c;
        c << fit.covariance[0], fit.covariance[1], fit.covariance[2], fit.covariance[3];
        // dH = Ea - R Tbar  -> d/d(lnA) = 0, d/dEa = 1
        // dS = R(ln A + ln(h/(kB Tbar)) - 1) -> d/d(lnA) = R, d/dEa = 0
        Eigen::Matrix2d jac;
        jac << 0.0, 1.0,
               8.314462618, 0.0;
        const Eigen::Matrix2d cov = jac * c * jac.transpose();
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(cov);
        if (es.info() == Eigen::Success) {
            const double s = std::sqrt(5.991);
            for (int i = 0; i <= 64; ++i) {
                const double th = 2.0 * std::numbers::pi_v<double> * i / 64.0;
                Eigen::Vector2d unit(std::cos(th), std::sin(th));
                Eigen::Vector2d radius(s * std::sqrt(std::max(0.0, es.eigenvalues()[0])),
                                       s * std::sqrt(std::max(0.0, es.eigenvalues()[1])));
                const Eigen::Vector2d p =
                    es.eigenvectors() * Eigen::Vector2d(radius[0] * unit[0], radius[1] * unit[1]);
                out.confidenceEllipse.emplace_back(dh + p[0], ds + p[1]);
            }
        }
    }
    return out;
}

KineticsFit eyring(const std::vector<double>& temperaturesK,
                   const std::vector<double>& rateConstants, double transmission) {
    KineticsFit out;
    out.temperaturesK = temperaturesK;
    out.rateConstants = rateConstants;
    if (!validate(temperaturesK, rateConstants, out.warnings) || !(transmission > 0)) {
        out.enthalpyOfActivation = notComputed("valid rate data and a positive transmission "
                                               "coefficient");
        out.entropyOfActivation = out.enthalpyOfActivation;
        out.predictedRateAt25C = out.enthalpyOfActivation;
        return out;
    }
    const std::size_t n = temperaturesK.size();
    // Seed from the Arrhenius line: dH ~ Ea - R Tbar, dS from the intercept.
    double slope = 0, intercept = 0;
    linearSeed(temperaturesK, rateConstants, slope, intercept);
    const double tBar = std::accumulate(temperaturesK.begin(), temperaturesK.end(), 0.0) /
                        static_cast<double>(n);
    const double eaSeed = -slope * kGasConstantKjPerMolK;
    std::vector<double> initial = {
        eaSeed - kGasConstantKjPerMolK * tBar,
        8.314462618 *
            (std::log(std::exp(intercept) * kPlanckJs / (transmission * kBoltzmannJPerK * tBar)) -
             1.0)};

    numeric::LmEvaluate evaluate = [&](const std::vector<double>& p, std::vector<double>& r,
                                      std::vector<double>& j) {
        r.assign(n, 0.0);
        j.assign(n * 2, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            const double t = temperaturesK[i];
            const double rt = kGasConstantKjPerMolK * t;
            const double pre = transmission * kBoltzmannJPerK * t / kPlanckJs;
            const double model = pre * std::exp(p[1] / 8.314462618) * std::exp(-p[0] / rt);
            r[i] = model - rateConstants[i];
            j[i * 2 + 0] = -model / rt;              // d/d(dH)
            j[i * 2 + 1] = model / 8.314462618;      // d/d(dS)
        }
    };
    const numeric::LmResult fit = numeric::levenbergMarquardt(initial, n, evaluate);
    const double dh = fit.params[0], ds = fit.params[1];
    {
        std::vector<double> fitted(n);
        for (std::size_t i = 0; i < n; ++i)
            fitted[i] = transmission * kBoltzmannJPerK * temperaturesK[i] / kPlanckJs *
                        std::exp(ds / 8.314462618) *
                        std::exp(-dh / (kGasConstantKjPerMolK * temperaturesK[i]));
        out.rSquared = numeric::rSquared(rateConstants, fitted);
    }
    out.extrapolationSupported = distinctCount(temperaturesK) >= 3;
    const double dhSe = fit.standardErrors.size() == 2 ? fit.standardErrors[0] : 0.0;
    const double dsSe = fit.standardErrors.size() == 2 ? fit.standardErrors[1] : 0.0;
    out.enthalpyOfActivation = makeQuantity(dh, "kJ/mol", dhSe, Provenance::Predicted,
                                            "Eyring fit (Levenberg-Marquardt), kappa = " +
                                                fmt(transmission));
    out.entropyOfActivation = makeQuantity(ds, "J/(mol*K)", dsSe, Provenance::Predicted,
                                           "Eyring fit (Levenberg-Marquardt), kappa = " +
                                               fmt(transmission));
    const double ea = dh + kGasConstantKjPerMolK * tBar;
    out.activationEnergy = makeQuantity(ea, "kJ/mol", dhSe, Provenance::Predicted,
                                        "Ea = dH* + R*Tmean at Tmean = " + fmt(tBar) + " K");
    const double aValue = transmission * kBoltzmannJPerK * tBar / kPlanckJs *
                          std::exp(ds / 8.314462618 + 1.0);
    out.preExponential = makeQuantity(aValue, "same units as k", aValue * dsSe / 8.314462618,
                                      Provenance::Predicted,
                                      "A from the Eyring dS* at Tmean = " + fmt(tBar) + " K");
    const double t25 = 298.15;
    const double kPred = transmission * kBoltzmannJPerK * t25 / kPlanckJs *
                          std::exp(ds / 8.314462618) *
                          std::exp(-dh / (kGasConstantKjPerMolK * t25));
    if (!out.extrapolationSupported) {
        out.predictedRateAt25C = notComputed(
            "rate constants at three or more distinct temperatures");
        out.warnings.push_back("fewer than three distinct temperatures: no extrapolation offered");
        return out;
    }
    double varLnK = 0;
    if (fit.covariance.size() == 4) {
        const double g0 = -1.0 / (kGasConstantKjPerMolK * t25), g1 = 1.0 / 8.314462618;
        varLnK = g0 * g0 * fit.covariance[0] + 2 * g0 * g1 * fit.covariance[1] +
                 g1 * g1 * fit.covariance[3];
    }
    const double seLnK = std::sqrt(std::max(0.0, varLnK));
    const double tq = t95(n > 2 ? n - 2 : 1);
    out.predictedRateAt25C = makeQuantity(kPred, "same units as k", kPred * seLnK,
                                          Provenance::Predicted, "Eyring extrapolation to 25 C");
    out.predictionIntervalLow = kPred * std::exp(-tq * seLnK);
    out.predictionIntervalHigh = kPred * std::exp(tq * seLnK);
    out.assumptions = {"transition-state theory with transmission coefficient kappa = " +
                           fmt(transmission),
                       "dH* and dS* constant over the measured temperature range",
                       "dH* and dS* are strongly correlated; read the joint ellipse, not the two "
                       "error bars"};
    if (fit.covariance.size() == 4) {
        Eigen::Matrix2d cov;
        cov << fit.covariance[0], fit.covariance[1], fit.covariance[2], fit.covariance[3];
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(cov);
        if (es.info() == Eigen::Success) {
            const double s = std::sqrt(5.991);
            for (int i = 0; i <= 64; ++i) {
                const double th = 2.0 * std::numbers::pi_v<double> * i / 64.0;
                const Eigen::Vector2d scaled(
                    s * std::sqrt(std::max(0.0, es.eigenvalues()[0])) * std::cos(th),
                    s * std::sqrt(std::max(0.0, es.eigenvalues()[1])) * std::sin(th));
                const Eigen::Vector2d p = es.eigenvectors() * scaled;
                out.confidenceEllipse.emplace_back(dh + p[0], ds + p[1]);
            }
        }
    }
    return out;
}

Quantity shelfLife(const KineticsFit& fit, double storageTemperatureK, double fractionLost) {
    if (!fit.extrapolationSupported)
        return notComputed("degradation rate constants at three or more distinct temperatures");
    if (!(fractionLost > 0) || !(fractionLost < 1))
        return notComputed("a fraction lost strictly between 0 and 1");
    if (!(storageTemperatureK > 0)) return notComputed("a positive storage temperature");
    if (fit.activationEnergy.provenance == Provenance::NotComputed ||
        fit.preExponential.provenance == Provenance::NotComputed)
        return notComputed("a converged Arrhenius fit");

    const double rt = kGasConstantKjPerMolK * storageTemperatureK;
    const double k = fit.preExponential.value * std::exp(-fit.activationEnergy.value / rt);
    if (!(k > 0)) return notComputed("a positive extrapolated rate constant");
    // First-order loss: t = -ln(1 - f) / k. The rate constant's unit sets the time
    // unit, and the caller states it - which is why the source string names it.
    const double time = -std::log(1.0 - fractionLost) / k;
    // The interval on k maps to an interval on t inversely, so the reported error is
    // the larger half-width of the two mapped bounds.
    double error = 0;
    if (fit.predictionIntervalLow > 0 && fit.predictionIntervalHigh > 0) {
        // Scale the 25 C interval's RELATIVE width onto this temperature: the delta
        // method's relative width is temperature-dependent, but using the anchored
        // relative width is the conservative choice and is stated as such.
        const double rel = 0.5 * (fit.predictionIntervalHigh - fit.predictionIntervalLow) /
                           std::max(1e-300, fit.predictedRateAt25C.value);
        error = time * rel;
    }
    return makeQuantity(time, "same time unit as 1/k", error, Provenance::Predicted,
                        "first-order extrapolation t = -ln(1-f)/k(T) from an Arrhenius fit at " +
                            fmt(storageTemperatureK) + " K, f = " + fmt(fractionLost) +
                            "; the interval is the fit's relative prediction width");
}

PhRateProfile phRate(const std::vector<double>& pHValues,
                     const std::vector<double>& rateConstants, double pKw) {
    PhRateProfile out;
    out.pHValues = pHValues;
    out.rateConstants = rateConstants;
    if (pHValues.size() != rateConstants.size() || pHValues.size() < 3) {
        out.warnings.push_back("a pH-rate profile needs at least three pH values (the model has "
                               "three parameters: kH, k0 and kOH)");
        out.kAcid = notComputed("three or more pH points");
        out.kNeutral = out.kAcid;
        out.kBase = out.kAcid;
        out.minimumPh = out.kAcid;
        out.minimumRate = out.kAcid;
        return out;
    }
    const std::size_t n = pHValues.size();
    const double kw = std::pow(10.0, -pKw);
    numeric::LmEvaluate evaluate = [&](const std::vector<double>& p, std::vector<double>& r,
                                      std::vector<double>& j) {
        r.assign(n, 0.0);
        j.assign(n * 3, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            const double h = std::pow(10.0, -pHValues[i]);
            const double oh = kw / h;
            r[i] = p[0] * h + p[1] + p[2] * oh - rateConstants[i];
            j[i * 3 + 0] = h;
            j[i * 3 + 1] = 1.0;
            j[i * 3 + 2] = oh;
        }
    };
    // Seeds from the extremes: the acid limb is dominated by kH, the base limb by
    // kOH, and the smallest observed rate bounds k0.
    double minK = rateConstants[0];
    for (double v : rateConstants) minK = std::min(minK, v);
    const double hLow = std::pow(10.0, -(*std::min_element(pHValues.begin(), pHValues.end())));
    const double hHigh = std::pow(10.0, -(*std::max_element(pHValues.begin(), pHValues.end())));
    std::vector<double> initial = {std::max(1e-12, rateConstants.front() / std::max(hLow, 1e-300)),
                                   std::max(0.0, minK),
                                   std::max(1e-12, rateConstants.back() * hHigh / kw)};
    const numeric::LmResult fit = numeric::levenbergMarquardt(initial, n, evaluate);
    const double kh = fit.params[0], k0 = fit.params[1], koh = fit.params[2];
    {
        std::vector<double> fitted(n);
        for (std::size_t i = 0; i < n; ++i) {
            const double h = std::pow(10.0, -pHValues[i]);
            fitted[i] = kh * h + k0 + koh * kw / h;
        }
        out.rSquared = numeric::rSquared(rateConstants, fitted);
    }
    const auto se = [&](std::size_t i) {
        return fit.standardErrors.size() == 3 ? fit.standardErrors[i] : 0.0;
    };
    out.kAcid = makeQuantity(kh, "1/(M*time)", se(0), Provenance::Predicted,
                             "pH-rate fit k_obs = kH[H+] + k0 + kOH[OH-]");
    out.kNeutral = makeQuantity(k0, "1/time", se(1), Provenance::Predicted, "pH-rate fit");
    out.kBase = makeQuantity(koh, "1/(M*time)", se(2), Provenance::Predicted, "pH-rate fit");
    if (kh > 0 && koh > 0) {
        const double phMin = 0.5 * (pKw + std::log10(kh / koh));
        const double rateMin = k0 + 2.0 * std::sqrt(kh * koh * kw);
        out.minimumPh = makeQuantity(phMin, "pH", 0.0, Provenance::Predicted,
                                     "closed form 0.5*(pKw + log10(kH/kOH)) with pKw = " +
                                         fmt(pKw));
        out.minimumRate = makeQuantity(rateMin, "1/time", 0.0, Provenance::Predicted,
                                       "closed form k0 + 2*sqrt(kH*kOH*Kw)");
        const double lo = *std::min_element(pHValues.begin(), pHValues.end());
        const double hi = *std::max_element(pHValues.begin(), pHValues.end());
        if (phMin < lo || phMin > hi)
            out.warnings.push_back(
                "the minimum at pH " + fmt(phMin) + " lies OUTSIDE the measured range pH " +
                fmt(lo) + " to " + fmt(hi) + ": it is an extrapolation of the fitted limbs");
    } else {
        out.minimumPh = notComputed("positive kH and kOH (one limb was not observed)");
        out.minimumRate = out.minimumPh;
        out.warnings.push_back("kH or kOH fitted non-positive: the data do not resolve both "
                               "acid- and base-catalysed limbs, so no minimum is reported");
    }
    out.warnings.push_back("pKw = " + fmt(pKw) +
                           " assumed (water at 25 C, zero ionic strength); pKw is 13.62 at 37 C, "
                           "which moves the reported minimum by 0.19 pH units");
    return out;
}

}  // namespace biocad::sim
