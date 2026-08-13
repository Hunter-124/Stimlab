#include "assay/Fits.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>

#include "numeric/Optimize.h"

namespace biocad::assay {
namespace {

constexpr double kLn10 = 2.302585092994045684;

std::string sci(double v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.6g", v);
    return buf;
}

FitResult failed(AssayModel model, const std::string& reason) {
    FitResult f;
    f.model = model;
    f.converged = false;
    f.derivedEc50 = notComputed(reason);
    f.derivedKd = notComputed(reason);
    f.note = reason;
    return f;
}

// Replicate weighting. 1/sd^2 is the statistically correct weight when the
// replicate SD is known; 1/y^2 is the constant-relative-error stand-in when it is
// not. A point whose weight is undefined gets weight 0 - excluded from the fit,
// still on the plot - and is counted so the note can say how many.
std::vector<double> buildWeights(const std::vector<double>& y, const std::vector<double>& sd,
                                 const FitOptions& options, std::size_t& skipped,
                                 std::string& label) {
    skipped = 0;
    label.clear();
    std::vector<double> w;
    if (options.sdWeighting) {
        label = "1/sd^2 replicate weighting";
        w.assign(y.size(), 0.0);
        for (std::size_t i = 0; i < y.size(); ++i) {
            if (i < sd.size() && sd[i] > 0.0) {
                w[i] = 1.0 / (sd[i] * sd[i]);
            } else {
                ++skipped;
            }
        }
    } else if (options.inverseSquareWeighting) {
        label = "1/y^2 weighting";
        w.assign(y.size(), 0.0);
        for (std::size_t i = 0; i < y.size(); ++i) {
            if (y[i] > 0.0) {
                w[i] = 1.0 / (y[i] * y[i]);
            } else {
                ++skipped;
            }
        }
    }
    return w;
}

// Two models that both reach double-precision residuals must tie, not be ranked
// by rounding error. The floor is (1e-12 * ||y||)^2: below that the residual is
// arithmetic noise, and AICc's parameter penalty is then the only thing left to
// separate the candidates - which is the correct answer.
double flooredSsr(double ssr, const std::vector<double>& y) {
    double scale = 0.0;
    for (double v : y) scale += v * v;
    const double floorSsr = 1e-24 * scale;
    return ssr > floorSsr ? ssr : floorSsr;
}

struct FitInputs {
    AssayModel                     model;
    std::vector<std::string>       names;
    std::vector<std::string>       units;
    std::vector<double>            y;         // observations, in fit order
    std::vector<double>            weights;
    std::string                    weightLabel;
    std::size_t                    weightSkipped = 0;
    std::vector<std::string>       assumptions;
};

// Runs the one fitter (optionally its robust wrapper) and packs everything the
// caller must be able to audit: rank, condition number, AICc, the residuals, and
// the reason standard errors were refused when they were.
FitResult runFit(const FitInputs& in, const numeric::LmEvaluate& evaluate,
                 const std::vector<double>& initial, const FitOptions& options) {
    const std::size_t n = in.y.size();
    const std::size_t k = initial.size();

    numeric::LmResult lm;
    bool robustUsed = false;
    std::string robustNote;
    if (options.robust) {
        const numeric::IrlsResult irls =
            numeric::tukeyBiweight(initial, n, evaluate, in.weights);
        lm = irls.fit;
        robustUsed = true;
        robustNote = "Tukey-biweight IRLS (tuning 4.685), scale " + sci(irls.scale);
        if (!irls.note.empty()) robustNote += "; " + irls.note;
    } else {
        lm = numeric::levenbergMarquardt(initial, n, evaluate, in.weights);
    }

    // Residuals at the returned parameters, from the same evaluator the fit used.
    std::vector<double> residuals;
    std::vector<double> jacobian;
    evaluate(lm.params, residuals, jacobian);
    std::vector<double> fitted(n);
    double ssr = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        fitted[i] = in.y[i] + residuals[i];
        const double w = in.weights.empty() ? 1.0 : in.weights[i];
        ssr += w * residuals[i] * residuals[i];
    }

    FitResult f;
    f.model = in.model;
    f.observations = n;
    f.converged = lm.converged;
    f.robust = robustUsed;
    f.rank = lm.rank;
    f.conditionNumber = lm.conditionNumber;
    f.rSquared = numeric::rSquared(in.y, fitted);
    f.aicc = numeric::aicc(flooredSsr(ssr, in.y), n, k);
    f.residuals = residuals;
    f.fittedY = fitted;
    f.assumptions = in.assumptions;

    const bool haveErrors = lm.standardErrors.size() == k;
    const std::string source = "least-squares fit, n=" + std::to_string(n) + " points";
    f.parameters.resize(k);
    for (std::size_t j = 0; j < k; ++j) {
        FittedParameter& p = f.parameters[j];
        p.name = j < in.names.size() ? in.names[j] : ("p" + std::to_string(j));
        const std::string unit = j < in.units.size() ? in.units[j] : std::string();
        // A fitted parameter is a constructed model artefact, not a measurement:
        // Provenance::Model, with the covariance standard error as its error bar.
        p.value = makeQuantity(lm.params[j], unit, haveErrors ? lm.standardErrors[j] : 0.0,
                               Provenance::Model, source);
    }

    f.note = lm.note;
    const auto append = [&f](const std::string& s) {
        if (s.empty()) return;
        if (!f.note.empty()) f.note += "; ";
        f.note += s;
    };
    append(robustNote);
    if (!in.weightLabel.empty()) {
        std::string wl = in.weightLabel;
        if (in.weightSkipped > 0) {
            wl += " (" + std::to_string(in.weightSkipped) +
                  " point(s) excluded: undefined weight)";
        }
        append(wl);
    }
    if (!haveErrors) {
        append("standard errors refused: rank " + std::to_string(lm.rank) + " of " +
               std::to_string(k) + " (rank-deficient Jacobian, condition number " +
               sci(lm.conditionNumber) + ")");
        f.warnings.push_back("parameter standard errors unavailable; the experiment did not "
                             "identify all parameters");
    }

    if (options.profileLikelihood) {
        for (std::size_t j = 0; j < k; ++j) {
            const numeric::ProfileInterval pi = numeric::profileLikelihood(
                lm.params, n, evaluate, j, options.confidence, in.weights);
            f.parameters[j].profileLower = pi.lower;
            f.parameters[j].profileUpper = pi.upper;
            f.parameters[j].profileComputed = pi.lowerFound && pi.upperFound;
            if (!f.parameters[j].profileComputed && !pi.note.empty()) {
                f.warnings.push_back("profile interval for " + f.parameters[j].name + ": " +
                                     pi.note);
            }
        }
    }
    return f;
}

void flagExtrapolated(FitResult& f, double derived, double lo, double hi) {
    if (!(derived > 0.0) || derived < lo || derived > hi) {
        f.extrapolated = true;
        f.warnings.push_back("the derived concentration " + sci(derived) +
                             " mol/L lies outside the tested range " + sci(lo) + " to " +
                             sci(hi) + " mol/L; it is an extrapolation, not a result");
    }
}

std::size_t distinctCount(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v.size();
}

}  // namespace

// ---------------------------------------------------------------------------
// 4PL
// ---------------------------------------------------------------------------

FitResult fitFourParameterLogistic(const std::vector<DosePoint>& points,
                                   const FitOptions& options) {
    if (points.size() < 4) {
        return failed(AssayModel::FourParameterLogistic,
                      "four-parameter logistic needs at least 4 points, got " +
                          std::to_string(points.size()));
    }
    for (const auto& p : points) {
        if (!(p.concentration > 0.0)) {
            return failed(AssayModel::FourParameterLogistic,
                          "a non-positive concentration cannot be log-transformed");
        }
    }

    const std::size_t n = points.size();
    std::vector<double> x(n), y(n), sd(n);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = std::log10(points[i].concentration);
        y[i] = points[i].response;
        sd[i] = points[i].sd;
    }

    FitInputs in;
    in.model = AssayModel::FourParameterLogistic;
    // Named by WHERE the asymptote sits on the concentration axis, not by where it
    // sits on the plot: in this parameterisation p[0] is always the response as
    // concentration goes to infinity. Calling it "top" reads inverted on a
    // descending inhibition curve, where the high-concentration plateau is the low
    // response.
    in.names = {"responseAtHighConc", "responseAtLowConc", "log10EC50", "hillSlope"};
    in.units = {"", "", "", ""};   // response scale arbitrary; log10 EC50 and slope unitless
    in.y = y;
    in.weights = buildWeights(y, sd, options, in.weightSkipped, in.weightLabel);
    in.assumptions = {
        "the response is monotonic in concentration over the tested range",
        "the slope is an EMPIRICAL slope; amplification and receptor reserve bend it",
        "response asymptotes are on the assay's arbitrary scale, so they carry no unit"};

    const double yMax = *std::max_element(y.begin(), y.end());
    const double yMin = *std::min_element(y.begin(), y.end());
    const double mid = 0.5 * (yMax + yMin);
    double logEc50Guess = x[0];
    double bestGap = std::abs(y[0] - mid);
    for (std::size_t i = 1; i < n; ++i) {
        const double gap = std::abs(y[i] - mid);
        if (gap < bestGap) {
            bestGap = gap;
            logEc50Guess = x[i];
        }
    }
    // Ascending or descending decides the sign of the initial slope; starting with
    // the wrong sign makes LM walk the asymptote flat instead of the transition.
    const bool ascending = y.back() >= y.front();
    std::vector<double> initial{ascending ? yMax : yMin, ascending ? yMin : yMax, logEc50Guess,
                                1.0};

    // u = 10^(B*(log10 C - x)); y = D + (A-D)/(1+u)
    //   dy/dA = 1/(1+u), dy/dD = u/(1+u)
    //   dy/dlog10C = -(A-D) * B * ln10 * u/(1+u)^2
    //   dy/dB      = -(A-D) * ln10 * (log10 C - x) * u/(1+u)^2
    auto evaluate = [&](const std::vector<double>& p, std::vector<double>& r,
                        std::vector<double>& jac) {
        const double a = p[0], d = p[1], logC = p[2], b = p[3];
        const double delta = a - d;
        r.resize(n);
        jac.resize(n * 4);
        for (std::size_t i = 0; i < n; ++i) {
            const double shift = logC - x[i];
            const double u = std::pow(10.0, b * shift);
            const double den = 1.0 + u;
            r[i] = d + delta / den - y[i];
            const double common = -delta * u / (den * den);
            jac[i * 4 + 0] = 1.0 / den;
            jac[i * 4 + 1] = u / den;
            jac[i * 4 + 2] = common * b * kLn10;
            jac[i * 4 + 3] = common * kLn10 * shift;
        }
    };

    FitResult f = runFit(in, evaluate, initial, options);
    f.fittedX = x;

    const double logEc50 = f.parameters[2].value.value;
    const double seLog = f.parameters[2].value.error;
    const double ec50 = std::pow(10.0, logEc50);
    f.derivedEc50 = makeQuantity(ec50, "mol/L", kLn10 * ec50 * seLog, Provenance::Model,
                                 "4PL fit; EC50 = 10^(fitted log10 EC50)");
    const double lo = std::pow(10.0, *std::min_element(x.begin(), x.end()));
    const double hi = std::pow(10.0, *std::max_element(x.begin(), x.end()));
    flagExtrapolated(f, ec50, lo, hi);
    return f;
}

// ---------------------------------------------------------------------------
// 5PL
// ---------------------------------------------------------------------------

FitResult fitFiveParameterLogistic(const std::vector<DosePoint>& points,
                                   const FitOptions& options) {
    for (const auto& p : points) {
        if (!(p.concentration > 0.0)) {
            return failed(AssayModel::FiveParameterLogistic,
                          "a non-positive concentration cannot be log-transformed");
        }
    }
    std::vector<double> concs;
    concs.reserve(points.size());
    for (const auto& p : points) concs.push_back(p.concentration);
    const std::size_t distinct = distinctCount(concs);
    if (distinct < 8) {
        return failed(AssayModel::FiveParameterLogistic,
                      "five-parameter logistic requires at least 8 distinct concentrations to "
                      "identify the asymmetry parameter G, got " + std::to_string(distinct));
    }

    const std::size_t n = points.size();
    std::vector<double> lx(n), y(n), sd(n);
    for (std::size_t i = 0; i < n; ++i) {
        lx[i] = std::log10(points[i].concentration);
        y[i] = points[i].response;
        sd[i] = points[i].sd;
    }

    FitInputs in;
    in.model = AssayModel::FiveParameterLogistic;
    // NOTE the index order is the OPPOSITE of the 4PL above: here p[0] = A is the
    // asymptote as concentration goes to ZERO ((x/C)^B -> 0), because that is the
    // published 5PL parameterisation. Position-explicit names are what stop that
    // from becoming a silent transcription error in a panel.
    in.names = {"responseAtLowConc", "responseAtHighConc", "log10C", "slopeB", "asymmetryG"};
    in.units = {"", "", "", "", ""};
    in.y = y;
    in.weights = buildWeights(y, sd, options, in.weightSkipped, in.weightLabel);
    in.assumptions = {
        "C is the curve's inflection scale, NOT the EC50 when G != 1",
        "the reported EC50 is C*(2^(1/G)-1)^(1/B), the true half-maximal concentration",
        "response asymptotes are on the assay's arbitrary scale, so they carry no unit"};

    const double yMax = *std::max_element(y.begin(), y.end());
    const double yMin = *std::min_element(y.begin(), y.end());
    const double mid = 0.5 * (yMax + yMin);
    double logCGuess = lx[0];
    double bestGap = std::abs(y[0] - mid);
    for (std::size_t i = 1; i < n; ++i) {
        const double gap = std::abs(y[i] - mid);
        if (gap < bestGap) {
            bestGap = gap;
            logCGuess = lx[i];
        }
    }
    // A is the low-concentration asymptote of this parameterisation ((x/C)^B -> 0),
    // D the high-concentration one.
    const bool ascending = y.back() >= y.front();
    std::vector<double> initial{ascending ? yMin : yMax, ascending ? yMax : yMin, logCGuess, 1.0,
                               1.0};

    // t = (x/C)^B = 10^(B*(log10 x - log10 C)); w = 1+t; F = w^-G; y = D + (A-D)F
    //   dy/dA = F, dy/dD = 1-F
    //   dy/dlog10C = (A-D) * G*B*ln10*t*w^(-G-1)
    //   dy/dB      = -(A-D) * G*t*ln(x/C)*w^(-G-1)
    //   dy/dG      = -(A-D) * F * ln(w)
    auto evaluate = [&](const std::vector<double>& p, std::vector<double>& r,
                        std::vector<double>& jac) {
        const double a = p[0], d = p[1], logC = p[2], b = p[3], g = p[4];
        const double delta = a - d;
        r.resize(n);
        jac.resize(n * 5);
        for (std::size_t i = 0; i < n; ++i) {
            const double lr = lx[i] - logC;             // log10(x/C)
            const double t = std::pow(10.0, b * lr);
            const double w = 1.0 + t;
            const double f = std::pow(w, -g);
            const double wg1 = std::pow(w, -g - 1.0);
            r[i] = d + delta * f - y[i];
            jac[i * 5 + 0] = f;
            jac[i * 5 + 1] = 1.0 - f;
            jac[i * 5 + 2] = delta * g * b * kLn10 * t * wg1;
            jac[i * 5 + 3] = -delta * g * t * (kLn10 * lr) * wg1;
            jac[i * 5 + 4] = -delta * f * std::log(w);
        }
    };

    FitResult f = runFit(in, evaluate, initial, options);
    f.fittedX = lx;

    const double c = std::pow(10.0, f.parameters[2].value.value);
    const double b = f.parameters[3].value.value;
    const double g = f.parameters[4].value.value;
    // The half-maximal point of an asymmetric logistic. Reporting C here instead
    // would be wrong by (2^(1/G)-1)^(1/B), which is an order of magnitude for
    // plausible G and B.
    const double ec50 = c * std::pow(std::pow(2.0, 1.0 / g) - 1.0, 1.0 / b);
    f.derivedEc50 = makeQuantity(ec50, "mol/L", 0.0, Provenance::Model,
                                 "5PL fit; EC50 = C*(2^(1/G)-1)^(1/B), not C");
    f.warnings.push_back("5PL EC50 has no propagated standard error: it is a nonlinear function "
                         "of three correlated parameters; use the profile interval of log10C");
    const double lo = std::pow(10.0, *std::min_element(lx.begin(), lx.end()));
    const double hi = std::pow(10.0, *std::max_element(lx.begin(), lx.end()));
    flagExtrapolated(f, ec50, lo, hi);
    return f;
}

// ---------------------------------------------------------------------------
// Enzyme kinetics
// ---------------------------------------------------------------------------

FitResult fitMichaelisMenten(const std::vector<KineticPoint>& points,
                             const FitOptions& options) {
    if (points.size() < 3) {
        return failed(AssayModel::MichaelisMenten,
                      "Michaelis-Menten needs at least 3 points, got " +
                          std::to_string(points.size()));
    }
    const std::size_t n = points.size();
    std::vector<double> s(n), y(n), sd(n);
    for (std::size_t i = 0; i < n; ++i) {
        s[i] = points[i].substrate;
        y[i] = points[i].velocity;
        sd[i] = points[i].sd;
    }

    FitInputs in;
    in.model = AssayModel::MichaelisMenten;
    in.names = {"Vmax", "Km"};
    in.units = {"", "mol/L"};
    in.y = y;
    in.weights = buildWeights(y, sd, options, in.weightSkipped, in.weightLabel);
    in.assumptions = {"initial rates, measured in the linear phase of product formation",
                      "[S] >> [E]t, so free substrate equals total substrate",
                      "fitted on the untransformed rates, never on a reciprocal plot"};

    const double vGuess = *std::max_element(y.begin(), y.end()) * 1.2;
    std::vector<double> sorted = s;
    std::sort(sorted.begin(), sorted.end());
    std::vector<double> initial{vGuess, sorted[sorted.size() / 2]};

    auto evaluate = [&](const std::vector<double>& p, std::vector<double>& r,
                        std::vector<double>& jac) {
        const double vmax = p[0], km = p[1];
        r.resize(n);
        jac.resize(n * 2);
        for (std::size_t i = 0; i < n; ++i) {
            const double den = km + s[i];
            r[i] = vmax * s[i] / den - y[i];
            jac[i * 2 + 0] = s[i] / den;
            jac[i * 2 + 1] = -vmax * s[i] / (den * den);
        }
    };

    FitResult f = runFit(in, evaluate, initial, options);
    f.fittedX = s;
    return f;
}

FitResult fitHill(const std::vector<KineticPoint>& points, const FitOptions& options) {
    if (points.size() < 4) {
        return failed(AssayModel::Hill,
                      "the Hill model needs at least 4 points, got " +
                          std::to_string(points.size()));
    }
    for (const auto& p : points) {
        if (!(p.substrate > 0.0)) {
            return failed(AssayModel::Hill, "a non-positive [S] cannot be raised to a power");
        }
    }
    const std::size_t n = points.size();
    std::vector<double> s(n), y(n), sd(n);
    for (std::size_t i = 0; i < n; ++i) {
        s[i] = points[i].substrate;
        y[i] = points[i].velocity;
        sd[i] = points[i].sd;
    }

    FitInputs in;
    in.model = AssayModel::Hill;
    in.names = {"Vmax", "log10K", "hillCoefficient"};
    in.units = {"", "", ""};
    in.y = y;
    in.weights = buildWeights(y, sd, options, in.weightSkipped, in.weightLabel);
    in.assumptions = {"the Hill coefficient is empirical and is not a subunit count",
                      "K is the half-saturating [S], reported as log10 for conditioning"};

    std::vector<double> sorted = s;
    std::sort(sorted.begin(), sorted.end());
    std::vector<double> initial{*std::max_element(y.begin(), y.end()) * 1.2,
                                std::log10(sorted[sorted.size() / 2]), 1.0};

    // a = S^h, b = K^h = 10^(h*log10K); v = Vmax*a/(a+b)
    //   dv/dVmax = a/(a+b)
    //   dv/dlog10K = -Vmax*a*b*h*ln10/(a+b)^2
    //   dv/dh      =  Vmax*a*b*(ln S - ln K)/(a+b)^2
    auto evaluate = [&](const std::vector<double>& p, std::vector<double>& r,
                        std::vector<double>& jac) {
        const double vmax = p[0], logK = p[1], h = p[2];
        const double lnK = logK * kLn10;
        r.resize(n);
        jac.resize(n * 3);
        for (std::size_t i = 0; i < n; ++i) {
            const double lnS = std::log(s[i]);
            const double a = std::exp(h * lnS);
            const double b = std::exp(h * lnK);
            const double den = a + b;
            r[i] = vmax * a / den - y[i];
            jac[i * 3 + 0] = a / den;
            jac[i * 3 + 1] = -vmax * a * b * h * kLn10 / (den * den);
            jac[i * 3 + 2] = vmax * a * b * (lnS - lnK) / (den * den);
        }
    };

    FitResult f = runFit(in, evaluate, initial, options);
    f.fittedX = s;
    return f;
}

FitResult fitSubstrateInhibition(const std::vector<KineticPoint>& points,
                                 const FitOptions& options) {
    if (points.size() < 4) {
        return failed(AssayModel::SubstrateInhibition,
                      "substrate inhibition has 3 parameters and needs at least 4 points, got " +
                          std::to_string(points.size()));
    }
    const std::size_t n = points.size();
    std::vector<double> s(n), y(n), sd(n);
    for (std::size_t i = 0; i < n; ++i) {
        s[i] = points[i].substrate;
        y[i] = points[i].velocity;
        sd[i] = points[i].sd;
    }

    FitInputs in;
    in.model = AssayModel::SubstrateInhibition;
    in.names = {"Vmax", "Km", "Ki"};
    in.units = {"", "mol/L", "mol/L"};
    in.y = y;
    in.weights = buildWeights(y, sd, options, in.weightSkipped, in.weightLabel);
    in.assumptions = {"one inhibitory substrate site; Vmax is the unattainable asymptote",
                      "the rate peaks at [S] = sqrt(Km*Ki), which is reported as derivedKd"};

    std::vector<double> sorted = s;
    std::sort(sorted.begin(), sorted.end());
    const double sMid = sorted[sorted.size() / 2];
    const double sHi = sorted.back();
    std::vector<double> initial{*std::max_element(y.begin(), y.end()) * 2.0, sMid, sHi};

    // v = Vmax*S/(Km + S + S^2/Ki)
    //   dv/dVmax = S/den
    //   dv/dKm   = -Vmax*S/den^2
    //   dv/dKi   =  Vmax*S^3/(Ki^2 * den^2)
    auto evaluate = [&](const std::vector<double>& p, std::vector<double>& r,
                        std::vector<double>& jac) {
        const double vmax = p[0], km = p[1], ki = p[2];
        r.resize(n);
        jac.resize(n * 3);
        for (std::size_t i = 0; i < n; ++i) {
            const double sq = s[i] * s[i];
            const double den = km + s[i] + sq / ki;
            r[i] = vmax * s[i] / den - y[i];
            jac[i * 3 + 0] = s[i] / den;
            jac[i * 3 + 1] = -vmax * s[i] / (den * den);
            jac[i * 3 + 2] = vmax * s[i] * sq / (ki * ki * den * den);
        }
    };

    FitResult f = runFit(in, evaluate, initial, options);
    f.fittedX = s;
    const double km = f.parameters[1].value.value;
    const double ki = f.parameters[2].value.value;
    if (km > 0.0 && ki > 0.0) {
        f.derivedKd = makeQuantity(std::sqrt(km * ki), "mol/L", 0.0, Provenance::Model,
                                   "substrate inhibition; [S] of maximal rate = sqrt(Km*Ki)");
    } else {
        f.derivedKd = notComputed("a positive Km and Ki (the fit returned a non-positive one)");
    }
    return f;
}

double morrisonFraction(double enzymeTotal, double inhibitor, double kiApp) {
    if (!(enzymeTotal > 0.0)) return classicInhibitionFraction(inhibitor, kiApp);
    const double q = enzymeTotal + inhibitor + kiApp;
    const double disc = std::max(0.0, q * q - 4.0 * enzymeTotal * inhibitor);
    return 1.0 - (q - std::sqrt(disc)) / (2.0 * enzymeTotal);
}

double classicInhibitionFraction(double inhibitor, double kiApp) {
    return kiApp / (kiApp + inhibitor);
}

FitResult fitMorrisonTightBinding(const std::vector<DosePoint>& points, double enzymeTotal,
                                  const FitOptions& options) {
    if (points.size() < 3) {
        return failed(AssayModel::MorrisonTightBinding,
                      "Morrison tight binding needs at least 3 points, got " +
                          std::to_string(points.size()));
    }
    if (!(enzymeTotal > 0.0)) {
        return failed(AssayModel::MorrisonTightBinding,
                      "[E]t (a measured total enzyme concentration; the quadratic solution is "
                      "defined by it and it is never fitted)");
    }

    const std::size_t n = points.size();
    std::vector<double> ic(n), y(n), sd(n);
    for (std::size_t i = 0; i < n; ++i) {
        ic[i] = points[i].concentration;
        y[i] = points[i].response;
        sd[i] = points[i].sd;
    }

    FitInputs in;
    in.model = AssayModel::MorrisonTightBinding;
    in.names = {"v0", "KiApp"};
    in.units = {"", "mol/L"};
    in.y = y;
    in.weights = buildWeights(y, sd, options, in.weightSkipped, in.weightLabel);
    in.assumptions = {
        "[E]t is a measured input, not a fitted parameter",
        "the quadratic (Morrison) solution is used; the classic 1/(1+I/Ki) form assumes "
        "[I] >> [E]t and fails when Ki approaches [E]t",
        "Ki_app = Ki*(1 + [S]/Km) for a competitive tight binder; converting needs [S] and Km"};

    std::vector<double> sortedI = ic;
    std::sort(sortedI.begin(), sortedI.end());
    std::vector<double> initial{*std::max_element(y.begin(), y.end()),
                                std::max(sortedI[sortedI.size() / 2], enzymeTotal)};

    // v = v0 * (1 - (Q - sqrt(Q^2 - 4 Et I)) / (2 Et)),  Q = Et + I + K
    //   dv/dv0 = the bracket
    //   dv/dK  = -v0 * (1 - Q/sqrt(Q^2-4 Et I)) / (2 Et)
    const double et = enzymeTotal;
    auto evaluate = [&](const std::vector<double>& p, std::vector<double>& r,
                        std::vector<double>& jac) {
        const double v0 = p[0], k = p[1];
        r.resize(n);
        jac.resize(n * 2);
        for (std::size_t i = 0; i < n; ++i) {
            const double q = et + ic[i] + k;
            const double disc = std::max(1e-300, q * q - 4.0 * et * ic[i]);
            const double root = std::sqrt(disc);
            const double frac = 1.0 - (q - root) / (2.0 * et);
            r[i] = v0 * frac - y[i];
            jac[i * 2 + 0] = frac;
            jac[i * 2 + 1] = -v0 * (1.0 - q / root) / (2.0 * et);
        }
    };

    FitResult f = runFit(in, evaluate, initial, options);
    f.fittedX = ic;
    const double kApp = f.parameters[1].value.value;
    f.derivedKd = makeQuantity(kApp, "mol/L", f.parameters[1].value.error, Provenance::Model,
                               "Morrison quadratic fit; apparent Ki at the assay's [S]");
    // The number a reader needs to judge whether tight binding mattered at all.
    if (kApp > 0.0) {
        f.warnings.push_back("Ki_app / [E]t = " + sci(kApp / et) +
                             "; below ~10 the classic approximation is invalid, which is why "
                             "the quadratic form is used");
    }
    return f;
}

// ---------------------------------------------------------------------------
// Global inhibition modality
// ---------------------------------------------------------------------------

namespace {

// Denominator of v = Vmax*S/den for each modality, with its derivative with
// respect to every fitted parameter. Km and the Ki values are fitted directly (not
// as logs) so the reported parameters are the ones a chemist quotes.
// Local label for the note strings. data/Domain.h has the JSON enum mapping but no
// C++ label function, and adding one there would collide with a sibling agent's
// file, so the string lives with its only consumer.
const char* modalityLabel(InhibitionModality m) {
    switch (m) {
        case InhibitionModality::Competitive: return "competitive";
        case InhibitionModality::Noncompetitive: return "noncompetitive";
        case InhibitionModality::Uncompetitive: return "uncompetitive";
        case InhibitionModality::Mixed: return "mixed";
        case InhibitionModality::RadioligandBinding: return "radioligand-binding";
        case InhibitionModality::Unknown: break;
    }
    return "unknown";
}

struct Denominator {
    double value = 0;
    double dKm = 0;
    double dKi1 = 0;   // Ki (or Kic for mixed)
    double dKi2 = 0;   // Kiu for mixed, unused otherwise
};

Denominator denominatorFor(InhibitionModality m, double s, double i, const std::vector<double>& p) {
    const double km = p[1];
    Denominator d;
    switch (m) {
        case InhibitionModality::Competitive: {
            const double ki = p[2];
            d.value = km * (1.0 + i / ki) + s;
            d.dKm = 1.0 + i / ki;
            d.dKi1 = -km * i / (ki * ki);
            break;
        }
        case InhibitionModality::Uncompetitive: {
            const double ki = p[2];
            d.value = km + s * (1.0 + i / ki);
            d.dKm = 1.0;
            d.dKi1 = -s * i / (ki * ki);
            break;
        }
        case InhibitionModality::Noncompetitive: {
            const double ki = p[2];
            d.value = (km + s) * (1.0 + i / ki);
            d.dKm = 1.0 + i / ki;
            d.dKi1 = -(km + s) * i / (ki * ki);
            break;
        }
        case InhibitionModality::Mixed: {
            const double kic = p[2], kiu = p[3];
            d.value = km * (1.0 + i / kic) + s * (1.0 + i / kiu);
            d.dKm = 1.0 + i / kic;
            d.dKi1 = -km * i / (kic * kic);
            d.dKi2 = -s * i / (kiu * kiu);
            break;
        }
        default:
            break;
    }
    return d;
}

}  // namespace

ModelComparison fitInhibitionModality(const std::vector<InhibitionPoint>& matrix,
                                      const FitOptions& options) {
    ModelComparison out;
    const std::size_t n = matrix.size();
    if (n < 8) {
        out.conclusion = "the [S] x [I] matrix needs at least 8 cells to fit four candidate "
                         "modalities, got " + std::to_string(n);
        return out;
    }
    std::vector<double> sVals, iVals;
    for (const auto& p : matrix) {
        sVals.push_back(p.substrate);
        iVals.push_back(p.inhibitor);
    }
    const std::size_t distinctS = distinctCount(sVals);
    const std::size_t distinctI = distinctCount(iVals);

    std::vector<double> s(n), ic(n), y(n), sd(n);
    for (std::size_t j = 0; j < n; ++j) {
        s[j] = matrix[j].substrate;
        ic[j] = matrix[j].inhibitor;
        y[j] = matrix[j].velocity;
        sd[j] = matrix[j].sd;
    }

    std::size_t skipped = 0;
    std::string weightLabel;
    const std::vector<double> weights = buildWeights(y, sd, options, skipped, weightLabel);

    std::vector<double> sortedS = s, sortedI = ic;
    std::sort(sortedS.begin(), sortedS.end());
    std::sort(sortedI.begin(), sortedI.end());
    double medianI = sortedI[sortedI.size() / 2];
    if (!(medianI > 0.0)) medianI = sortedI.back() > 0.0 ? sortedI.back() : 1.0;
    const double vGuess = *std::max_element(y.begin(), y.end()) * 1.5;
    const double kmGuess = sortedS[sortedS.size() / 2];

    const InhibitionModality kCandidates[4] = {
        InhibitionModality::Competitive, InhibitionModality::Uncompetitive,
        InhibitionModality::Noncompetitive, InhibitionModality::Mixed};

    for (InhibitionModality m : kCandidates) {
        const bool mixed = m == InhibitionModality::Mixed;
        const std::size_t k = mixed ? 4u : 3u;

        FitInputs in;
        in.model = AssayModel::MichaelisMenten;   // the substrate law is MM in all four
        in.names = mixed ? std::vector<std::string>{"Vmax", "Km", "Kic", "Kiu"}
                         : std::vector<std::string>{"Vmax", "Km", "Ki"};
        in.units = mixed ? std::vector<std::string>{"", "mol/L", "mol/L", "mol/L"}
                         : std::vector<std::string>{"", "mol/L", "mol/L"};
        in.y = y;
        in.weights = weights;
        in.weightLabel = weightLabel;
        in.weightSkipped = skipped;
        in.assumptions = {
            "one global fit over the whole [S] x [I] matrix; per-[I] IC50s cannot "
            "distinguish modality",
            "initial rates under rapid equilibrium",
            std::to_string(distinctS) + " distinct [S] x " + std::to_string(distinctI) +
                " distinct [I]"};
        if (distinctS < 2) {
            in.assumptions.push_back(
                "only one [S] was tested: competitive and uncompetitive coincide at [S] = Km, "
                "so the matrix cannot identify the modality");
        }

        auto evaluate = [&, m](const std::vector<double>& p, std::vector<double>& r,
                               std::vector<double>& jac) {
            r.resize(n);
            jac.resize(n * k);
            const double vmax = p[0];
            for (std::size_t j = 0; j < n; ++j) {
                const Denominator d = denominatorFor(m, s[j], ic[j], p);
                const double v = vmax * s[j] / d.value;
                r[j] = v - y[j];
                const double common = -vmax * s[j] / (d.value * d.value);
                jac[j * k + 0] = s[j] / d.value;
                jac[j * k + 1] = common * d.dKm;
                jac[j * k + 2] = common * d.dKi1;
                if (mixed) jac[j * k + 3] = common * d.dKi2;
            }
        };

        std::vector<double> initial{vGuess, kmGuess, medianI};
        if (mixed) initial.push_back(medianI);

        FitResult f = runFit(in, evaluate, initial, options);
        f.modality = m;
        f.fittedX = s;
        f.seriesId = std::string("global [S]x[I]: ") + modalityLabel(m);
        out.candidates.push_back(std::move(f));
    }

    std::sort(out.candidates.begin(), out.candidates.end(),
              [](const FitResult& a, const FitResult& b) { return a.aicc < b.aicc; });

    out.deltaAicc = out.candidates.size() >= 2 ? out.candidates[1].aicc - out.candidates[0].aicc
                                               : 0.0;
    // Under 2 AICc units the data do not choose between the top two models. Saying
    // "competitive" there would be inventing a mechanism out of rounding.
    out.decisive = out.deltaAicc >= 2.0;
    const InhibitionModality winner = out.candidates.front().modality;
    if (!out.decisive) {
        out.candidates.front().modality = InhibitionModality::Unknown;
        out.candidates.front().warnings.push_back(
            "top-two AICc difference " + sci(out.deltaAicc) +
            " < 2: the modality is Unknown, not the best-ranked candidate");
        out.conclusion = std::string("inconclusive: ") + modalityLabel(winner) +
                         " ranks first but only by " + sci(out.deltaAicc) +
                         " AICc units over " +
                         modalityLabel(out.candidates[1].modality) +
                         "; modality reported as unknown";
    } else {
        out.conclusion = std::string(modalityLabel(winner)) + " inhibition, " +
                         sci(out.deltaAicc) + " AICc units ahead of " +
                         modalityLabel(out.candidates[1].modality);
    }
    return out;
}

// ---------------------------------------------------------------------------
// Inverse prediction and diagnostics
// ---------------------------------------------------------------------------

InversePrediction inversePredict(const FitResult& fit, double response,
                                 const std::vector<DosePoint>& points) {
    InversePrediction out;
    if (!fit.converged) {
        out.note = "the fit did not converge, so it cannot be inverted";
        return out;
    }
    const bool isFivePl = fit.model == AssayModel::FiveParameterLogistic;
    if (fit.model != AssayModel::FourParameterLogistic && !isFivePl) {
        out.note = "inverse prediction is defined for the 4PL and 5PL curves only";
        return out;
    }
    if (fit.parameters.size() < 4) {
        out.note = "the fit carries no parameters";
        return out;
    }

    const double a = fit.parameters[0].value.value;
    const double d = fit.parameters[1].value.value;
    const double logLoc = fit.parameters[2].value.value;
    const double b = fit.parameters[3].value.value;
    double logX = 0.0;
    if (isFivePl) {
        const double g = fit.parameters[4].value.value;
        // y = D + (A-D) (1+t)^-G  =>  t = ((A-D)/(y-D))^(1/G) - 1, t = (x/C)^B
        const double ratio = (a - d) / (response - d);
        if (!(ratio > 0.0)) {
            out.note = "the requested response lies outside the fitted asymptotes";
            return out;
        }
        const double t = std::pow(ratio, 1.0 / g) - 1.0;
        if (!(t > 0.0)) {
            out.note = "the requested response lies outside the fitted asymptotes";
            return out;
        }
        logX = logLoc + std::log10(t) / b;
    } else {
        // y = D + (A-D)/(1+u), u = 10^(B(logC - x))  =>  u = (A-y)/(y-D)
        const double u = (a - response) / (response - d);
        if (!(u > 0.0)) {
            out.note = "the requested response lies outside the fitted asymptotes";
            return out;
        }
        logX = logLoc - std::log10(u) / b;
    }
    out.concentration = std::pow(10.0, logX);

    // The horizontal uncertainty of a logistic is dominated by its location
    // parameter, so the profile interval of that parameter is translated onto the
    // requested response. This is an approximation and is labelled as one.
    const FittedParameter& loc = fit.parameters[2];
    if (loc.profileComputed) {
        out.lower = std::pow(10.0, logX + (loc.profileLower - logLoc));
        out.upper = std::pow(10.0, logX + (loc.profileUpper - logLoc));
        out.intervalFound = true;
        out.note = "interval translated from the profile-likelihood interval of " + loc.name +
                   "; approximate, it ignores slope and asymptote uncertainty";
    } else {
        out.note = "no interval: refit with FitOptions::profileLikelihood to get one";
    }

    if (!points.empty()) {
        double lo = points.front().concentration, hi = points.front().concentration;
        for (const auto& p : points) {
            lo = std::min(lo, p.concentration);
            hi = std::max(hi, p.concentration);
        }
        out.extrapolated = out.concentration < lo || out.concentration > hi;
        if (out.extrapolated) {
            out.note += (out.note.empty() ? "" : "; ");
            out.note += "outside the tested range " + sci(lo) + " to " + sci(hi) + " mol/L";
        }
    }
    return out;
}

std::vector<LineweaverBurkPoint> lineweaverBurk(const std::vector<KineticPoint>& points) {
    std::vector<LineweaverBurkPoint> out;
    out.reserve(points.size());
    for (const auto& p : points) {
        if (!(p.substrate > 0.0) || !(p.velocity > 0.0)) continue;
        out.push_back({1.0 / p.substrate, 1.0 / p.velocity});
    }
    return out;
}

// ---------------------------------------------------------------------------
// Well adapters
// ---------------------------------------------------------------------------

namespace {

// Per-concentration replicate SD, so 1/sd^2 weighting has real numbers behind it.
// A single replicate has no SD; it gets -1 (absent) rather than 0, which would
// otherwise mean "infinitely precise".
std::map<double, double> replicateSds(const std::vector<Well>& wells) {
    std::map<double, std::vector<double>> byConc;
    for (const auto& w : wells) {
        if (w.excluded) continue;
        byConc[w.concentration].push_back(w.readout);
    }
    std::map<double, double> sds;
    for (const auto& [conc, vals] : byConc) {
        if (vals.size() < 2) {
            sds[conc] = -1.0;
            continue;
        }
        double mean = 0.0;
        for (double v : vals) mean += v;
        mean /= static_cast<double>(vals.size());
        double ss = 0.0;
        for (double v : vals) ss += (v - mean) * (v - mean);
        sds[conc] = std::sqrt(ss / static_cast<double>(vals.size() - 1));
    }
    return sds;
}

}  // namespace

std::vector<DosePoint> doseSeriesFromWells(const std::vector<Well>& wells) {
    const std::map<double, double> sds = replicateSds(wells);
    std::vector<DosePoint> out;
    out.reserve(wells.size());
    for (const auto& w : wells) {
        if (w.excluded) continue;
        const auto it = sds.find(w.concentration);
        out.push_back({w.concentration, w.readout, it == sds.end() ? -1.0 : it->second});
    }
    return out;
}

std::vector<KineticPoint> kineticSeriesFromWells(const std::vector<Well>& wells) {
    const std::map<double, double> sds = replicateSds(wells);
    std::vector<KineticPoint> out;
    out.reserve(wells.size());
    for (const auto& w : wells) {
        if (w.excluded) continue;
        const auto it = sds.find(w.concentration);
        out.push_back({w.concentration, w.readout, it == sds.end() ? -1.0 : it->second});
    }
    return out;
}

std::vector<InhibitionPoint> inhibitionMatrixFromWells(const std::vector<Well>& wells) {
    std::vector<InhibitionPoint> out;
    out.reserve(wells.size());
    for (const auto& w : wells) {
        if (w.excluded) continue;
        // seriesId carries the inhibitor concentration for this column; a series id
        // that is not a number is skipped rather than silently treated as [I] = 0,
        // which would fabricate an uninhibited control.
        char* end = nullptr;
        const double inhibitor = std::strtod(w.seriesId.c_str(), &end);
        if (end == w.seriesId.c_str() || *end != '\0') continue;
        out.push_back({w.concentration, inhibitor, w.readout, -1.0});
    }
    return out;
}

FitResult fitSeries(const std::vector<Well>& series, AssayModel model,
                    const FitOptions& options) {
    switch (model) {
        case AssayModel::FourParameterLogistic:
            return fitFourParameterLogistic(doseSeriesFromWells(series), options);
        case AssayModel::FiveParameterLogistic:
            return fitFiveParameterLogistic(doseSeriesFromWells(series), options);
        case AssayModel::MichaelisMenten:
            return fitMichaelisMenten(kineticSeriesFromWells(series), options);
        case AssayModel::Hill:
            return fitHill(kineticSeriesFromWells(series), options);
        case AssayModel::SubstrateInhibition:
            return fitSubstrateInhibition(kineticSeriesFromWells(series), options);
        default:
            break;
    }
    // Morrison needs a measured [E]t and the biophysics models live in
    // assay/Biophysics.*; refusing by name beats guessing a missing input.
    return failed(model, "this model is not a plate-series fit in assay/Fits.cpp (Morrison "
                         "needs a measured [E]t; SPR, DSF and ITC live in assay/Biophysics)");
}

ModelComparison compareModels(const std::vector<Well>& series,
                              const std::vector<AssayModel>& models,
                              const FitOptions& options) {
    ModelComparison out;
    for (AssayModel m : models) out.candidates.push_back(fitSeries(series, m, options));
    std::sort(out.candidates.begin(), out.candidates.end(),
              [](const FitResult& a, const FitResult& b) {
                  if (a.converged != b.converged) return a.converged;
                  return a.aicc < b.aicc;
              });
    if (out.candidates.size() >= 2 && out.candidates[0].converged && out.candidates[1].converged) {
        out.deltaAicc = out.candidates[1].aicc - out.candidates[0].aicc;
        out.decisive = out.deltaAicc >= 2.0;
        out.conclusion = out.decisive
                             ? "the best model leads by " + sci(out.deltaAicc) + " AICc units"
                             : "AICc difference " + sci(out.deltaAicc) +
                                   " < 2: the data do not choose between the top two models";
    } else {
        out.conclusion = "fewer than two converged candidates; no comparison is possible";
    }
    return out;
}

ModelComparison fitInhibitionModalityFromWells(const std::vector<Well>& matrix,
                                               const FitOptions& options) {
    return fitInhibitionModality(inhibitionMatrixFromWells(matrix), options);
}

}  // namespace biocad::assay
