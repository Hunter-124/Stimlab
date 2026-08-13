#include "assay/Biophysics.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "numeric/Ode.h"
#include "numeric/Optimize.h"

namespace biocad::assay {
namespace {

// Gas constant in kcal/(mol K). Same value as chem::kGasConstantKcal; repeated here
// rather than linking biocad_chem into biocad_assay for one scalar.
constexpr double kGasConstantKcal = 1.987204259e-3;
constexpr double kKelvinOffset = 273.15;
// 1 kcal = 1e9 microcalories. ITC heats are entered in ucal because that is what the
// instrument integrates, while enthalpies are kcal/mol; the conversion lives in one
// place so the two never get multiplied in the wrong direction.
constexpr double kUcalPerKcal = 1e9;

// -------------------------------------------------------------------------------
// Shared helpers
// -------------------------------------------------------------------------------

double covAt(const std::vector<double>& cov, std::size_t n, std::size_t i, std::size_t j) {
    if (cov.size() != n * n) {
        return 0.0;
    }
    return cov[i * n + j];
}

double standardError(const std::vector<double>& errors, std::size_t i) {
    return i < errors.size() ? errors[i] : 0.0;
}

FittedParameter modelParameter(std::string name, double value, std::string unit, double error,
                               std::string source) {
    FittedParameter p;
    p.name = std::move(name);
    p.value = makeQuantity(value, std::move(unit), error, Provenance::Model, std::move(source));
    return p;
}

FittedParameter missingParameter(std::string name, std::string reason) {
    FittedParameter p;
    p.name = std::move(name);
    p.value = notComputed(std::move(reason));
    return p;
}

FitResult refused(const std::string& seriesId, AssayModel model, std::string note) {
    FitResult out;
    out.seriesId = seriesId;
    out.model = model;
    out.converged = false;
    out.note = std::move(note);
    out.derivedKd = notComputed(out.note);
    out.warnings.push_back(out.note);
    return out;
}

void finish(FitResult& out, const numeric::LmResult& fit, const std::vector<double>& observed,
            std::size_t parameterCount) {
    double ssr = 0.0;
    for (double r : out.residuals) {
        ssr += r * r;
    }
    out.observations = observed.size();
    out.rSquared = numeric::rSquared(observed, out.fittedY);
    out.aicc = numeric::aicc(ssr, observed.size(), parameterCount);
    out.conditionNumber = fit.conditionNumber;
    out.rank = fit.rank;
    out.converged = fit.converged;
    // The optimizer's stopping reason is appended, never substituted: the model-level
    // note (Wiseman c, the derivative Tm, the Rmax ratio) is the part a reader needs.
    if (!fit.note.empty()) {
        out.note = out.note.empty() ? fit.note : out.note + "; " + fit.note;
    }
}

// -------------------------------------------------------------------------------
// Kinetics: one integrator, two rate laws
// -------------------------------------------------------------------------------

// f = dR/dt at the given bulk concentration, plus df/dR and df/dp for each kinetic
// parameter. Written as a callback so the Langmuir and the mass-transport model share
// the integration, the sensitivity propagation and the whole result-assembly path.
using RateFn = std::function<void(double conc, double response, const double* p, double& f,
                                 double& dfdR, double* dfdp)>;

void langmuirRate(double conc, double R, const double* p, double& f, double& dfdR, double* dfdp) {
    const double ka = p[0];
    const double kd = p[1];
    const double rmax = p[2];
    f = ka * conc * (rmax - R) - kd * R;
    dfdR = -ka * conc - kd;
    dfdp[0] = conc * (rmax - R);
    dfdp[1] = -R;
    dfdp[2] = ka * conc;
}

void massTransportRate(double conc, double R, const double* p, double& f, double& dfdR,
                       double* dfdp) {
    const double ka = p[0];
    const double kd = p[1];
    const double rmax = p[2];
    const double kt = p[3];
    const double free = rmax - R;
    const double g = ka * conc * free - kd * R;      // the intrinsic 1:1 rate
    const double den = kt + ka * free;               // transport in series with binding
    const double den2 = den * den;
    f = kt * g / den;
    const double dgdR = -ka * conc - kd;
    const double ddendR = -ka;
    dfdR = kt * (dgdR * den - g * ddendR) / den2;
    dfdp[0] = kt * (conc * free * den - g * free) / den2;
    dfdp[1] = -kt * R / den;
    dfdp[2] = kt * (ka * conc * den - g * ka) / den2;
    dfdp[3] = g * ka * free / den2;
}

// Integrates R and dR/d(ln p_i) along one sensorgram, calling `sample` at every
// observed time. Sensitivities are taken with respect to LOG parameters because ka
// and kd differ by eight decades in a normal SPR run; a Jacobian in linear space is
// then numerically rank-deficient long before the experiment is.
void integrateCurve(const KineticCurve& curve, const RateFn& rate, const double* p,
                    std::size_t nKinetic, double maxStep,
                    const std::function<void(std::size_t, double, const double*)>& sample) {
    std::vector<double> y(1 + nKinetic, 0.0);   // R = 0 and zero sensitivities at t = 0
    std::vector<double> dfdp(nKinetic, 0.0);
    double conc = curve.concentrationM;

    numeric::OdeDerivative deriv = [&](double, const std::vector<double>& state,
                                       std::vector<double>& dydt) {
        double f = 0.0;
        double dfdR = 0.0;
        rate(conc, state[0], p, f, dfdR, dfdp.data());
        dydt[0] = f;
        for (std::size_t i = 0; i < nKinetic; ++i) {
            // d/dt (dR/d ln p_i) = (df/dR) * s_i + p_i * df/dp_i
            dydt[1 + i] = dfdR * state[1 + i] + p[i] * dfdp[i];
        }
    };
    const auto noObserve = [](double, const std::vector<double>&) {};

    double t = 0.0;
    for (std::size_t i = 0; i < curve.timeS.size(); ++i) {
        double target = curve.timeS[i];
        while (t < target - 1e-12) {
            double segmentEnd = target;
            // Never step across the buffer switch: the rate law changes there.
            if (t < curve.dissociationStartS - 1e-12 && curve.dissociationStartS < target) {
                segmentEnd = curve.dissociationStartS;
            }
            conc = t < curve.dissociationStartS - 1e-12 ? curve.concentrationM : 0.0;
            const double span = segmentEnd - t;
            const double substeps = std::max(1.0, std::ceil(span / std::max(maxStep, 1e-9)));
            numeric::rk4Integrate(t, segmentEnd, span / substeps, y, deriv, noObserve);
            t = segmentEnd;
        }
        sample(i, y[0], y.data() + 1);
    }
}

std::size_t countObservations(const KineticExperiment& experiment) {
    std::size_t n = 0;
    for (const KineticCurve& c : experiment.curves) {
        n += std::min(c.timeS.size(), c.responseRu.size());
    }
    return n;
}

// Steady-state analysis: Req = Rmax*C/(C + KD) over the end-of-association responses.
// A separate two-parameter fit, not a rearrangement of the kinetic constants, because
// the whole point of quoting a steady-state KD is that it is independent evidence.
struct SteadyState {
    double kd = 0;
    double kdError = 0;
    double rmax = 0;
    bool   converged = false;
};

SteadyState fitSteadyState(const std::vector<double>& conc, const std::vector<double>& req) {
    SteadyState out;
    if (conc.size() < 2) {
        return out;
    }
    const double rmaxGuess = *std::max_element(req.begin(), req.end()) * 1.2;
    const double kdGuess = *std::max_element(conc.begin(), conc.end());
    numeric::LmResult fit = numeric::levenbergMarquardt(
        {std::log(rmaxGuess), std::log(kdGuess)}, conc.size(),
        [&](const std::vector<double>& lp, std::vector<double>& r, std::vector<double>& j) {
            const double rmax = std::exp(lp[0]);
            const double kd = std::exp(lp[1]);
            for (std::size_t i = 0; i < conc.size(); ++i) {
                const double denom = conc[i] + kd;
                r[i] = rmax * conc[i] / denom - req[i];
                j[i * 2 + 0] = rmax * conc[i] / denom;                       // d/d ln Rmax
                j[i * 2 + 1] = -rmax * conc[i] * kd / (denom * denom);       // d/d ln KD
            }
        });
    out.rmax = std::exp(fit.params[0]);
    out.kd = std::exp(fit.params[1]);
    out.kdError = out.kd * standardError(fit.standardErrors, 1);
    out.converged = fit.converged;
    return out;
}

FitResult fitKinetics(const KineticExperiment& experiment, const KineticFitOptions& options,
                      AssayModel model, const RateFn& rate, std::size_t nKinetic) {
    if (experiment.curves.empty()) {
        return refused(experiment.seriesId, model, "no sensorgram supplied");
    }
    for (const KineticCurve& c : experiment.curves) {
        if (c.timeS.size() != c.responseRu.size() || c.timeS.size() < 4) {
            return refused(experiment.seriesId, model,
                           "every sensorgram needs at least 4 matched time/response points");
        }
        if (c.concentrationM <= 0.0) {
            return refused(experiment.seriesId, model,
                           "analyte concentration must be positive for every curve");
        }
    }

    const std::size_t nCurves = experiment.curves.size();
    const std::size_t nObs = countObservations(experiment);
    const bool offsets = options.fitOffsets;
    const std::size_t nParams = nKinetic + (offsets ? nCurves : 0);

    double maxObserved = 0.0;
    for (const KineticCurve& c : experiment.curves) {
        for (double r : c.responseRu) {
            maxObserved = std::max(maxObserved, r);
        }
    }
    const double rmaxInit = options.rmaxInitial > 0.0 ? options.rmaxInitial
                                                      : std::max(maxObserved * 1.1, 1e-6);

    // Starting values come from the trace, not from a constant. A fixed ka = 1e5 is a
    // guess about the interaction; the initial association slope and the dissociation
    // tail are measurements of it, and starting from them is what keeps the optimizer
    // out of the kd -> 0 / Rmax -> inf trough that a cold start falls into.
    const KineticCurve& top = *std::max_element(
        experiment.curves.begin(), experiment.curves.end(),
        [](const KineticCurve& a, const KineticCurve& b) {
            return a.concentrationM < b.concentrationM;
        });
    double kdInit = options.kdInitial;
    if (kdInit <= 0.0) {
        // ln R against t over the dissociation phase: slope = -kd for a 1:1 decay.
        double sx = 0, sy = 0, sxx = 0, sxy = 0, n = 0;
        for (std::size_t i = 0; i < top.timeS.size(); ++i) {
            if (top.timeS[i] > top.dissociationStartS + 1e-12 && top.responseRu[i] > 0.0) {
                const double x = top.timeS[i] - top.dissociationStartS;
                const double y = std::log(top.responseRu[i]);
                sx += x;
                sy += y;
                sxx += x * x;
                sxy += x * y;
                n += 1.0;
            }
        }
        const double denom = n * sxx - sx * sx;
        const double slope = n >= 3 && denom > 0.0 ? (n * sxy - sx * sy) / denom : 0.0;
        kdInit = slope < 0.0 ? -slope : 1e-3;
    }
    double kaInit = options.kaInitial;
    if (kaInit <= 0.0) {
        // dR/dt at t = 0 is ka*C*Rmax, before any appreciable occupancy.
        double slope = 0.0;
        for (std::size_t i = 1; i < top.timeS.size(); ++i) {
            if (top.timeS[i] > top.dissociationStartS) {
                break;
            }
            slope = (top.responseRu[i] - top.responseRu[0]) / (top.timeS[i] - top.timeS[0]);
            break;
        }
        kaInit = slope > 0.0 ? slope / (top.concentrationM * rmaxInit) : 1e5;
    }

    std::vector<double> initial(nParams, 0.0);
    initial[0] = std::log(kaInit);
    initial[1] = std::log(kdInit);
    initial[2] = std::log(rmaxInit);
    if (nKinetic > 3) {
        // kt only matters when it is comparable to ka*Rmax; start it there.
        initial[3] = std::log(options.ktInitial > 0.0 ? options.ktInitial : kaInit * rmaxInit);
    }

    std::vector<double> observed;
    observed.reserve(nObs);
    for (const KineticCurve& c : experiment.curves) {
        observed.insert(observed.end(), c.responseRu.begin(), c.responseRu.end());
    }

    std::vector<double> linear(nKinetic, 0.0);
    auto evaluate = [&](const std::vector<double>& params, std::vector<double>& residuals,
                        std::vector<double>& jacobian) {
        for (std::size_t i = 0; i < nKinetic; ++i) {
            linear[i] = std::exp(params[i]);
        }
        std::fill(jacobian.begin(), jacobian.end(), 0.0);
        std::size_t row = 0;
        for (std::size_t c = 0; c < nCurves; ++c) {
            const KineticCurve& curve = experiment.curves[c];
            const double offset = offsets ? params[nKinetic + c] : 0.0;
            integrateCurve(curve, rate, linear.data(), nKinetic, options.maxStepS,
                           [&](std::size_t i, double R, const double* sens) {
                               const std::size_t r = row + i;
                               residuals[r] = R + offset - curve.responseRu[i];
                               for (std::size_t k = 0; k < nKinetic; ++k) {
                                   jacobian[r * nParams + k] = sens[k];
                               }
                               if (offsets) {
                                   jacobian[r * nParams + nKinetic + c] = 1.0;
                               }
                           });
            row += curve.timeS.size();
        }
    };

    numeric::LmResult fit = numeric::levenbergMarquardt(initial, nObs, evaluate);

    const double ka = std::exp(fit.params[0]);
    const double kd = std::exp(fit.params[1]);
    const double rmax = std::exp(fit.params[2]);

    FitResult out;
    out.seriesId = experiment.seriesId;
    out.model = model;
    out.parameters.push_back(modelParameter(
        "ka", ka, "1/(M s)", ka * standardError(fit.standardErrors, 0), "global 1:1 fit"));
    out.parameters.push_back(modelParameter(
        "kd", kd, "1/s", kd * standardError(fit.standardErrors, 1), "global 1:1 fit"));
    out.parameters.push_back(modelParameter(
        "Rmax", rmax, "RU", rmax * standardError(fit.standardErrors, 2), "global 1:1 fit"));
    if (nKinetic > 3) {
        const double kt = std::exp(fit.params[3]);
        out.parameters.push_back(modelParameter(
            "kt", kt, "RU/(M s)", kt * standardError(fit.standardErrors, 3),
            "two-compartment mass transport, surface compartment at quasi steady state"));
    }
    for (std::size_t c = 0; offsets && c < nCurves; ++c) {
        const double v = fit.params[nKinetic + c];
        out.parameters.push_back(modelParameter(
            "baseline offset " + std::to_string(c), v, "RU",
            standardError(fit.standardErrors, nKinetic + c), "per-curve baseline"));
    }

    // KD = kd/ka. In log space the propagation is exact and needs no derivative:
    // var(ln KD) = var(ln kd) + var(ln ka) - 2 cov(ln ka, ln kd).
    const double kdKinetic = kd / ka;
    const double varLn = covAt(fit.covariance, nParams, 1, 1) + covAt(fit.covariance, nParams, 0, 0)
                         - 2.0 * covAt(fit.covariance, nParams, 0, 1);
    out.derivedKd = makeQuantity(kdKinetic, "M", kdKinetic * std::sqrt(std::max(varLn, 0.0)),
                                 Provenance::Model, "kinetic KD = kd/ka from the global fit");

    // Equilibrium test. Req is the response the association phase would reach at this
    // concentration; a trace that stopped at 40% of it has no plateau to read a
    // steady-state KD from, whatever a curve through the points would say.
    double worstFraction = 1.0;
    std::size_t worstCurve = 0;
    std::vector<double> ssConc;
    std::vector<double> ssReq;
    for (std::size_t c = 0; c < nCurves; ++c) {
        const KineticCurve& curve = experiment.curves[c];
        const double offset = offsets ? fit.params[nKinetic + c] : 0.0;
        double lastResponse = 0.0;
        bool found = false;
        for (std::size_t i = 0; i < curve.timeS.size(); ++i) {
            if (curve.timeS[i] <= curve.dissociationStartS + 1e-12) {
                lastResponse = curve.responseRu[i] - offset;
                found = true;
            }
        }
        if (!found) {
            continue;
        }
        const double req = rmax * curve.concentrationM / (curve.concentrationM + kdKinetic);
        const double fraction = req > 0.0 ? lastResponse / req : 0.0;
        if (fraction < worstFraction) {
            worstFraction = fraction;
            worstCurve = c;
        }
        ssConc.push_back(curve.concentrationM);
        ssReq.push_back(lastResponse);
    }

    if (worstFraction < options.equilibriumFraction) {
        char buffer[256];
        std::snprintf(buffer, sizeof buffer,
                      "steady-state KD withheld: curve %zu (%.3g M) reached only %.1f%% of the "
                      "fitted equilibrium response, below the required %.0f%% of Req",
                      worstCurve, experiment.curves[worstCurve].concentrationM,
                      100.0 * worstFraction, 100.0 * options.equilibriumFraction);
        out.parameters.push_back(missingParameter("KD (steady state)", buffer));
        out.warnings.emplace_back(buffer);
    } else {
        const SteadyState ss = fitSteadyState(ssConc, ssReq);
        if (ss.converged) {
            out.parameters.push_back(modelParameter(
                "KD (steady state)", ss.kd, "M", ss.kdError,
                "Langmuir isotherm through the end-of-association responses"));
            out.parameters.push_back(modelParameter(
                "Rmax (steady state)", ss.rmax, "RU", 0.0,
                "from the same isotherm; compare with the kinetic Rmax"));
        } else {
            out.parameters.push_back(missingParameter(
                "KD (steady state)", "steady-state isotherm fit did not converge"));
        }
    }

    if (experiment.theoreticalRmaxRu > 0.0) {
        out.parameters.push_back(modelParameter(
            "Rmax (theoretical)", experiment.theoreticalRmaxRu, "RU", 0.0,
            "from the immobilisation level and the analyte/ligand mass ratio"));
        const double ratio = rmax / experiment.theoreticalRmaxRu;
        char buffer[192];
        std::snprintf(buffer, sizeof buffer, "fitted Rmax is %.2fx the theoretical Rmax", ratio);
        out.note = buffer;
        if (ratio > 1.5 || ratio < 0.5) {
            out.warnings.emplace_back(
                std::string(buffer) +
                "; a 1:1 monolayer cannot do that, so the surface is heterogeneous or the "
                "immobilisation level is misstated");
        }
    }

    out.fittedX.reserve(nObs);
    out.fittedY.assign(nObs, 0.0);
    out.residuals.assign(nObs, 0.0);
    for (const KineticCurve& c : experiment.curves) {
        out.fittedX.insert(out.fittedX.end(), c.timeS.begin(), c.timeS.end());
    }
    {
        std::vector<double> residuals(nObs, 0.0);
        std::vector<double> jacobian(nObs * nParams, 0.0);
        evaluate(fit.params, residuals, jacobian);
        for (std::size_t i = 0; i < nObs; ++i) {
            out.residuals[i] = residuals[i];
            out.fittedY[i] = observed[i] + residuals[i];
        }
    }

    out.assumptions.emplace_back("one homogeneous 1:1 interaction; ka, kd and Rmax are global "
                                 "across the concentration series");
    out.assumptions.emplace_back("bulk analyte concentration constant during association and "
                                 "zero during dissociation");
    if (nKinetic > 3) {
        out.assumptions.emplace_back("the surface compartment is at quasi steady state, so kt is "
                                     "an instrument-scaled transport coefficient in RU/(M s)");
    }
    finish(out, fit, observed, nParams);
    return out;
}

// -------------------------------------------------------------------------------
// Melt curves
// -------------------------------------------------------------------------------

struct MeltData {
    std::vector<double> temperatureC;
    std::vector<double> signal;
    bool                truncated = false;
    double              truncationC = 0;
};

// SYPRO Orange traces peak and then fall as the dye leaves aggregating protein. The
// falling limb is dye behaviour, not a second unfolding transition, and leaving it in
// pulls the fitted unfolded baseline down and Tm with it.
MeltData prepareMelt(const MeltCurve& curve) {
    MeltData out;
    const std::size_t n = std::min(curve.temperatureC.size(), curve.signal.size());
    std::size_t last = n;
    if (curve.syproOrange && n > 0) {
        const std::size_t peak = static_cast<std::size_t>(
            std::max_element(curve.signal.begin(), curve.signal.begin() + static_cast<long>(n)) -
            curve.signal.begin());
        if (peak + 1 < n) {
            last = peak + 1;
            out.truncated = true;
            out.truncationC = curve.temperatureC[peak];
        }
    }
    out.temperatureC.assign(curve.temperatureC.begin(),
                            curve.temperatureC.begin() + static_cast<long>(last));
    out.signal.assign(curve.signal.begin(), curve.signal.begin() + static_cast<long>(last));
    return out;
}

// Savitzky-Golay first derivative, window 9, order 2. For a first derivative the
// quadratic and the linear least-squares filters give the SAME coefficients - the
// quadratic term is even and drops out of the odd derivative - so the convolution
// weights are c_i = i / (h * sum i^2), i = -4..4, sum i^2 = 60.
constexpr int kSgHalfWidth = 4;
constexpr double kSgNormalization = 60.0;

double sgDerivative(const std::vector<double>& y, std::size_t centre, double spacing) {
    double acc = 0.0;
    for (int i = -kSgHalfWidth; i <= kSgHalfWidth; ++i) {
        acc += static_cast<double>(i) * y[centre + static_cast<std::size_t>(i)];
    }
    return acc / (kSgNormalization * spacing);
}

struct Spacing {
    double mean = 0;
    double cv = 0;
    bool   uniform = false;
};

Spacing gridSpacing(const std::vector<double>& t) {
    Spacing out;
    if (t.size() < 2) {
        return out;
    }
    double sum = 0.0;
    for (std::size_t i = 1; i < t.size(); ++i) {
        sum += t[i] - t[i - 1];
    }
    out.mean = sum / static_cast<double>(t.size() - 1);
    double var = 0.0;
    for (std::size_t i = 1; i < t.size(); ++i) {
        const double d = (t[i] - t[i - 1]) - out.mean;
        var += d * d;
    }
    var /= static_cast<double>(t.size() - 1);
    out.cv = out.mean > 0.0 ? std::sqrt(var) / out.mean : 1.0;
    out.uniform = out.mean > 0.0 && out.cv <= 0.05;
    return out;
}

// Tm from the derivative extremum, refined by the vertex of the parabola through the
// peak and its two neighbours so the answer is not quantised to the sample grid.
struct DerivativePeak {
    double tm = 0;
    bool   found = false;
    bool   atEdge = false;
    double spacing = 0;
    std::string reason;
};

DerivativePeak derivativePeak(const MeltData& data) {
    DerivativePeak out;
    const std::size_t n = data.temperatureC.size();
    if (n < 2 * kSgHalfWidth + 3) {
        out.reason = "Savitzky-Golay window 9 needs at least 11 points";
        return out;
    }
    const Spacing spacing = gridSpacing(data.temperatureC);
    if (!spacing.uniform) {
        out.reason = "temperature spacing is not uniform to within 5%; the window-9 order-2 "
                     "convolution coefficients assume a uniform grid";
        return out;
    }
    out.spacing = spacing.mean;
    const std::size_t first = kSgHalfWidth;
    const std::size_t last = n - kSgHalfWidth - 1;
    std::vector<double> d(n, 0.0);
    for (std::size_t i = first; i <= last; ++i) {
        d[i] = sgDerivative(data.signal, i, spacing.mean);
    }
    std::size_t peak = first;
    for (std::size_t i = first; i <= last; ++i) {
        if (std::abs(d[i]) > std::abs(d[peak])) {
            peak = i;
        }
    }
    if (peak == first || peak == last) {
        out.atEdge = true;
        out.reason = "the derivative extremum sits at the edge of the scanned range";
    }
    double tm = data.temperatureC[peak];
    if (peak > first && peak < last) {
        const double a = d[peak - 1];
        const double b = d[peak];
        const double c = d[peak + 1];
        const double denom = a - 2.0 * b + c;
        if (std::abs(denom) > 0.0) {
            tm += 0.5 * spacing.mean * (a - c) / denom;
        }
    }
    out.tm = tm;
    out.found = true;
    return out;
}

struct MeltInitial {
    double fn = 0;
    double sn = 0;
    double fu = 0;
    double su = 0;
    double tm = 0;
    double width = 1.5;
};

MeltInitial meltInitial(const MeltData& data) {
    MeltInitial out;
    const std::size_t n = data.temperatureC.size();
    const std::size_t edge = std::max<std::size_t>(2, n / 10);
    double lowSum = 0.0;
    double highSum = 0.0;
    for (std::size_t i = 0; i < edge; ++i) {
        lowSum += data.signal[i];
        highSum += data.signal[n - 1 - i];
    }
    out.fn = lowSum / static_cast<double>(edge);
    out.fu = highSum / static_cast<double>(edge);
    const DerivativePeak peak = derivativePeak(data);
    if (peak.found) {
        out.tm = peak.tm;
    } else {
        // Fall back to the half-signal crossing, which needs no derivative.
        const double mid = 0.5 * (out.fn + out.fu);
        out.tm = data.temperatureC[n / 2];
        for (std::size_t i = 1; i < n; ++i) {
            if ((data.signal[i - 1] - mid) * (data.signal[i] - mid) <= 0.0) {
                out.tm = data.temperatureC[i];
                break;
            }
        }
    }
    return out;
}

FitResult fitMelt(const MeltCurve& curve, const MeltFitOptions& options, AssayModel model) {
    const MeltData data = prepareMelt(curve);
    const std::size_t n = data.temperatureC.size();
    if (n < 12) {
        return refused(curve.seriesId, model,
                       "a melt fit needs at least 12 usable points after truncation");
    }

    const MeltInitial guess = meltInitial(data);
    const bool twoState = model == AssayModel::TwoStateThermodynamic;
    const std::size_t nParams = 6;
    const double dCp = options.deltaCpKcalPerMolK;

    std::vector<double> initial(nParams, 0.0);
    initial[0] = guess.fn;
    initial[1] = 0.0;
    initial[2] = guess.fu;
    initial[3] = 0.0;
    initial[4] = guess.tm;
    // Boltzmann's sixth parameter is the transition width in degC; the two-state
    // model's is the van't Hoff enthalpy at Tm in kcal/mol.
    initial[5] = twoState ? 100.0 : guess.width;

    const std::vector<double>& T = data.temperatureC;
    const std::vector<double>& y = data.signal;

    auto evaluate = [&](const std::vector<double>& p, std::vector<double>& residuals,
                        std::vector<double>& jacobian) {
        const double fn = p[0];
        const double sn = p[1];
        const double fu = p[2];
        const double su = p[3];
        const double tm = p[4];
        const double sixth = p[5];
        for (std::size_t i = 0; i < n; ++i) {
            const double t = T[i];
            const double bn = fn + sn * t;
            const double bu = fu + su * t;
            const double delta = bu - bn;
            double frac = 0.0;
            double dFracdTm = 0.0;
            double dFracdSixth = 0.0;
            if (twoState) {
                const double tk = t + kKelvinOffset;
                const double tmk = tm + kKelvinOffset;
                const double dG = sixth * (1.0 - tk / tmk)
                                  - dCp * ((tmk - tk) + tk * std::log(tk / tmk));
                const double rt = kGasConstantKcal * tk;
                const double u = dG / rt;
                frac = 1.0 / (1.0 + std::exp(u));
                const double dFracdu = -frac * (1.0 - frac);
                const double ddGdH = 1.0 - tk / tmk;
                const double ddGdTm = sixth * tk / (tmk * tmk) - dCp * (1.0 - tk / tmk);
                dFracdTm = dFracdu * ddGdTm / rt;
                dFracdSixth = dFracdu * ddGdH / rt;
            } else {
                const double u = (tm - t) / sixth;
                frac = 1.0 / (1.0 + std::exp(u));
                const double g = frac * (1.0 - frac);
                dFracdTm = -g / sixth;
                dFracdSixth = g * (tm - t) / (sixth * sixth);
            }
            const double fitted = bn + frac * delta;
            residuals[i] = fitted - y[i];
            double* row = jacobian.data() + i * nParams;
            row[0] = 1.0 - frac;
            row[1] = t * (1.0 - frac);
            row[2] = frac;
            row[3] = t * frac;
            row[4] = delta * dFracdTm;
            row[5] = delta * dFracdSixth;
        }
    };

    numeric::LmResult fit = numeric::levenbergMarquardt(initial, n, evaluate);

    FitResult out;
    out.seriesId = curve.seriesId;
    out.model = model;
    // Baselines are Heuristic with an EMPTY unit: DSF signal is RFU or a 350/330
    // ratio, an arbitrary scale, so attaching a physical unit would be a lie.
    out.parameters.push_back(FittedParameter{
        "native baseline intercept",
        makeQuantity(fit.params[0], "", standardError(fit.standardErrors, 0), Provenance::Heuristic,
                     "arbitrary fluorescence scale"),
        0, 0, false});
    out.parameters.push_back(FittedParameter{
        "native baseline slope",
        makeQuantity(fit.params[1], "", standardError(fit.standardErrors, 1), Provenance::Heuristic,
                     "arbitrary fluorescence scale per degC"),
        0, 0, false});
    out.parameters.push_back(FittedParameter{
        "unfolded baseline intercept",
        makeQuantity(fit.params[2], "", standardError(fit.standardErrors, 2), Provenance::Heuristic,
                     "arbitrary fluorescence scale"),
        0, 0, false});
    out.parameters.push_back(FittedParameter{
        "unfolded baseline slope",
        makeQuantity(fit.params[3], "", standardError(fit.standardErrors, 3), Provenance::Heuristic,
                     "arbitrary fluorescence scale per degC"),
        0, 0, false});
    out.parameters.push_back(modelParameter(
        "Tm", fit.params[4], "degC", standardError(fit.standardErrors, 4),
        twoState ? "two-state thermodynamic fit" : "Boltzmann fit with sloping baselines"));
    if (twoState) {
        out.parameters.push_back(modelParameter(
            "dHm", fit.params[5], "kcal/mol", standardError(fit.standardErrors, 5),
            "van't Hoff enthalpy at Tm from the two-state fit"));
        char buffer[160];
        std::snprintf(buffer, sizeof buffer,
                      "dCp fixed at %.3g kcal/(mol K); it is not identifiable from a single "
                      "transition and is echoed, not fitted",
                      dCp);
        out.assumptions.emplace_back(buffer);
        out.assumptions.emplace_back("two-state unfolding with no intermediate and no aggregation");
    } else {
        out.parameters.push_back(modelParameter(
            "transition width", fit.params[5], "degC", standardError(fit.standardErrors, 5),
            "Boltzmann width; deliberately NOT converted into an enthalpy"));
        out.assumptions.emplace_back("the Boltzmann form is empirical: its width carries no "
                                     "thermodynamics, so no dH is reported from it");
    }
    out.assumptions.emplace_back("independent linear native and unfolded baselines");

    if (data.truncated) {
        char buffer[224];
        std::snprintf(buffer, sizeof buffer,
                      "SYPRO Orange trace truncated at its fluorescence maximum (%.2f degC): the "
                      "post-peak decay is dye behaviour, not further unfolding",
                      data.truncationC);
        out.warnings.emplace_back(buffer);
    }

    const DerivativePeak peak = derivativePeak(data);
    if (peak.found) {
        const double gap = std::abs(fit.params[4] - peak.tm);
        char buffer[224];
        std::snprintf(buffer, sizeof buffer,
                      "Savitzky-Golay (window 9, order 2) derivative Tm = %.3f degC; model Tm "
                      "differs by %.3f degC",
                      peak.tm, gap);
        out.note = out.note.empty() ? buffer : out.note + "; " + buffer;
        if (gap > options.derivativeToleranceC) {
            out.warnings.emplace_back(
                std::string(buffer) +
                " - above the stated tolerance, so the transition is not cleanly two-state");
        }
        if (peak.atEdge) {
            out.warnings.emplace_back(peak.reason);
        }
    } else {
        out.warnings.emplace_back("derivative Tm not computed: " + peak.reason);
    }

    out.fittedX = T;
    out.fittedY.assign(n, 0.0);
    out.residuals.assign(n, 0.0);
    {
        std::vector<double> residuals(n, 0.0);
        std::vector<double> jacobian(n * nParams, 0.0);
        evaluate(fit.params, residuals, jacobian);
        for (std::size_t i = 0; i < n; ++i) {
            out.residuals[i] = residuals[i];
            out.fittedY[i] = y[i] + residuals[i];
        }
    }
    finish(out, fit, y, nParams);
    return out;
}

}  // namespace

FitResult fitLangmuirKinetics(const KineticExperiment& experiment,
                              const KineticFitOptions& options) {
    return fitKinetics(experiment, options, AssayModel::LangmuirKinetics, &langmuirRate, 3);
}

FitResult fitMassTransportKinetics(const KineticExperiment& experiment,
                                   const KineticFitOptions& options) {
    // Seed from the plain 1:1 fit. The two-compartment model is the 1:1 model with one
    // extra parameter, and starting kt from a cold guess while ka, kd and Rmax are also
    // cold lets the optimizer trade transport against affinity and stall in a local
    // minimum with kt -> 0. This is the same staged approach the instrument software
    // takes, and it costs one cheap extra fit.
    KineticFitOptions seeded = options;
    const FitResult plain = fitLangmuirKinetics(experiment, options);
    if (plain.converged) {
        for (const FittedParameter& p : plain.parameters) {
            if (p.name == "ka") {
                seeded.kaInitial = p.value.value;
            } else if (p.name == "kd") {
                seeded.kdInitial = p.value.value;
            } else if (p.name == "Rmax") {
                seeded.rmaxInitial = p.value.value;
            }
        }
    }
    return fitKinetics(experiment, seeded, AssayModel::MassTransportKinetics, &massTransportRate,
                       4);
}

FitResult fitBoltzmannMelt(const MeltCurve& curve, const MeltFitOptions& options) {
    return fitMelt(curve, options, AssayModel::BoltzmannMelt);
}

FitResult fitTwoStateMelt(const MeltCurve& curve, const MeltFitOptions& options) {
    return fitMelt(curve, options, AssayModel::TwoStateThermodynamic);
}

Quantity derivativeTm(const MeltCurve& curve) {
    const MeltData data = prepareMelt(curve);
    const DerivativePeak peak = derivativePeak(data);
    if (!peak.found) {
        return notComputed(peak.reason);
    }
    std::string source = "Savitzky-Golay first derivative, window 9, order 2";
    if (data.truncated) {
        source += "; trace truncated at its fluorescence maximum";
    }
    // The error bar is half a grid step: that is the resolution of a peak position
    // read off a sampled derivative, parabolic refinement notwithstanding.
    return makeQuantity(peak.tm, "degC", 0.5 * peak.spacing, Provenance::Measured,
                        std::move(source));
}

// -------------------------------------------------------------------------------
// ITC
// -------------------------------------------------------------------------------

FitResult fitWisemanIsotherm(const ItcExperiment& experiment) {
    const std::size_t n = experiment.injections.size();
    if (n < 4) {
        return refused(experiment.seriesId, AssayModel::WisemanIsotherm,
                       "an isotherm fit needs at least 4 injections");
    }
    if (experiment.cellVolumeL <= 0.0 || experiment.macromoleculeM <= 0.0
        || experiment.titrantM <= 0.0) {
        return refused(experiment.seriesId, AssayModel::WisemanIsotherm,
                       "cell volume, cell macromolecule concentration and syringe titrant "
                       "concentration must all be positive");
    }
    if (experiment.blankHeatUcal.size() != n) {
        return refused(experiment.seriesId, AssayModel::WisemanIsotherm,
                       "missing blank heat-of-dilution run: one blank heat per injection is "
                       "required, because the dilution heat is the same order as the late "
                       "injections that determine n and K");
    }

    const double v0 = experiment.cellVolumeL;
    const double m0 = experiment.macromoleculeM;
    const double x0 = experiment.titrantM;
    const double temperature = experiment.temperatureK;

    // Displaced-volume correction for an overfilled cell: injected volume pushes cell
    // contents out through the access tube, which is an exponential dilution of the
    // macromolecule and an exponential approach to x0 for the titrant.
    std::vector<double> mt(n, 0.0);
    std::vector<double> xt(n, 0.0);
    std::vector<double> dv(n, 0.0);
    std::vector<double> observed(n, 0.0);
    double injected = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        dv[i] = experiment.injections[i].volumeL;
        injected += dv[i];
        const double dilution = std::exp(-injected / v0);
        mt[i] = m0 * dilution;
        xt[i] = x0 * (1.0 - dilution);
        observed[i] = experiment.injections[i].heatUcal - experiment.blankHeatUcal[i];
    }

    const std::size_t nParams = 3;   // n, ln K, dH
    std::vector<double> initial{1.0, std::log(1e6), -10.0};

    // Total heat after injection i (kcal), Wiseman one set of sites, and its exact
    // derivatives. K is fitted as ln K so the Jacobian column is scale-free.
    auto totalHeat = [&](std::size_t i, const std::vector<double>& p, double* dQ) {
        const double nSites = p[0];
        const double k = std::exp(p[1]);
        const double dH = p[2];
        const double xr = xt[i] / mt[i];
        const double a = 0.5 * mt[i] * dH * v0;
        const double b = 1.0 + xr / nSites + 1.0 / (nSites * k * mt[i]);
        const double disc = std::max(b * b - 4.0 * xr / nSites, 0.0);
        const double root = std::sqrt(disc);
        const double f = b - root;
        const double q = a * nSites * f;
        // dB/dn and dB/dK, then dF through the square root.
        const double dBdn = -(xr + 1.0 / (k * mt[i])) / (nSites * nSites);
        const double dDdn = 2.0 * b * dBdn + 4.0 * xr / (nSites * nSites);
        const double dFdn = root > 0.0 ? dBdn - dDdn / (2.0 * root) : dBdn;
        const double dBdK = -1.0 / (nSites * k * k * mt[i]);
        const double dFdK = root > 0.0 ? dBdK - (2.0 * b * dBdK) / (2.0 * root) : dBdK;
        dQ[0] = a * (f + nSites * dFdn);
        dQ[1] = a * nSites * dFdK * k;     // chain rule for ln K
        dQ[2] = q / dH;                    // Q is linear in dH
        return q;
    };

    auto evaluate = [&](const std::vector<double>& p, std::vector<double>& residuals,
                        std::vector<double>& jacobian) {
        double prevQ = 0.0;
        double prevD[3] = {0.0, 0.0, 0.0};
        for (std::size_t i = 0; i < n; ++i) {
            double d[3] = {0.0, 0.0, 0.0};
            const double q = totalHeat(i, p, d);
            // Observed injection heat with the displaced-volume term: the fraction
            // dv/V0 of the cell contents leaves at the average of the pre- and
            // post-injection heat content.
            const double w = dv[i] / v0;
            const double modelUcal = (q - prevQ + w * 0.5 * (q + prevQ)) * kUcalPerKcal;
            residuals[i] = modelUcal - observed[i];
            for (std::size_t j = 0; j < nParams; ++j) {
                jacobian[i * nParams + j] =
                    (d[j] - prevD[j] + w * 0.5 * (d[j] + prevD[j])) * kUcalPerKcal;
            }
            prevQ = q;
            for (std::size_t j = 0; j < nParams; ++j) {
                prevD[j] = d[j];
            }
        }
    };

    numeric::LmResult fit = numeric::levenbergMarquardt(initial, n, evaluate);

    const double nSites = fit.params[0];
    const double k = std::exp(fit.params[1]);
    const double dH = fit.params[2];
    const double seN = standardError(fit.standardErrors, 0);
    const double seLnK = standardError(fit.standardErrors, 1);
    const double seH = standardError(fit.standardErrors, 2);

    FitResult out;
    out.seriesId = experiment.seriesId;
    out.model = AssayModel::WisemanIsotherm;

    // c FIRST. c decides whether this titration could have determined K at all: below
    // ~1 the isotherm is a straight line with no curvature to fit K to, above ~1000 it
    // is a step whose position gives n and whose corner gives nothing.
    const double c = nSites * k * m0;   // Mt(0), the cell concentration before injection 1
    const double cError = c * std::sqrt(seLnK * seLnK + (nSites > 0 ? (seN / nSites) * (seN / nSites)
                                                                   : 0.0));
    out.parameters.push_back(modelParameter("c (Wiseman, n*K*Mt)", c, "", cError,
                                            "n * K * Mt(0); read before the parameters"));
    out.parameters.push_back(
        modelParameter("n", nSites, "", seN, "sites per macromolecule, one set of sites"));
    out.parameters.push_back(modelParameter("K", k, "1/M", k * seLnK, "association constant"));
    out.parameters.push_back(modelParameter("dH", dH, "kcal/mol", seH, "binding enthalpy"));

    const double rt = kGasConstantKcal * temperature;
    const double dG = -rt * fit.params[1];
    const double seG = rt * seLnK;
    out.parameters.push_back(modelParameter("dG", dG, "kcal/mol", seG, "-RT ln K at the run "
                                                                      "temperature"));
    const double minusTdS = dG - dH;
    const double varTdS = rt * rt * covAt(fit.covariance, nParams, 1, 1)
                          + covAt(fit.covariance, nParams, 2, 2)
                          + 2.0 * rt * covAt(fit.covariance, nParams, 1, 2);
    out.parameters.push_back(modelParameter("-T dS", minusTdS, "kcal/mol",
                                            std::sqrt(std::max(varTdS, 0.0)), "dG - dH"));

    out.derivedKd = makeQuantity(1.0 / k, "M", seLnK / k, Provenance::Model,
                                 "1/K from the Wiseman isotherm");

    {
        char buffer[192];
        std::snprintf(buffer, sizeof buffer, "Wiseman c = n*K*Mt = %.4g at %.2f K", c, temperature);
        out.note = buffer;
    }
    if (c < 1.0 || c > 1000.0) {
        char buffer[256];
        std::snprintf(buffer, sizeof buffer,
                      "Wiseman c = %.4g is outside c ~ 1-1000: %s, so K from this titration is "
                      "not determined by the data even if the fit converged",
                      c,
                      c < 1.0 ? "the isotherm has no curvature to fit K to"
                              : "the isotherm is a step whose corner carries no K information");
        out.warnings.emplace_back(buffer);
    }

    out.assumptions.emplace_back("one set of identical, independent sites");
    out.assumptions.emplace_back("displaced-volume correction for an overfilled cell: "
                                 "Mt = Mt(0)*exp(-V/V0), Xt = Xt(0)*(1 - exp(-V/V0))");
    out.assumptions.emplace_back("blank heat of dilution subtracted per injection from a "
                                 "titrant-into-buffer run");

    out.fittedX.assign(n, 0.0);
    out.fittedY.assign(n, 0.0);
    out.residuals.assign(n, 0.0);
    {
        std::vector<double> residuals(n, 0.0);
        std::vector<double> jacobian(n * nParams, 0.0);
        evaluate(fit.params, residuals, jacobian);
        for (std::size_t i = 0; i < n; ++i) {
            out.fittedX[i] = xt[i] / mt[i];   // molar ratio, the conventional abscissa
            out.residuals[i] = residuals[i];
            out.fittedY[i] = observed[i] + residuals[i];
        }
    }
    finish(out, fit, observed, nParams);
    return out;
}

}  // namespace biocad::assay
