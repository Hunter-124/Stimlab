#include "pkpd/Fits.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include "numeric/Optimize.h"

namespace biocad::pkpd {
namespace {

constexpr double kLn10 = 2.302585092994045684;


// Molar quantities span many decades, and std::to_string prints 5e-10 as "0.000000",
// which would hide the very number the tight-binding note exists to show.
std::string sci(double v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.4g", v);
    return buf;
}

std::string pointCount(std::size_t n) { return "n=" + std::to_string(n) + " points"; }

// Every Quantity NotComputed, with one shared reason. Used by the precondition guards
// so a failed fit can never be mistaken for a fit that happened to return zeros.
CurveFit failedFit(const std::string& reason) {
    CurveFit f;
    f.top = notComputed(reason);
    f.bottom = notComputed(reason);
    f.ec50 = notComputed(reason);
    f.hillSlope = notComputed(reason);
    f.converged = false;
    f.note = reason;
    return f;
}

SchildResult failedSchild(const std::string& reason) {
    SchildResult r;
    r.pA2 = notComputed(reason);
    r.slope = notComputed(reason);
    r.kb = notComputed(reason);
    r.kbUsable = false;
    r.note = reason;
    return r;
}

// Two-sided 95% Student t quantiles for df 1..30, then 1.96 beyond. A short table is
// honest here: the incomplete beta function is not otherwise needed anywhere in the
// tree, and for df > 30 the normal quantile is within 2% of the exact value.
double tQuantile95(int df) {
    static const double kTable[31] = {
        0.0,
        12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228,
        2.201,  2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086,
        2.080,  2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042};
    if (df < 1) return 0.0;
    if (df <= 30) return kTable[df];
    return 1.96;
}

}  // namespace

CurveFit fitFourParameterLogistic(const std::vector<DoseResponsePoint>& points,
                                  bool inverseSquareWeighting) {
    // Four parameters means four observations is the bare minimum for a determined system.
    if (points.size() < 4) {
        return failedFit("four-parameter logistic needs at least 4 points, got " +
                         std::to_string(points.size()));
    }
    for (const auto& p : points) {
        if (!(p.concentration > 0.0)) {
            return failedFit("a non-positive concentration cannot be log-transformed");
        }
    }

    const std::size_t n = points.size();
    std::vector<double> x(n), y(n);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = std::log10(points[i].concentration);
        y[i] = points[i].effect;
    }

    std::vector<double> weights;
    std::size_t skipped = 0;
    if (inverseSquareWeighting) {
        // 1/y^2 is the constant-relative-error weighting. A non-positive response has no
        // defined relative error, so it gets weight 0 (excluded) and is counted in the note
        // rather than silently reweighted to 1.
        weights.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            if (y[i] > 0.0) {
                weights[i] = 1.0 / (y[i] * y[i]);
            } else {
                weights[i] = 0.0;
                ++skipped;
            }
        }
    }

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

    // params = [Top, Bottom, log10 EC50, nH]
    std::vector<double> initial{yMax, yMin, logEc50Guess, 1.0};

    // With u = (EC50/[A])^nH = 10^(nH*(log10 EC50 - x)):
    //   E        = Bottom + (Top - Bottom) / (1 + u)
    //   dE/dTop      = 1/(1+u)
    //   dE/dBottom   = u/(1+u)
    //   dE/dlogEC50  = -delta * nH * ln10 * u / (1+u)^2
    //   dE/dnH       = -delta * ln10 * (log10 EC50 - x) * u / (1+u)^2
    // The last two carry the minus sign because raising EC50 or steepening the slope
    // pushes the curve to the right, which lowers E at a fixed x on an ascending curve.
    auto evaluate = [&](const std::vector<double>& p, std::vector<double>& residuals,
                        std::vector<double>& jacobian) {
        const double top = p[0];
        const double bottom = p[1];
        const double logEc50 = p[2];
        const double nH = p[3];
        const double delta = top - bottom;
        residuals.resize(n);
        jacobian.resize(n * 4);
        for (std::size_t i = 0; i < n; ++i) {
            const double shift = logEc50 - x[i];
            const double u = std::pow(10.0, nH * shift);
            const double denom = 1.0 + u;
            const double model = bottom + delta / denom;
            residuals[i] = model - y[i];
            const double common = -delta * u / (denom * denom);
            jacobian[i * 4 + 0] = 1.0 / denom;
            jacobian[i * 4 + 1] = u / denom;
            jacobian[i * 4 + 2] = common * nH * kLn10;
            jacobian[i * 4 + 3] = common * kLn10 * shift;
        }
    };

    const numeric::LmResult lm = numeric::levenbergMarquardt(initial, n, evaluate, weights);

    CurveFit fit;
    fit.iterations = lm.iterations;
    fit.converged = lm.converged;

    // The fitter never sees the raw observations, so R^2 is computed here from the model
    // evaluated at the returned parameters.
    std::vector<double> fitted(n);
    {
        const double top = lm.params[0];
        const double bottom = lm.params[1];
        const double logEc50 = lm.params[2];
        const double nH = lm.params[3];
        for (std::size_t i = 0; i < n; ++i) {
            const double u = std::pow(10.0, nH * (logEc50 - x[i]));
            fitted[i] = bottom + (top - bottom) / (1.0 + u);
        }
    }
    fit.rSquared = numeric::rSquared(y, fitted);

    const bool haveErrors = lm.standardErrors.size() == 4;
    const double seTop = haveErrors ? lm.standardErrors[0] : 0.0;
    const double seBottom = haveErrors ? lm.standardErrors[1] : 0.0;
    const double seLogEc50 = haveErrors ? lm.standardErrors[2] : 0.0;
    const double seSlope = haveErrors ? lm.standardErrors[3] : 0.0;

    const std::string source = "4PL fit, " + pointCount(n);
    // Top and Bottom live on the assay's own arbitrary response scale, so they are
    // Heuristic with an empty unit; makeQuantity() would throw on any unit here.
    fit.top = makeQuantity(lm.params[0], "", seTop, Provenance::Heuristic,
                           source + "; assay response scale, arbitrary units");
    fit.bottom = makeQuantity(lm.params[1], "", seBottom, Provenance::Heuristic,
                              source + "; assay response scale, arbitrary units");

    const double ec50 = std::pow(10.0, lm.params[2]);
    // Error propagation of EC50 = 10^p: dEC50 = ln10 * EC50 * d(log10 EC50).
    fit.ec50 = makeQuantity(ec50, "mol/L", kLn10 * ec50 * seLogEc50, Provenance::Predicted,
                            source + "; EC50 from the fitted log10 EC50");
    // This is an EMPIRICAL SLOPE and nothing more: in a functional assay signal
    // amplification and receptor reserve bend it independently of binding stoichiometry,
    // so a value above 1 is not evidence about the number of binding sites.
    fit.hillSlope = makeQuantity(lm.params[3], "", seSlope, Provenance::Predicted,
                                 source + "; empirical slope (dimensionless)");

    fit.note = lm.note;
    if (inverseSquareWeighting) {
        fit.note += fit.note.empty() ? "" : "; ";
        fit.note += "1/y^2 weighting";
        if (skipped > 0) {
            fit.note += " (" + std::to_string(skipped) +
                        " non-positive response(s) excluded: undefined weight)";
        }
    }
    if (!haveErrors) {
        fit.note += fit.note.empty() ? "" : "; ";
        fit.note += "standard errors unavailable (rank-deficient Jacobian)";
    }
    return fit;
}

Quantity kiFromIc50(const ChengPrusoffInput& in) {
    if (!(in.ic50 > 0.0)) return notComputed("ic50 (must be a positive molar IC50)");

    double ki = 0.0;
    std::string equation;
    switch (in.modality) {
        case InhibitionModality::Competitive: {
            if (in.substrate < 0.0) return notComputed("substrate [S] (competitive Cheng-Prusoff)");
            if (in.km < 0.0) return notComputed("km (competitive Cheng-Prusoff)");
            if (!(in.km > 0.0)) return notComputed("km (must be positive)");
            ki = in.ic50 / (1.0 + in.substrate / in.km);
            equation = "competitive: Ki = IC50 / (1 + [S]/Km)";
            break;
        }
        case InhibitionModality::RadioligandBinding: {
            if (in.radioligand < 0.0) return notComputed("radioligand [L*] (Cheng-Prusoff binding)");
            if (in.kdRadioligand < 0.0) return notComputed("kdRadioligand (Cheng-Prusoff binding)");
            if (!(in.kdRadioligand > 0.0)) return notComputed("kdRadioligand (must be positive)");
            ki = in.ic50 / (1.0 + in.radioligand / in.kdRadioligand);
            equation = "radioligand binding: Ki = IC50 / (1 + [L*]/Kd)";
            break;
        }
        case InhibitionModality::Noncompetitive: {
            // A noncompetitive inhibitor binds equally well with and without substrate,
            // so IC50 is already substrate-independent and needs no correction.
            ki = in.ic50;
            equation = "noncompetitive: Ki = IC50 (no substrate correction)";
            break;
        }
        case InhibitionModality::Uncompetitive: {
            if (in.substrate < 0.0) return notComputed("substrate [S] (uncompetitive Cheng-Prusoff)");
            if (in.km < 0.0) return notComputed("km (uncompetitive Cheng-Prusoff)");
            if (!(in.substrate > 0.0)) return notComputed("substrate [S] (must be positive)");
            ki = in.ic50 / (1.0 + in.km / in.substrate);
            equation = "uncompetitive: Ki = IC50 / (1 + Km/[S])";
            break;
        }
    }

    std::string source = equation;
    if (in.enzymeConc >= 0.0 && ki <= 10.0 * in.enzymeConc) {
        // Tight binding: the classic equation assumes free [I] equals total [I], i.e.
        // [I] >> [E]t. Once Ki approaches [E]t roughly half the enzyme sequesters the
        // inhibitor, so the Morrison depletion correction Ki - [E]t/2 is shown beside it.
        const double corrected = std::max(ki - 0.5 * in.enzymeConc, 1e-15);
        source += "; TIGHT BINDING: the classic value assumes [I] >> [E]t, which fails here"
                  " ([E]t = " + sci(in.enzymeConc) + " mol/L); classic Ki = " +
                  sci(ki) + " mol/L, depletion-corrected Ki - [E]t/2 = " +
                  sci(corrected) + " mol/L";
    }
    return makeQuantity(ki, "mol/L", 0.0, Provenance::Predicted, source);
}

SchildResult schild(const std::vector<SchildPoint>& points) {
    if (points.size() < 3) {
        return failedSchild("Schild regression needs at least 3 points, got " +
                            std::to_string(points.size()));
    }
    for (const auto& p : points) {
        if (!(p.antagonist > 0.0)) {
            return failedSchild("a non-positive antagonist concentration cannot be log-transformed");
        }
        if (!(p.doseRatio > 1.0)) {
            return failedSchild("every dose ratio must exceed 1: log10(DR - 1) is undefined otherwise");
        }
    }

    const std::size_t n = points.size();
    double sx = 0.0, sy = 0.0;
    std::vector<double> bx(n), by(n);
    for (std::size_t i = 0; i < n; ++i) {
        bx[i] = std::log10(points[i].antagonist);
        by[i] = std::log10(points[i].doseRatio - 1.0);
        sx += bx[i];
        sy += by[i];
    }
    const double mx = sx / static_cast<double>(n);
    const double my = sy / static_cast<double>(n);
    double sxx = 0.0, sxy = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        sxx += (bx[i] - mx) * (bx[i] - mx);
        sxy += (bx[i] - mx) * (by[i] - my);
    }
    if (!(sxx > 0.0)) {
        return failedSchild("all antagonist concentrations are identical: the slope is undefined");
    }

    const double slope = sxy / sxx;
    const double intercept = my - slope * mx;

    double sse = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double r = by[i] - (intercept + slope * bx[i]);
        sse += r * r;
    }
    const int df = static_cast<int>(n) - 2;
    const double slopeSe = df > 0 ? std::sqrt((sse / static_cast<double>(df)) / sxx) : 0.0;
    const double halfWidth = tQuantile95(df) * slopeSe;

    SchildResult out;
    const std::string source = "Schild regression of log10(DR-1) on log10[B], " + pointCount(n);
    // x-intercept -intercept/slope is log10 of the [B] giving DR = 2; pA2 is its negative.
    const double pA2 = intercept / slope;
    out.pA2 = makeQuantity(pA2, "", 0.0, Provenance::Predicted, source + "; pA2 = -log10[B] at DR = 2");
    out.slope = makeQuantity(slope, "", slopeSe, Provenance::Predicted,
                             source + "; slope with a 95% CI from n-2 degrees of freedom");
    out.slopeCiLow = slope - halfWidth;
    out.slopeCiHigh = slope + halfWidth;

    const bool unitSlope = out.slopeCiLow <= 1.0 && out.slopeCiHigh >= 1.0;
    if (unitSlope) {
        out.kbUsable = true;
        out.kb = makeQuantity(std::pow(10.0, -pA2), "mol/L", 0.0, Provenance::Predicted,
                              source + "; KB = 10^(-pA2), valid because the slope CI includes 1");
        out.note = "slope CI includes 1: consistent with simple competitive antagonism";
    } else {
        out.kbUsable = false;
        out.kb = notComputed("unit Schild slope");
        out.note = "slope CI [" + sci(out.slopeCiLow) + ", " +
                   sci(out.slopeCiHigh) +
                   "] excludes 1: the antagonism is not simple competitive, so a KB from this"
                   " plot would be meaningless";
    }
    return out;
}

}  // namespace biocad::pkpd
