#include "sim/Nca.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace biocad::sim {
namespace {

constexpr const char* kAucSource = "linear-up/log-down trapezoid over the observed series";

// Log-linear ("log-down") trapezoid over one falling interval:
//   integral = (c1 - c2) * dt / ln(c1/c2)
// which is EXACT for a monoexponential decay, whereas the linear trapezoid always
// overestimates it. The rule falls back to linear whenever the segment rises or an
// endpoint is non-positive, because the logarithm is undefined there.
double segmentArea(double t1, double c1, double t2, double c2) {
    const double dt = t2 - t1;
    if (dt <= 0) return 0.0;
    if (c2 < c1 && c1 > 0 && c2 > 0) return dt * (c1 - c2) / std::log(c1 / c2);
    return dt * 0.5 * (c1 + c2);
}

// First moment over the same interval, with the matching log-down form
//   integral t*c dt = dt*(t1*c1 - t2*c2)/ln(c1/c2) + dt^2*(c1 - c2)/ln(c1/c2)^2.
double segmentMoment(double t1, double c1, double t2, double c2) {
    const double dt = t2 - t1;
    if (dt <= 0) return 0.0;
    if (c2 < c1 && c1 > 0 && c2 > 0) {
        const double l = std::log(c1 / c2);
        return dt * (t1 * c1 - t2 * c2) / l + dt * dt * (c1 - c2) / (l * l);
    }
    return dt * 0.5 * (t1 * c1 + t2 * c2);
}

struct TerminalFit {
    bool   ok = false;
    double lambdaZ = 0;
    double intercept = 0;   // ln C0 of the terminal line
    double adjustedR2 = 0;
    int    points = 0;
    double stdError = 0;    // of lambda_z
};

// Regress ln(C) on t over [begin, end). Reported adjusted R-squared is
// 1 - (1 - R^2)*(k-1)/(k-2): one predictor, so a point that does not improve the
// fit more than its degree of freedom costs LOWERS the criterion.
TerminalFit regress(const std::vector<double>& t, const std::vector<double>& c,
                    std::size_t begin, std::size_t end) {
    TerminalFit f;
    const std::size_t k = end - begin;
    if (k < 3) return f;
    double sx = 0, sy = 0;
    for (std::size_t i = begin; i < end; ++i) {
        if (c[i] <= 0) return f;
        sx += t[i];
        sy += std::log(c[i]);
    }
    const double mx = sx / static_cast<double>(k), my = sy / static_cast<double>(k);
    double sxx = 0, sxy = 0, syy = 0;
    for (std::size_t i = begin; i < end; ++i) {
        const double dx = t[i] - mx, dy = std::log(c[i]) - my;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    if (sxx <= 0 || syy <= 0) return f;
    const double slope = sxy / sxx;
    if (slope >= 0) return f;   // a rising terminal phase is not an elimination slope
    const double r2 = (sxy * sxy) / (sxx * syy);
    f.ok = true;
    f.lambdaZ = -slope;
    f.intercept = my - slope * mx;
    f.points = static_cast<int>(k);
    f.adjustedR2 = 1.0 - (1.0 - r2) * (static_cast<double>(k) - 1.0) /
                             (static_cast<double>(k) - 2.0);
    const double residual = std::max(syy - slope * sxy, 0.0);
    const double s2 = residual / (static_cast<double>(k) - 2.0);
    f.stdError = std::sqrt(s2 / sxx);
    return f;
}

}  // namespace

NcaResult noncompartmental(const ConcentrationSeries& observed) {
    NcaResult r;
    r.subjectId = observed.subjectId;

    const auto& t = observed.timeH;
    const auto& c = observed.concentration;
    if (t.size() != c.size() || t.size() < 2) {
        r.cmax = notComputed("at least two concentration-time points");
        r.tmax = notComputed("at least two concentration-time points");
        r.aucLast = notComputed("at least two concentration-time points");
        r.aucInfinity = notComputed("at least two concentration-time points");
        r.warnings.push_back("the series has fewer than two usable points");
        return r;
    }
    for (std::size_t i = 1; i < t.size(); ++i) {
        if (t[i] <= t[i - 1]) {
            r.warnings.push_back("the time column is not strictly increasing; NCA needs an "
                                 "ordered series and the result is not computed");
            r.cmax = notComputed("a strictly increasing time column");
            return r;
        }
    }

    // Cmax/Tmax are DATA, not a fit: the observed maximum, with no interpolation.
    std::size_t iMax = 0;
    for (std::size_t i = 1; i < c.size(); ++i)
        if (c[i] > c[iMax]) iMax = i;
    r.cmax = makeQuantity(c[iMax], "mg/L", 0, Provenance::Measured,
                          "observed maximum concentration");
    r.tmax = makeQuantity(t[iMax], "h", 0, Provenance::Measured,
                          "time of the observed maximum, not interpolated");

    double aucLast = 0, aumcLast = 0;
    for (std::size_t i = 1; i < t.size(); ++i) {
        aucLast += segmentArea(t[i - 1], c[i - 1], t[i], c[i]);
        aumcLast += segmentMoment(t[i - 1], c[i - 1], t[i], c[i]);
    }
    r.aucLast = makeQuantity(aucLast, "mg*h/L", 0, Provenance::Measured, kAucSource);

    // lambda_z: every window that starts STRICTLY after Tmax and holds at least
    // three points is tried, and the best adjusted R-squared wins.
    TerminalFit best;
    const std::size_t n = t.size();
    for (std::size_t begin = iMax + 1; begin + 3 <= n; ++begin) {
        const TerminalFit f = regress(t, c, begin, n);
        if (f.ok && (!best.ok || f.adjustedR2 > best.adjustedR2)) best = f;
    }

    const double tLast = t.back();
    if (!best.ok) {
        const char* missing = "at least three positive points strictly after Tmax";
        r.lambdaZ = notComputed(missing);
        r.halfLife = notComputed(missing);
        r.aucInfinity = notComputed(missing);
        r.percentExtrapolated = notComputed(missing);
        r.clearance = notComputed(missing);
        r.volumeZ = notComputed(missing);
        r.volumeSteadyState = notComputed(missing);
        r.aumc = makeQuantity(aumcLast, "mg*h^2/L", 0, Provenance::Measured, kAucSource);
        r.meanResidenceTime = notComputed(missing);
        r.warnings.push_back("no terminal slope could be fitted: NCA requires three or more "
                             "positive concentrations strictly after Tmax, so half-life and "
                             "everything extrapolated from it are not computed");
    } else {
        const double lz = best.lambdaZ;
        // Clast is the REGRESSION-PREDICTED value, not the last observed point: the
        // last point carries the most assay noise, and using it makes the whole
        // extrapolated tail inherit that one measurement's error.
        const double clastPred = std::exp(best.intercept - lz * tLast);
        const double tail = clastPred / lz;
        const double aucInf = aucLast + tail;
        const double pctExtrap = 100.0 * tail / aucInf;
        const double aumcInf = aumcLast + tLast * clastPred / lz + clastPred / (lz * lz);

        r.lambdaZPointCount = best.points;
        r.adjustedRSquared = best.adjustedR2;
        r.lambdaZ = makeQuantity(lz, "1/h", best.stdError, Provenance::Measured,
                                 "log-linear regression over the " +
                                     std::to_string(best.points) +
                                     " terminal points after Tmax, chosen by adjusted "
                                     "R-squared = " +
                                     std::to_string(best.adjustedR2));
        r.halfLife = makeQuantity(std::log(2.0) / lz, "h", 0, Provenance::Measured,
                                  "ln(2)/lambda_z");
        r.aucInfinity = makeQuantity(aucInf, "mg*h/L", 0, Provenance::Measured,
                                     "AUClast + Clast,pred/lambda_z");
        r.percentExtrapolated = makeQuantity(pctExtrap, "%", 0, Provenance::Measured,
                                             "100 * (Clast,pred/lambda_z) / AUCinf");
        r.aumc = makeQuantity(aumcInf, "mg*h^2/L", 0, Provenance::Measured,
                              "AUMClast + tlast*Clast,pred/lambda_z + Clast,pred/lambda_z^2");
        r.meanResidenceTime = makeQuantity(aumcInf / aucInf, "h", 0, Provenance::Measured,
                                           "AUMCinf/AUCinf");
        r.extrapolationUnreliable = pctExtrap > 20.0;

        // Exposure-normalised parameters. The formula is the same with and without
        // an intravenous route; the LABEL is not, because without an IV reference
        // the F cannot be separated and CL/F is what was actually computed.
        const bool ss = observed.tauH > 0;
        double clearanceDenominator = aucInf;
        std::string clSource = "dose/AUCinf";
        if (ss) {
            // At steady state the interval area is the exposure per dose, and using
            // AUCinf there would double-count the accumulated tail.
            double aucTau = 0;
            double cMin = c.front();
            for (std::size_t i = 1; i < t.size() && t[i] <= observed.tauH; ++i)
                aucTau += segmentArea(t[i - 1], c[i - 1], t[i], c[i]);
            for (std::size_t i = 0; i < c.size() && t[i] <= observed.tauH; ++i)
                cMin = std::min(cMin, c[i]);
            if (aucTau > 0) {
                r.aucTau = makeQuantity(aucTau, "mg*h/L", 0, Provenance::Measured,
                                        "linear-up/log-down over one dosing interval");
                r.cAverage = makeQuantity(aucTau / observed.tauH, "mg/L", 0,
                                          Provenance::Measured, "AUCtau/tau");
                if (cMin > 0)
                    r.swing = makeQuantity((r.cmax.value - cMin) / cMin, "", 0,
                                           Provenance::Measured, "(Cmax - Cmin)/Cmin");
                else
                    r.swing = notComputed("a positive trough concentration");
                clearanceDenominator = aucTau;
                clSource = "dose/AUCtau at steady state";
            } else {
                r.aucTau = notComputed("samples covering one full dosing interval");
                r.cAverage = notComputed("samples covering one full dosing interval");
                r.swing = notComputed("samples covering one full dosing interval");
            }
        } else {
            r.aucTau = notComputed("a steady-state dosing interval (tau)");
            r.cAverage = notComputed("a steady-state dosing interval (tau)");
            r.swing = notComputed("a steady-state dosing interval (tau)");
        }

        if (observed.dose > 0 && clearanceDenominator > 0) {
            const double cl = observed.dose / clearanceDenominator;
            const char* label = observed.intravenous ? "CL" : "CL/F";
            r.clearance = makeQuantity(cl, "L/h", 0, Provenance::Measured,
                                       std::string(label) + " = " + clSource);
            r.volumeZ = makeQuantity(cl / lz, "L", 0, Provenance::Measured,
                                     std::string(observed.intravenous ? "Vz" : "Vz/F") +
                                         " = CL/lambda_z");
            if (observed.intravenous)
                // Vss = CL * MRT is only a volume when the whole dose entered the
                // systemic circulation; after an extravascular dose MRT contains the
                // absorption time and the product is not a volume at all.
                r.volumeSteadyState = makeQuantity(cl * r.meanResidenceTime.value, "L", 0,
                                                   Provenance::Measured, "Vss = CL * MRT");
            else
                r.volumeSteadyState = notComputed(
                    "an intravenous reference: after an extravascular dose MRT includes "
                    "mean absorption time, so CL*MRT is not Vss");
        } else {
            const char* missing = "a positive dose";
            r.clearance = notComputed(missing);
            r.volumeZ = notComputed(missing);
            r.volumeSteadyState = notComputed(missing);
        }

        if (r.extrapolationUnreliable) {
            const std::string pct = std::to_string(pctExtrap);
            r.warnings.push_back("percent of AUCinf extrapolated beyond the last sample is " +
                                 pct + "%, above the conventional 20% limit: AUCinf is mostly "
                                       "the terminal fit rather than the data");
            for (const char* derived :
                 {"AUCinf", "CL (or CL/F)", "Vz (or Vz/F)", "Vss", "MRT", "AUMC"})
                r.warnings.push_back(std::string(derived) +
                                     " inherits that extrapolation and is unreliable");
        }
        if (best.points < 3)
            r.warnings.push_back("lambda_z was fitted to fewer than three points");
    }

    r.warnings.push_back("noncompartmental analysis describes the exposure that was "
                         "observed; it is not a dose, a dose adjustment or a prediction "
                         "about anyone who was not sampled");
    return r;
}

}  // namespace biocad::sim
