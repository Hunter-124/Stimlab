#include "sim/Population.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "numeric/Ode.h"
#include "pkpd/PkEngine.h"
#include "sim/Random.h"

namespace biocad::sim {
namespace {

// Which PkModelSpec field a named random effect multiplies. Parameters are named by
// the caller as strings because Omega comes from a fit or a pack, not from C++, and
// an unrecognised name has to be REPORTED rather than silently ignored - a warning
// that says "CLx was not applied" is the difference between a band that is wrong and
// a band that is known to be wrong.
Quantity* parameterSlot(PkModelSpec& spec, const std::string& name) {
    if (name == "CL") return &spec.clearance;
    if (name == "V" || name == "V1") return &spec.volume;
    if (name == "V2") return &spec.volumePeripheral;
    if (name == "Q") return &spec.intercompartmental;
    if (name == "ka") return &spec.absorptionRate;
    if (name == "F") return &spec.bioavailability;
    if (name == "fu") return &spec.unboundFraction;
    if (name == "Vmax") return &spec.vmax;
    if (name == "Km") return &spec.km;
    return nullptr;
}

// Linear-interpolated percentile (the "type 7" definition, which is what R's
// quantile() and NumPy's percentile() both default to) over an already-sorted
// sample. Stated explicitly because a percentile is not one number: nearest-rank
// and interpolated definitions differ by a whole order statistic at small N.
double percentileSorted(const std::vector<double>& sorted, double q) {
    if (sorted.empty()) return 0.0;
    if (sorted.size() == 1) return sorted.front();
    const double pos = q * static_cast<double>(sorted.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
    const std::size_t hi = std::min(lo + 1, sorted.size() - 1);
    const double frac = pos - static_cast<double>(lo);
    return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
}

}  // namespace

PopulationProfile simulatePopulation(const PkModelSpec& model, const DoseRegimen& regimen,
                                     const VariabilitySpec& variability) {
    PopulationProfile out;
    out.spec = variability;
    out.provenanceStatement =
        "This band is the variability that was entered - Omega, the fit covariance and the "
        "residual error above - propagated through the stated PK model with seed " +
        std::to_string(variability.seed) +
        "; it is not a prediction about any individual, and it is not a dose.";

    const int nSubjects = std::max(1, variability.subjects);
    const std::size_t nParams = variability.parameters.size();

    // The typical-value profile also fixes the time grid: every subject shares it,
    // which is what makes a percentile across subjects at a given time meaningful.
    const PkProfile typical = pkpd::simulate(model, regimen);
    out.times = typical.timeH;
    const std::size_t nTimes = out.times.size();
    if (nTimes == 0) {
        out.medianAuc = notComputed("a non-empty simulated time grid");
        out.medianCmax = notComputed("a non-empty simulated time grid");
        out.aucCv = notComputed("a non-empty simulated time grid");
        out.warnings.push_back("the PK engine returned an empty profile");
        return out;
    }

    // ONE allocation, indexed [timeIndex * nSubjects + subject]: time is the slow
    // axis, so all subjects at a given time are contiguous and the percentile pass
    // is a straight sort over a contiguous span instead of a strided gather. At
    // 1000 subjects x 4000 time points this is a single 32 MB block; the same data
    // as a vector-of-vectors would be 1000 separate allocations, would touch four
    // cache lines per percentile element, and could not be reduced in one sweep.
    std::vector<double> conc(nTimes * static_cast<std::size_t>(nSubjects), 0.0);

    Pcg64Dxsm rng(variability.seed);

    // Layer 1: between-subject variability. Omega is the covariance of the LOG-scale
    // random effects, so Pi = theta_i * exp(eta_i) keeps every parameter positive -
    // an additive eta would eventually hand the integrator a negative clearance.
    std::vector<double> omegaChol;
    bool bsv = variability.betweenSubject && nParams > 0;
    if (bsv) {
        if (variability.omega.size() != nParams * nParams ||
            !choleskyLower(variability.omega, nParams, omegaChol)) {
            bsv = false;
            out.warnings.push_back(
                "Omega is not a positive-definite covariance over the named parameters, so "
                "between-subject variability was NOT applied; the matrix is reported back "
                "unchanged rather than repaired, because a nudged Omega is no longer the "
                "variability that was entered");
        }
    }

    // Layer 2: parameter uncertainty. Latin hypercube over the unit hypercube, rank-
    // correlated by Iman-Conover to the correlation implied by the fit covariance,
    // then mapped to log-normal multipliers. LHS rather than plain Monte Carlo
    // because the uncertainty layer is usually the smaller sample and stratification
    // is what stops a 200-draw uncertainty band from being lumpy.
    std::vector<double> lhs;
    std::vector<double> uncertaintySd(nParams, 0.0);
    bool uncertainty = variability.parameterUncertainty && nParams > 0 &&
                       variability.parameterCovariance.size() == nParams * nParams;
    if (variability.parameterUncertainty && !uncertainty)
        out.warnings.push_back("parameter uncertainty was requested but the fit covariance is "
                               "not square over the named parameters, so it was NOT applied");
    if (uncertainty) {
        std::vector<double> corr(nParams * nParams, 0.0);
        for (std::size_t i = 0; i < nParams; ++i)
            uncertaintySd[i] = std::sqrt(std::max(
                variability.parameterCovariance[i * nParams + i], 0.0));
        for (std::size_t i = 0; i < nParams; ++i) {
            for (std::size_t j = 0; j < nParams; ++j) {
                const double d = uncertaintySd[i] * uncertaintySd[j];
                corr[i * nParams + j] =
                    (i == j) ? 1.0
                             : (d > 0 ? variability.parameterCovariance[i * nParams + j] / d
                                      : 0.0);
            }
        }
        lhs = latinHypercube(static_cast<std::size_t>(nSubjects), nParams, rng);
        if (nParams >= 2 && nSubjects >= 3 && !imanConover(lhs, static_cast<std::size_t>(nSubjects),
                                                           nParams, corr))
            out.warnings.push_back("the correlation implied by the fit covariance is not "
                                   "positive definite, so the uncertainty draws are "
                                   "stratified but UNCORRELATED");
    }

    const bool residual = variability.residualError &&
                          (variability.proportionalResidualCv > 0 ||
                           variability.additiveResidualSd > 0);

    std::vector<double> subjectAuc(static_cast<std::size_t>(nSubjects), 0.0);
    std::vector<double> subjectCmax(static_cast<std::size_t>(nSubjects), 0.0);
    std::vector<std::string> unknownNames;
    std::vector<double> eta;

    for (int s = 0; s < nSubjects; ++s) {
        PkModelSpec subject = model;

        if (bsv) {
            multivariateNormal(omegaChol, nParams, rng, eta);
            for (std::size_t p = 0; p < nParams; ++p) {
                Quantity* slot = parameterSlot(subject, variability.parameters[p]);
                if (!slot) {
                    if (s == 0) unknownNames.push_back(variability.parameters[p]);
                    continue;
                }
                slot->value *= std::exp(eta[p]);
            }
        }
        if (uncertainty) {
            for (std::size_t p = 0; p < nParams; ++p) {
                Quantity* slot = parameterSlot(subject, variability.parameters[p]);
                if (!slot) continue;
                const double u = lhs[static_cast<std::size_t>(s) * nParams + p];
                slot->value *= std::exp(inverseNormalCdf(u) * uncertaintySd[p]);
            }
        }

        const PkProfile profile = pkpd::simulate(subject, regimen);
        const std::size_t m = std::min(nTimes, profile.concentrationMgPerL.size());
        double cmax = 0.0;
        for (std::size_t t = 0; t < m; ++t) {
            double c = profile.concentrationMgPerL[t];
            if (residual) {
                // Proportional then additive, and the result is CLAMPED at zero: a
                // negative observed concentration is not a measurement, and letting
                // one into the band drags the 5th percentile below the axis.
                if (variability.proportionalResidualCv > 0)
                    c *= 1.0 + variability.proportionalResidualCv * rng.nextNormal();
                if (variability.additiveResidualSd > 0)
                    c += variability.additiveResidualSd * rng.nextNormal();
                c = std::max(c, 0.0);
            }
            conc[t * static_cast<std::size_t>(nSubjects) + static_cast<std::size_t>(s)] = c;
            cmax = std::max(cmax, c);
        }
        subjectCmax[static_cast<std::size_t>(s)] = cmax;

        std::vector<double> series(nTimes);
        for (std::size_t t = 0; t < nTimes; ++t)
            series[t] = conc[t * static_cast<std::size_t>(nSubjects) + static_cast<std::size_t>(s)];
        subjectAuc[static_cast<std::size_t>(s)] = numeric::trapezoid(out.times, series);

        if (static_cast<int>(out.sampleTrajectories.size()) < kMaxStoredTrajectories)
            out.sampleTrajectories.push_back(std::move(series));
    }

    for (const auto& n : unknownNames)
        out.warnings.push_back("random effect '" + n +
                               "' does not name a PK parameter (CL, V/V1, V2, Q, ka, F, fu, "
                               "Vmax, Km) and was NOT applied");

    out.bands.resize(nTimes);
    std::vector<double> column(static_cast<std::size_t>(nSubjects));
    for (std::size_t t = 0; t < nTimes; ++t) {
        const double* row = &conc[t * static_cast<std::size_t>(nSubjects)];
        column.assign(row, row + nSubjects);
        std::sort(column.begin(), column.end());
        out.bands[t].timeH = out.times[t];
        out.bands[t].p5 = percentileSorted(column, 0.05);
        out.bands[t].p50 = percentileSorted(column, 0.50);
        out.bands[t].p95 = percentileSorted(column, 0.95);
    }

    std::vector<double> sortedAuc = subjectAuc, sortedCmax = subjectCmax;
    std::sort(sortedAuc.begin(), sortedAuc.end());
    std::sort(sortedCmax.begin(), sortedCmax.end());
    // Model, not Predicted: this is a constructed artefact from entered parameters,
    // and no benchmark error exists for it because no benchmark was run.
    out.medianAuc = makeQuantity(percentileSorted(sortedAuc, 0.5), "mg*h/L", 0,
                                 Provenance::Model,
                                 "median over " + std::to_string(nSubjects) +
                                     " simulated subjects, trapezoid over the horizon");
    out.medianCmax = makeQuantity(percentileSorted(sortedCmax, 0.5), "mg/L", 0,
                                  Provenance::Model,
                                  "median over " + std::to_string(nSubjects) +
                                      " simulated subjects");

    double mean = 0;
    for (double a : subjectAuc) mean += a;
    mean /= static_cast<double>(nSubjects);
    double var = 0;
    for (double a : subjectAuc) var += (a - mean) * (a - mean);
    if (nSubjects > 1 && mean > 0) {
        var /= static_cast<double>(nSubjects - 1);
        out.aucCv = makeQuantity(100.0 * std::sqrt(var) / mean, "%", 0, Provenance::Model,
                                 "coefficient of variation of AUC across simulated subjects");
    } else {
        out.aucCv = notComputed("more than one subject with a positive mean AUC");
    }

    out.assumptions.push_back(
        "eta ~ MVN(0, Omega) on the log scale, so Pi = theta_i * exp(eta_i) and every "
        "parameter stays positive");
    out.assumptions.push_back(
        std::string("between-subject variability: ") + (bsv ? "ON" : "OFF") +
        "; parameter uncertainty: " + (uncertainty ? "ON (Latin hypercube + Iman-Conover)" : "OFF") +
        "; residual error: " + (residual ? "ON" : "OFF"));
    out.assumptions.push_back("sampler: PCG64-DXSM with Wichura AS241 normals, seed " +
                              std::to_string(variability.seed) +
                              " - the same seed reproduces this band exactly");
    out.assumptions.push_back(std::to_string(nSubjects) + " subjects; at most " +
                              std::to_string(kMaxStoredTrajectories) +
                              " individual trajectories are stored for the overlay, while the "
                              "percentiles use every subject");
    for (const auto& a : typical.assumptions) out.assumptions.push_back(a);
    if (typical.flipFlop)
        out.warnings.push_back("ka < ke in the typical parameters: the terminal phase reflects "
                               "absorption, not elimination");
    return out;
}

}  // namespace biocad::sim
