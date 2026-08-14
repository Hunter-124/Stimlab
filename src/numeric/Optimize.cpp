#include "numeric/Optimize.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <Eigen/Dense>

namespace biocad::numeric {
namespace {

// Weighted cost 0.5 * sum(w_i * r_i^2). Returns false when any residual is
// non-finite: an unusable model evaluation must stop the fit, not poison it with
// NaN that silently propagates into a "converged" answer.
bool weightedCost(const std::vector<double>& residuals, const std::vector<double>& weights,
                  double& out) {
    double sum = 0.0;
    for (std::size_t i = 0; i < residuals.size(); ++i) {
        const double r = residuals[i];
        if (!std::isfinite(r)) {
            return false;
        }
        const double w = weights.empty() ? 1.0 : weights[i];
        sum += w * r * r;
    }
    out = 0.5 * sum;
    return std::isfinite(out);
}

}  // namespace

LmResult levenbergMarquardt(std::vector<double> initial, std::size_t observations,
                            const LmEvaluate& evaluate, const std::vector<double>& weights,
                            const LmOptions& options) {
    LmResult result;
    result.params = initial;

    const std::size_t nParams = initial.size();
    if (nParams == 0) {
        result.note = "no parameters to fit";
        return result;
    }
    if (observations == 0) {
        result.note = "no observations supplied";
        return result;
    }
    if (observations < nParams) {
        result.note = "fewer observations (" + std::to_string(observations) + ") than parameters ("
                      + std::to_string(nParams) + "); the fit is underdetermined";
        return result;
    }
    if (!weights.empty() && weights.size() != observations) {
        result.note = "weight count does not match the observation count";
        return result;
    }

    std::vector<double> residuals(observations, 0.0);
    std::vector<double> jacobian(observations * nParams, 0.0);

    evaluate(result.params, residuals, jacobian);
    double cost = 0.0;
    if (!weightedCost(residuals, weights, cost)) {
        result.note = "the model produced a non-finite residual at the initial parameters";
        return result;
    }
    result.cost = cost;

    Eigen::MatrixXd jtwj(static_cast<Eigen::Index>(nParams), static_cast<Eigen::Index>(nParams));
    Eigen::VectorXd jtwr(static_cast<Eigen::Index>(nParams));
    double lambda = options.lambdaInit;
    bool converged = false;
    std::string note;

    for (int iter = 0; iter < options.maxIterations; ++iter) {
        result.iterations = iter + 1;

        jtwj.setZero();
        jtwr.setZero();
        for (std::size_t i = 0; i < observations; ++i) {
            const double w = weights.empty() ? 1.0 : weights[i];
            const double* row = &jacobian[i * nParams];
            for (std::size_t j = 0; j < nParams; ++j) {
                const double wj = w * row[j];
                jtwr(static_cast<Eigen::Index>(j)) += wj * residuals[i];
                for (std::size_t k = j; k < nParams; ++k) {
                    jtwj(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(k)) += wj * row[k];
                }
            }
        }
        for (std::size_t j = 0; j < nParams; ++j) {
            for (std::size_t k = 0; k < j; ++k) {
                jtwj(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(k)) =
                    jtwj(static_cast<Eigen::Index>(k), static_cast<Eigen::Index>(j));
            }
        }

        // Marquardt scaling: damp along the curvature of each parameter rather than
        // an isotropic identity, so parameters on wildly different scales (log EC50
        // against a Top of 100) are damped comparably.
        Eigen::VectorXd diag = jtwj.diagonal();
        for (Eigen::Index j = 0; j < diag.size(); ++j) {
            if (!(diag(j) > 0.0)) {
                diag(j) = 1.0;
            }
        }

        bool accepted = false;
        for (int attempt = 0; attempt < 60 && !accepted; ++attempt) {
            Eigen::MatrixXd damped = jtwj;
            damped.diagonal() += lambda * diag;

            Eigen::LDLT<Eigen::MatrixXd> solver(damped);
            if (solver.info() != Eigen::Success) {
                lambda *= options.lambdaUp;
                continue;
            }
            const Eigen::VectorXd step = solver.solve(-jtwr);
            if (!step.allFinite()) {
                lambda *= options.lambdaUp;
                continue;
            }

            std::vector<double> trial(nParams);
            for (std::size_t j = 0; j < nParams; ++j) {
                trial[j] = result.params[j] + step(static_cast<Eigen::Index>(j));
            }

            std::vector<double> trialResiduals(observations, 0.0);
            std::vector<double> trialJacobian(observations * nParams, 0.0);
            evaluate(trial, trialResiduals, trialJacobian);
            double trialCost = 0.0;
            if (!weightedCost(trialResiduals, weights, trialCost) || trialCost >= cost) {
                lambda *= options.lambdaUp;
                continue;
            }

            double stepNorm = 0.0;
            double paramNorm = 0.0;
            for (std::size_t j = 0; j < nParams; ++j) {
                stepNorm += step(static_cast<Eigen::Index>(j)) * step(static_cast<Eigen::Index>(j));
                paramNorm += result.params[j] * result.params[j];
            }
            stepNorm = std::sqrt(stepNorm);
            paramNorm = std::sqrt(paramNorm);

            const double improvement = (cost - trialCost) / std::max(cost, 1e-300);
            result.params = trial;
            residuals = trialResiduals;
            jacobian = trialJacobian;
            cost = trialCost;
            lambda = std::max(lambda / options.lambdaDown, 1e-300);
            accepted = true;

            if (improvement <= options.tolCost) {
                converged = true;
                note = "relative cost improvement below tolCost";
            } else if (stepNorm <= options.tolStep * (paramNorm + options.tolStep)) {
                converged = true;
                note = "relative parameter step below tolStep";
            }
        }

        if (!accepted) {
            // No downhill step exists at any damping: this is a local minimum in
            // every direction the model can move, which is convergence.
            converged = true;
            note = "no cost-reducing step available at any damping";
        }
        if (converged) {
            break;
        }
    }

    result.cost = cost;
    result.converged = converged;
    if (!converged) {
        note = "hit maxIterations (" + std::to_string(options.maxIterations)
               + ") without meeting a convergence tolerance";
    }

    // Standard errors only when the weighted Jacobian is genuinely full rank and
    // there are residual degrees of freedom. The gate is an SVD, not an inversion
    // attempt: a nearly-singular normal-equation matrix can be "successfully"
    // factorized and still hand back error bars that are pure noise amplification,
    // and the condition number is the number that says so. A fabricated error bar
    // is worse than no error bar, so the failure mode is an empty vector plus a
    // reason.
    {
        std::vector<double> finalResiduals(observations, 0.0);
        std::vector<double> finalJacobian(observations * nParams, 0.0);
        evaluate(result.params, finalResiduals, finalJacobian);

        Eigen::MatrixXd jw(static_cast<Eigen::Index>(observations),
                           static_cast<Eigen::Index>(nParams));
        for (std::size_t i = 0; i < observations; ++i) {
            const double sw = weights.empty() ? 1.0 : std::sqrt(std::max(0.0, weights[i]));
            for (std::size_t j = 0; j < nParams; ++j) {
                jw(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
                    sw * finalJacobian[i * nParams + j];
            }
        }

        if (!jw.allFinite()) {
            note += "; standard errors unavailable: the Jacobian is not finite";
        } else {
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(jw, Eigen::ComputeThinU | Eigen::ComputeThinV);
            const Eigen::VectorXd& s = svd.singularValues();
            const double sMax = s.size() > 0 ? s(0) : 0.0;
            const double sMin = s.size() > 0 ? s(s.size() - 1) : 0.0;
            // The standard numerical-rank tolerance: max(m,n) * eps * largest
            // singular value. Anything below it is indistinguishable from zero at
            // double precision, so a direction the data does not constrain cannot
            // masquerade as a constrained one.
            const double tol = static_cast<double>(std::max(observations, nParams))
                               * std::numeric_limits<double>::epsilon() * sMax;
            std::size_t rank = 0;
            for (Eigen::Index k = 0; k < s.size(); ++k) {
                if (s(k) > tol) {
                    ++rank;
                }
            }
            result.rank = rank;
            result.conditionNumber = (sMin > 0.0) ? sMax / sMin
                                                  : std::numeric_limits<double>::infinity();

            if (observations <= nParams) {
                note += "; standard errors unavailable: no residual degrees of freedom";
            } else if (rank < nParams) {
                note += "; standard errors unavailable: the weighted Jacobian is rank-deficient (rank "
                        + std::to_string(rank) + " of " + std::to_string(nParams)
                        + "; the data do not identify every parameter)";
            } else {
                // cov = (J'WJ)^-1 * variance = V S^-2 V' * variance.
                Eigen::MatrixXd sInv2 = Eigen::MatrixXd::Zero(
                    static_cast<Eigen::Index>(nParams), static_cast<Eigen::Index>(nParams));
                for (Eigen::Index k = 0; k < s.size(); ++k) {
                    sInv2(k, k) = 1.0 / (s(k) * s(k));
                }
                const double dof = static_cast<double>(observations - nParams);
                const double variance = cost * 2.0 / dof;
                const Eigen::MatrixXd cov = svd.matrixV() * sInv2
                                            * svd.matrixV().transpose() * variance;
                std::vector<double> errors(nParams, 0.0);
                std::vector<double> flat(nParams * nParams, 0.0);
                bool ok = cov.allFinite();
                for (std::size_t j = 0; ok && j < nParams; ++j) {
                    const double v = cov(static_cast<Eigen::Index>(j),
                                         static_cast<Eigen::Index>(j));
                    if (!(v >= 0.0) || !std::isfinite(v)) {
                        ok = false;
                        break;
                    }
                    errors[j] = std::sqrt(v);
                    for (std::size_t k = 0; k < nParams; ++k) {
                        flat[j * nParams + k] = cov(static_cast<Eigen::Index>(j),
                                                    static_cast<Eigen::Index>(k));
                    }
                }
                if (!ok) {
                    note += "; standard errors unavailable: non-positive covariance diagonal";
                } else {
                    result.standardErrors = errors;
                    result.covariance = flat;
                }
            }
        }
    }

    result.note = note;
    return result;
}

IrlsResult tukeyBiweight(std::vector<double> initial, std::size_t observations,
                         const LmEvaluate& evaluate, const std::vector<double>& priorWeights,
                         const IrlsOptions& irls, const LmOptions& options) {
    IrlsResult out;
    if (!priorWeights.empty() && priorWeights.size() != observations) {
        out.note = "priorWeights size does not match the observation count";
        return out;
    }
    if (irls.tuning <= 0.0) {
        out.note = "the biweight tuning constant must be positive";
        return out;
    }

    std::vector<double> params = std::move(initial);
    std::vector<double> weights = priorWeights.empty()
                                      ? std::vector<double>(observations, 1.0)
                                      : priorWeights;
    std::vector<double> residuals(observations, 0.0);
    std::vector<double> jacobian(observations * params.size(), 0.0);
    double previousScale = 0.0;

    for (int iter = 0; iter < irls.maxIterations; ++iter) {
        out.fit = levenbergMarquardt(params, observations, evaluate, weights, options);
        params = out.fit.params;
        out.iterations = iter + 1;

        evaluate(params, residuals, jacobian);

        // Scale from the median absolute deviation of the residuals about zero:
        // the regression residuals are already centred by the fit, so re-centring
        // them on their own median would hide a systematic offset.
        std::vector<double> absResiduals(observations, 0.0);
        for (std::size_t i = 0; i < observations; ++i) {
            absResiduals[i] = std::abs(residuals[i]);
        }
        std::vector<double> sorted = absResiduals;
        std::sort(sorted.begin(), sorted.end());
        const std::size_t mid = sorted.size() / 2;
        const double mad = sorted.empty() ? 0.0
                           : (sorted.size() % 2 == 1 ? sorted[mid]
                                                     : 0.5 * (sorted[mid - 1] + sorted[mid]));
        const double scale = 1.4826 * mad;
        out.scale = scale;

        if (!(scale > 0.0)) {
            // At least half the residuals are exactly zero, so the MAD is zero and
            // there is no scale left to downweight against - which is the normal end
            // state on clean data once an outlier has already been rejected. The
            // weights that produced THIS fit are the answer; overwriting them with
            // ones would report "nothing was downweighted" about a fit that only
            // looks exact because something was.
            if (out.robustWeights.empty()) {
                out.robustWeights = std::vector<double>(observations, 1.0);
            }
            out.converged = true;
            out.note = "residual scale reached zero: the remaining observations are fit exactly "
                       "under the weights reported here";
            return out;
        }

        std::vector<double> robust(observations, 0.0);
        for (std::size_t i = 0; i < observations; ++i) {
            const double u = residuals[i] / (irls.tuning * scale);
            const double t = 1.0 - u * u;
            robust[i] = (std::abs(u) < 1.0) ? t * t : 0.0;   // Tukey biweight
            weights[i] = robust[i] * (priorWeights.empty() ? 1.0 : priorWeights[i]);
        }
        out.robustWeights = robust;

        if (iter > 0 && previousScale > 0.0
            && std::abs(scale - previousScale) <= irls.tolScale * previousScale) {
            out.converged = true;
            return out;
        }
        previousScale = scale;
    }

    out.note = "hit IrlsOptions::maxIterations (" + std::to_string(irls.maxIterations)
               + ") without a stable residual scale";
    return out;
}

namespace {

// The chi-square(1 df) quantiles for the confidence levels this function accepts.
// A table, not an interpolation: an interval whose threshold was interpolated
// between tabulated quantiles is not the interval it claims to be.
bool chiSquareOneDf(double confidence, double& out) {
    if (std::abs(confidence - 0.68) < 1e-9) { out = 0.9891; return true; }
    if (std::abs(confidence - 0.90) < 1e-9) { out = 2.7055; return true; }
    if (std::abs(confidence - 0.95) < 1e-9) { out = 3.8415; return true; }
    if (std::abs(confidence - 0.99) < 1e-9) { out = 6.6349; return true; }
    return false;
}

// Weighted sum of squared residuals - the quantity the profile threshold is set on.
bool sumOfSquares(const std::vector<double>& params, std::size_t observations,
                  const LmEvaluate& evaluate, const std::vector<double>& weights, double& out) {
    std::vector<double> residuals(observations, 0.0);
    std::vector<double> jacobian(observations * params.size(), 0.0);
    evaluate(params, residuals, jacobian);
    double ssr = 0.0;
    for (std::size_t i = 0; i < observations; ++i) {
        if (!std::isfinite(residuals[i])) {
            return false;
        }
        const double w = weights.empty() ? 1.0 : weights[i];
        ssr += w * residuals[i] * residuals[i];
    }
    out = ssr;
    return true;
}

}  // namespace

ProfileInterval profileLikelihood(const std::vector<double>& best, std::size_t observations,
                                  const LmEvaluate& evaluate, std::size_t index,
                                  double confidence, const std::vector<double>& weights,
                                  const LmOptions& options) {
    ProfileInterval out;
    out.confidence = confidence;
    const std::size_t nParams = best.size();
    if (index >= nParams) {
        out.note = "parameter index out of range";
        return out;
    }
    if (observations <= nParams) {
        out.note = "no residual degrees of freedom: a profile interval is undefined";
        return out;
    }
    double quantile = 0.0;
    if (!chiSquareOneDf(confidence, quantile)) {
        out.note = "unsupported confidence level: use 0.68, 0.90, 0.95 or 0.99";
        return out;
    }

    double ssrBest = 0.0;
    if (!sumOfSquares(best, observations, evaluate, weights, ssrBest)) {
        out.note = "the model does not evaluate finitely at the reported optimum";
        return out;
    }
    const double dof = static_cast<double>(observations - nParams);
    const double variance = ssrBest / dof;
    out.threshold = ssrBest + quantile * variance;
    out.lower = best[index];
    out.upper = best[index];
    if (!(variance > 0.0)) {
        out.note = "the fit is exact: the profile is flat and the interval is the point estimate";
        out.lowerFound = out.upperFound = true;
        return out;
    }

    // Profiling parameter `index` means refitting the OTHERS with it held fixed,
    // so build a reduced evaluate that inserts the fixed value and drops the
    // corresponding Jacobian column.
    const auto profileAt = [&](double value) -> double {
        std::vector<double> reduced;
        reduced.reserve(nParams - 1);
        for (std::size_t j = 0; j < nParams; ++j) {
            if (j != index) {
                reduced.push_back(best[j]);
            }
        }
        if (reduced.empty()) {
            std::vector<double> full = best;
            full[index] = value;
            double ssr = 0.0;
            return sumOfSquares(full, observations, evaluate, weights, ssr)
                       ? ssr : std::numeric_limits<double>::infinity();
        }
        const std::size_t nFree = nParams - 1;
        LmEvaluate wrapped = [&, value](const std::vector<double>& p,
                                        std::vector<double>& residuals,
                                        std::vector<double>& jac) {
            std::vector<double> full(nParams, 0.0);
            for (std::size_t j = 0, k = 0; j < nParams; ++j) {
                full[j] = (j == index) ? value : p[k++];
            }
            std::vector<double> fullJac(observations * nParams, 0.0);
            evaluate(full, residuals, fullJac);
            for (std::size_t i = 0; i < observations; ++i) {
                for (std::size_t j = 0, k = 0; j < nParams; ++j) {
                    if (j != index) {
                        jac[i * nFree + k++] = fullJac[i * nParams + j];
                    }
                }
            }
        };
        const LmResult refit = levenbergMarquardt(reduced, observations, wrapped, weights,
                                                  options);
        std::vector<double> full(nParams, 0.0);
        for (std::size_t j = 0, k = 0; j < nParams; ++j) {
            full[j] = (j == index) ? value : refit.params[k++];
        }
        double ssr = 0.0;
        return sumOfSquares(full, observations, evaluate, weights, ssr)
                   ? ssr : std::numeric_limits<double>::infinity();
    };

    // Step out geometrically until the profile crosses the threshold, then bisect.
    // The initial step is 10% of the estimate (or 1e-3 absolute for a parameter
    // sitting at zero), doubled up to 40 times, which spans nine orders of
    // magnitude before giving up and saying it never crossed.
    const double centre = best[index];
    const double step0 = (std::abs(centre) > 0.0) ? 0.1 * std::abs(centre) : 1e-3;
    for (int direction = -1; direction <= 1; direction += 2) {
        double step = step0;
        double outside = 0.0;
        bool crossed = false;
        for (int k = 0; k < 40; ++k) {
            const double trial = centre + direction * step;
            if (profileAt(trial) > out.threshold) {
                outside = trial;
                crossed = true;
                break;
            }
            step *= 2.0;
        }
        if (!crossed) {
            continue;
        }
        double inside = centre;
        for (int k = 0; k < 200; ++k) {
            const double mid = 0.5 * (inside + outside);
            if (profileAt(mid) > out.threshold) {
                outside = mid;
            } else {
                inside = mid;
            }
            if (std::abs(outside - inside) <= 1e-12 * std::max(1.0, std::abs(mid))) {
                break;
            }
        }
        if (direction < 0) {
            out.lower = 0.5 * (inside + outside);
            out.lowerFound = true;
        } else {
            out.upper = 0.5 * (inside + outside);
            out.upperFound = true;
        }
    }

    if (!out.lowerFound || !out.upperFound) {
        out.note = "the profile did not cross the threshold on at least one side within the "
                   "search span: that side of the interval is unbounded by these data";
    }
    return out;
}

double aicc(double sumOfSquaredResiduals, std::size_t observations, std::size_t parameters) {
    const double n = static_cast<double>(observations);
    const double k = static_cast<double>(parameters);
    if (observations == 0 || !(sumOfSquaredResiduals > 0.0) || n <= k + 1.0) {
        return std::numeric_limits<double>::infinity();
    }
    return n * std::log(sumOfSquaredResiduals / n) + 2.0 * k + 2.0 * k * (k + 1.0) / (n - k - 1.0);
}

double rSquared(const std::vector<double>& observed, const std::vector<double>& fitted) {
    if (observed.empty() || observed.size() != fitted.size()) {
        return 0.0;
    }
    double mean = 0.0;
    for (double v : observed) {
        mean += v;
    }
    mean /= static_cast<double>(observed.size());

    double ssRes = 0.0;
    double ssTot = 0.0;
    for (std::size_t i = 0; i < observed.size(); ++i) {
        const double d = observed[i] - fitted[i];
        const double t = observed[i] - mean;
        ssRes += d * d;
        ssTot += t * t;
    }
    if (ssTot == 0.0) {
        return 0.0;
    }
    return 1.0 - ssRes / ssTot;
}

}  // namespace biocad::numeric
