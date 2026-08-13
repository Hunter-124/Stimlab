#include "numeric/Optimize.h"

#include <algorithm>
#include <cmath>
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

    // Standard errors only when the normal-equation matrix is genuinely invertible
    // and there are residual degrees of freedom. A fabricated error bar is worse
    // than no error bar, so the failure mode is an empty vector plus a reason.
    if (observations <= nParams) {
        note += "; standard errors unavailable: no residual degrees of freedom";
    } else {
        Eigen::LDLT<Eigen::MatrixXd> solver(jtwj);
        const double dof = static_cast<double>(observations - nParams);
        bool ok = solver.info() == Eigen::Success && jtwj.allFinite();
        Eigen::MatrixXd cov;
        if (ok) {
            cov = solver.solve(Eigen::MatrixXd::Identity(static_cast<Eigen::Index>(nParams),
                                                         static_cast<Eigen::Index>(nParams)));
            ok = solver.info() == Eigen::Success && cov.allFinite();
            if (ok) {
                const Eigen::MatrixXd check = jtwj * cov;
                const double err = (check - Eigen::MatrixXd::Identity(
                                                static_cast<Eigen::Index>(nParams),
                                                static_cast<Eigen::Index>(nParams)))
                                       .cwiseAbs()
                                       .maxCoeff();
                ok = err < 1e-6;
            }
        }
        if (!ok) {
            note += "; standard errors unavailable: the normal-equation matrix is rank-deficient";
        } else {
            const double variance = cost * 2.0 / dof;
            std::vector<double> errors(nParams, 0.0);
            for (std::size_t j = 0; j < nParams; ++j) {
                const double v = cov(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(j))
                                 * variance;
                if (!(v >= 0.0) || !std::isfinite(v)) {
                    errors.clear();
                    note += "; standard errors unavailable: non-positive covariance diagonal";
                    break;
                }
                errors[j] = std::sqrt(v);
            }
            result.standardErrors = errors;
        }
    }

    result.note = note;
    return result;
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
