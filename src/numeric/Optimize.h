// numeric/Optimize.h - Levenberg-Marquardt least squares with an analytic Jacobian.
//
// One fitter for the whole application. Every curve fit (dose-response, Schild,
// binding kinetics) goes through this rather than growing a second implementation
// with subtly different convergence behaviour.
//
// The caller supplies residuals AND their analytic Jacobian: finite differences on
// a four-parameter logistic in log space lose several digits near the asymptotes,
// and the analytic derivatives are short closed forms anyway.
#pragma once

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace biocad::numeric {

struct LmOptions {
    int    maxIterations = 200;
    double lambdaInit = 1e-3;   // damping; grows on a rejected step, shrinks on accept
    double lambdaUp = 10.0;
    double lambdaDown = 10.0;
    double tolCost = 1e-14;     // relative cost improvement to declare convergence
    double tolStep = 1e-14;     // relative parameter step to declare convergence
};

struct LmResult {
    std::vector<double> params;
    std::vector<double> standardErrors;  // from the covariance matrix; empty if rank-deficient
    double cost = 0.0;                   // 0.5 * sum(residual^2) at the returned params
    double rSquared = 0.0;
    int    iterations = 0;
    bool   converged = false;
    std::string note;                    // why it stopped, or why errors are unavailable
};

// residuals[i] = model(params, i) - observed[i]
// jacobian is row-major, nObs x nParams: jacobian[i * nParams + j] = d residual[i] / d params[j]
using LmEvaluate = std::function<void(const std::vector<double>& params,
                                      std::vector<double>& residuals,
                                      std::vector<double>& jacobian)>;

// `weights` may be empty (unweighted) or hold one non-negative weight per observation;
// the fit then minimises sum(w_i * r_i^2), which is how 1/y^2 weighting is applied.
LmResult levenbergMarquardt(std::vector<double> initial, std::size_t observations,
                            const LmEvaluate& evaluate,
                            const std::vector<double>& weights = {},
                            const LmOptions& options = {});

// Coefficient of determination of a fitted series against observations. Returns 0
// when the observations have no variance, rather than dividing by zero.
double rSquared(const std::vector<double>& observed, const std::vector<double>& fitted);

}  // namespace biocad::numeric
