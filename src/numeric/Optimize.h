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
    std::vector<double> covariance;      // row-major nParams x nParams; empty with the errors
    double cost = 0.0;                   // 0.5 * sum(residual^2) at the returned params
    double rSquared = 0.0;
    int    iterations = 0;
    bool   converged = false;
    // Numerical rank of the weighted Jacobian and its condition number, from the
    // SVD that gates the standard errors. A caller that wants to know whether the
    // experiment identified the parameters reads these rather than guessing from
    // an empty error vector.
    std::size_t rank = 0;
    double      conditionNumber = 0.0;
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

// Tukey-biweight iteratively reweighted least squares. Robust regression is the
// honest answer to one bad well or one pipetting error: it downweights an outlier
// instead of deleting it, so the point stays on the plot and in the record.
//
// The scale is re-estimated each iteration as 1.4826 * MAD of the current
// residuals (1.4826 makes the MAD a consistent estimator of sigma for Gaussian
// data). `tuning` = 4.685 gives 95% asymptotic efficiency at the normal model,
// which is the published default and is stated here so it is not mistaken for a
// tuned-to-taste constant.
struct IrlsOptions {
    int    maxIterations = 30;
    double tuning = 4.685;
    double tolScale = 1e-12;   // relative change in scale to declare convergence
};

struct IrlsResult {
    LmResult            fit;            // the final weighted fit
    std::vector<double> robustWeights;  // 0..1 per observation; 0 = fully rejected
    double              scale = 0.0;    // 1.4826 * MAD of the final residuals
    int                 iterations = 0;
    bool                converged = false;
    std::string         note;
};

// `priorWeights` multiply the robust weights, so replicate weighting and
// robustness compose rather than competing.
IrlsResult tukeyBiweight(std::vector<double> initial, std::size_t observations,
                         const LmEvaluate& evaluate,
                         const std::vector<double>& priorWeights = {},
                         const IrlsOptions& irls = {},
                         const LmOptions& options = {});

// A profile-likelihood confidence interval for one parameter.
//
// WHY not params +/- 1.96 * SE: the Wald interval assumes the cost surface is a
// paraboloid in the parameter, which an EC50 or a kd is not. The profile interval
// refits every other parameter at each trial value and reports where the sum of
// squares rises past the likelihood-ratio threshold, so an asymmetric interval
// comes out asymmetric.
struct ProfileInterval {
    double      lower = 0.0;
    double      upper = 0.0;
    bool        lowerFound = false;   // false when the profile never crossed within the search span
    bool        upperFound = false;
    double      threshold = 0.0;      // the sum-of-squares level that defines the bound
    double      confidence = 0.0;
    std::string note;
};

// `confidence` must be one of 0.68, 0.90, 0.95, 0.99: the threshold is the
// chi-square(1 df) quantile times the residual variance, and shipping an
// interpolated quantile table would be inventing statistics. Any other value is
// refused through ProfileInterval::note.
ProfileInterval profileLikelihood(const std::vector<double>& best, std::size_t observations,
                                  const LmEvaluate& evaluate, std::size_t index,
                                  double confidence = 0.95,
                                  const std::vector<double>& weights = {},
                                  const LmOptions& options = {});

// Second-order-corrected Akaike information criterion for a least-squares fit:
// AICc = n*ln(SSR/n) + 2k + 2k(k+1)/(n-k-1). Returns infinity when n <= k+1,
// where the correction term is undefined - a model comparison that cannot be made
// must not return a finite number that looks like one.
double aicc(double sumOfSquaredResiduals, std::size_t observations, std::size_t parameters);

// Coefficient of determination of a fitted series against observations. Returns 0
// when the observations have no variance, rather than dividing by zero.
double rSquared(const std::vector<double>& observed, const std::vector<double>& fitted);

}  // namespace biocad::numeric
