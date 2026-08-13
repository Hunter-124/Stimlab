#include "sim/Solvers.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>

#include <Eigen/Dense>

#include "numeric/Ode.h"

namespace biocad::sim {
namespace {

// ROS3P (Lang & Verwer 2001), three stages, third order, L-stable, in Hairer &
// Wanner's convention:
//
//   (1/(gamma*h) I - J) k_i = f(y + sum_{j<i} a_ij k_j) + (1/h) sum_{j<i} c_ij k_j
//   y_new = y + sum m_i k_i
//
// so ONE LU factorization of (1/(gamma*h) I - J) serves all three stages. That is
// the whole reason a Rosenbrock method is used here instead of a Newton-iterated BDF.
// The third-order accuracy of this table is verified numerically by a step-halving
// convergence study in the test suite, because a mistyped coefficient yields a
// method that still converges, just at the wrong order - which no smoke test catches.
//
// WHY THE ERROR ESTIMATE IS AGAINST THE FIRST-ORDER SOLUTION, NOT A SECOND-ORDER
// EMBEDDED ONE. Applied to y' = J y, this table's stage vectors satisfy the identity
//
//   (a21 - 1) k1 + k2 == 0    exactly, for every h,
//
// (the coefficient (3 - sqrt 3) + gamma*c21 is exactly zero) and every set of
// weights that is second-order accurate is forced to the same m3. Any second-order
// embedded solution therefore differs from m only along a direction that vanishes
// identically on a LINEAR problem, so its error estimate is exactly zero for
// A -> B and A <=> B - and the step size would then grow without bound on precisely
// the problems whose analytic solutions we check. The estimate is taken against the
// provably first-order solution y + (1/gamma) k1 instead, which is the
// linearly-implicit Euler step embedded in stage one. It is O(h^2) on linear and
// nonlinear problems alike, so the step control is conservative rather than blind.
constexpr double kGamma = 7.886751345948129e-01;
constexpr double kA21 = 1.2679491924311228;   // = 1 / gamma = 3 - sqrt(3)
constexpr double kA31 = 1.2679491924311228;
constexpr double kA32 = 0.0;
constexpr double kC21 = -1.6076951545867362;  // = 6*sqrt(3) - 12
constexpr double kC31 = -3.4641016151377544;  // = -2*sqrt(3)
constexpr double kC32 = -1.7320508075688772;  // = -sqrt(3)
constexpr double kM1 = 2.0;
constexpr double kM2 = 5.7735026918962584e-01;   // = sqrt(3)/3
constexpr double kM3 = 4.2264973081037416e-01;   // = 1 - sqrt(3)/3
constexpr double kMh1 = 1.2679491924311228;      // = 1 / gamma
constexpr double kMh2 = 0.0;
constexpr double kMh3 = 0.0;
// Order of the ERROR ESTIMATE, which is what sets the step-selection exponent.
constexpr double kEstimateOrder = 1.0;

// WHY A ROSENBROCK STEP CONSERVES EXACTLY. With M = (1/(gamma*h)) I - J and a
// conservation law y (so y^T S = 0, hence y^T J = y^T S dv/dc = 0), y^T M = a y^T
// with a = 1/(gamma*h), so y^T M^-1 = y^T / a. Every stage right-hand side is a
// combination of f(c) (which satisfies y^T f = 0) and earlier stage vectors, so
// y^T k_i = 0 for every stage and the conserved quantity is preserved to roundoff,
// not to the local error. The reported drift is therefore a genuine audit of the
// arithmetic rather than a restatement of the tolerance.

double now() {
    using clock = std::chrono::steady_clock;
    return std::chrono::duration<double>(clock::now().time_since_epoch()).count();
}

std::vector<double> outputGrid(double horizon, int points) {
    const int n = std::max(2, points);
    std::vector<double> t(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) t[static_cast<std::size_t>(i)] = horizon * i / (n - 1);
    return t;
}

struct Recorder {
    std::vector<std::string>         speciesIds;
    std::vector<std::vector<double>> traj;
    void reserve(std::size_t ns, std::size_t nt) {
        traj.assign(ns, {});
        for (auto& row : traj) row.reserve(nt);
    }
    void push(const std::vector<double>& y) {
        for (std::size_t i = 0; i < y.size(); ++i) traj[i].push_back(y[i]);
    }
};

void fillConservation(const NetworkSpec& spec, TimeCourse& tc) {
    if (spec.conservationLaws.empty()) return;
    const std::size_t nt = tc.times.size();
    for (const auto& law : spec.conservationLaws) {
        std::vector<double> series(nt, 0.0);
        for (std::size_t k = 0; k < nt; ++k) {
            double sum = 0;
            for (std::size_t i = 0; i < law.size() && i < tc.trajectories.size(); ++i)
                if (law[i] != 0) sum += law[i] * tc.trajectories[i][k];
            series[k] = sum;
        }
        const double base = series.empty() ? 0.0 : series.front();
        for (double v : series) {
            const double drift = std::abs(base) > 0 ? std::abs(v - base) / std::abs(base)
                                                    : std::abs(v - base);
            tc.worstConservationDrift = std::max(tc.worstConservationDrift, drift);
        }
        tc.conservedQuantities.push_back(std::move(series));
    }
}

// One stochastic reaction channel: a propensity constant, its reactant multiset and
// the state change it applies. A reversible mass-action step becomes TWO channels,
// because the SSA fires elementary events, not net rates.
struct Channel {
    double                                     k = 0;
    std::vector<std::pair<std::size_t, int>>   reactants;   // species, order
    std::vector<std::pair<std::size_t, double>> delta;      // species, change
    int                                        order = 0;
};

double propensity(const Channel& ch, const std::vector<double>& n) {
    double a = ch.k;
    for (const auto& [i, order] : ch.reactants) {
        const double x = n[i];
        // prod (n choose a) * a! = falling factorial, divided by a! for a >= 2:
        // the standard deterministic-to-stochastic conversion for identical reactants.
        double term = 1;
        for (int q = 0; q < order; ++q) term *= (x - q);
        for (int q = 2; q <= order; ++q) term /= q;
        a *= term;
        if (a <= 0) return 0;
    }
    return a > 0 ? a : 0.0;
}

}  // namespace

double rosenbrockStep(const Network& net, std::vector<double>& y, double h, double relTol,
                      double absTol, SolverReport& report) {
    const std::size_t n = y.size();
    std::vector<double> jflat;
    net.jacobian(y, jflat);
    ++report.jacobianEvaluations;

    Eigen::MatrixXd j(static_cast<Eigen::Index>(n), static_cast<Eigen::Index>(n));
    for (std::size_t r = 0; r < n; ++r)
        for (std::size_t c = 0; c < n; ++c)
            j(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c)) = jflat[r * n + c];

    Eigen::MatrixXd m = -j;
    const double scale = 1.0 / (kGamma * h);
    for (std::size_t i = 0; i < n; ++i) m(static_cast<Eigen::Index>(i),
                                          static_cast<Eigen::Index>(i)) += scale;
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(m);   // the single factorization per step

    auto toEigen = [&](const std::vector<double>& v) {
        Eigen::VectorXd e(static_cast<Eigen::Index>(n));
        for (std::size_t i = 0; i < n; ++i) e[static_cast<Eigen::Index>(i)] = v[i];
        return e;
    };
    std::vector<double> tmp(n), f(n);
    // The lambda returns an EAGER vector: Eigen's solve() is a lazy expression that
    // would reference `rhs` after it went out of scope.
    auto stage = [&](const Eigen::VectorXd& shift,
                     const Eigen::VectorXd& cterm) -> Eigen::VectorXd {
        for (std::size_t i = 0; i < n; ++i)
            tmp[i] = y[i] + shift[static_cast<Eigen::Index>(i)];
        net.derivatives(tmp, f);
        const Eigen::VectorXd rhs = toEigen(f) + cterm / h;
        return lu.solve(rhs);
    };

    const Eigen::VectorXd zero = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(n));
    const Eigen::VectorXd k1 = stage(zero, zero);
    const Eigen::VectorXd k2 = stage(kA21 * k1, kC21 * k1);
    const Eigen::VectorXd k3 = stage(kA31 * k1 + kA32 * k2, kC31 * k1 + kC32 * k2);

    const Eigen::VectorXd step = kM1 * k1 + kM2 * k2 + kM3 * k3;
    const Eigen::VectorXd err = (kM1 - kMh1) * k1 + (kM2 - kMh2) * k2 + (kM3 - kMh3) * k3;

    double norm = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const double yn = y[i] + step[static_cast<Eigen::Index>(i)];
        const double tol = absTol + relTol * std::max(std::abs(y[i]), std::abs(yn));
        norm = std::max(norm, std::abs(err[static_cast<Eigen::Index>(i)]) / tol);
    }
    for (std::size_t i = 0; i < n; ++i) y[i] += step[static_cast<Eigen::Index>(i)];
    return norm;
}

TimeCourse integrate(const NetworkSpec& spec, const IntegrationOptions& options) {
    TimeCourse tc;
    Network net;
    std::string error;
    if (!Network::compile(spec, net, &error)) {
        tc.warnings.push_back("network does not compile: " + error);
        tc.solver.note = error;
        return tc;
    }
    if (!spec.wegscheiderViolations.empty()) {
        // Integrating a thermodynamically inconsistent cycle produces a steady state
        // that consumes nothing and does work forever. Refusing is the honest answer.
        for (const auto& w : spec.wegscheiderViolations)
            tc.warnings.push_back("refused: " + w);
        tc.solver.note = "network violates a Wegscheider cycle condition";
        return tc;
    }

    tc.times = outputGrid(options.horizon, options.outputPoints);
    tc.speciesIds = net.speciesIds();
    Recorder rec;
    rec.reserve(net.speciesCount(), tc.times.size());

    std::vector<double> y = net.initialState();
    SolverReport& rep = tc.solver;
    rep.relativeTolerance = options.relativeTolerance;
    rep.absoluteTolerance = options.absoluteTolerance;
    const double t0 = now();

    auto clip = [&] {
        for (double& v : y) {
            if (v < 0) {
                v = 0;
                ++rep.nonNegativityClips;
            }
        }
    };

    if (options.method == "rk4") {
        rep.method = "rk4";
        const double span = tc.times.size() > 1 ? tc.times[1] - tc.times[0] : options.horizon;
        const double h = options.fixedStep > 0 ? options.fixedStep : span / 10.0;
        numeric::OdeDerivative f = [&net](double, const std::vector<double>& c,
                                          std::vector<double>& d) { net.derivatives(c, d); };
        rec.push(y);
        for (std::size_t k = 1; k < tc.times.size(); ++k) {
            std::int64_t before = rep.acceptedSteps;
            numeric::rk4Integrate(tc.times[k - 1], tc.times[k], h, y, f,
                                  [&](double, const std::vector<double>&) { ++rep.acceptedSteps; });
            if (rep.acceptedSteps == before) ++rep.acceptedSteps;
            clip();
            rec.push(y);
        }
        rep.note = "fixed-step explicit RK4: stability, not accuracy, sets the step on a stiff "
                   "network";
    } else {
        rep.method = "rosenbrock";
        const double maxStep = options.maxStep > 0 ? options.maxStep : options.horizon;
        double h = std::min(options.initialStep > 0 ? options.initialStep : 1e-6, maxStep);
        rec.push(y);
        double t = 0;
        bool failed = false;
        for (std::size_t k = 1; k < tc.times.size() && !failed; ++k) {
            const double target = tc.times[k];
            while (t < target - 1e-15 * std::max(1.0, target)) {
                double trial = std::min({h, maxStep, target - t});
                std::vector<double> candidate = y;
                const double norm = rosenbrockStep(net, candidate, trial,
                                                   options.relativeTolerance,
                                                   options.absoluteTolerance, rep);
                if (norm <= 1.0 || trial <= 1e-14 * std::max(1.0, target)) {
                    y = candidate;
                    t += trial;
                    ++rep.acceptedSteps;
                    clip();
                    // Grow only the proposed step, never the one just taken to hit an
                    // output point, or the grid would throttle the whole integration.
                    const double exponent = -1.0 / (kEstimateOrder + 1.0);
                    const double factor = norm > 0 ? 0.9 * std::pow(norm, exponent) : 5.0;
                    h = std::min(maxStep, h * std::clamp(factor, 0.2, 5.0));
                } else {
                    ++rep.rejectedSteps;
                    h = trial * std::clamp(
                                    0.9 * std::pow(norm, -1.0 / (kEstimateOrder + 1.0)), 0.1, 0.9);
                    if (h < 1e-15 * std::max(1.0, target)) {
                        tc.warnings.push_back("step size collapsed below 1e-15 at t = " +
                                              std::to_string(t));
                        rep.note = "step size collapsed";
                        failed = true;
                        break;
                    }
                }
                if (rep.acceptedSteps + rep.rejectedSteps > options.maxSteps) {
                    tc.warnings.push_back("step budget exhausted at t = " + std::to_string(t));
                    rep.note = "step budget exhausted";
                    failed = true;
                    break;
                }
            }
            rec.push(y);
        }
        if (rep.note.empty())
            rep.note = "ROS3P, L-stable, one LU factorization per step";
    }

    rep.cpuSeconds = now() - t0;
    tc.trajectories = std::move(rec.traj);
    // Pad a truncated run so every trajectory has one value per reported time: a
    // ragged result would silently misalign a plot.
    for (auto& row : tc.trajectories)
        while (row.size() < tc.times.size()) row.push_back(row.empty() ? 0.0 : row.back());
    fillConservation(spec, tc);
    if (rep.nonNegativityClips > 0)
        tc.warnings.push_back(
            "WARNING: " + std::to_string(rep.nonNegativityClips) +
            " negative concentration(s) were clipped to zero. That is a step-size failure "
            "being reported, not repaired - tighten the tolerance or use the Rosenbrock "
            "solver before trusting this curve.");
    return tc;
}

std::uint64_t poissonVariate(Pcg64Dxsm& rng, double lambda) {
    if (!(lambda > 0)) return 0;
    if (lambda < 30.0) {
        // Knuth's product method: exact, and its expected iteration count is
        // lambda + 1, which is why it is abandoned above 30.
        const double limit = std::exp(-lambda);
        double p = 1.0;
        std::uint64_t k = 0;
        while (true) {
            p *= rng.nextDouble();
            if (p <= limit) return k;
            ++k;
            if (k > 1'000'000) return k;   // defensive; unreachable for lambda < 30
        }
    }
    // Hormann (1993) PTRS transformed rejection: constant expected work in lambda.
    const double b = 0.931 + 2.53 * std::sqrt(lambda);
    const double a = -0.059 + 0.02483 * b;
    const double invAlpha = 1.1239 + 1.1328 / (b - 3.4);
    const double vr = 0.9277 - 3.6224 / (b - 2.0);
    const double logLambda = std::log(lambda);
    while (true) {
        const double u = rng.nextDouble() - 0.5;
        const double v = rng.nextDouble();
        const double us = 0.5 - std::abs(u);
        const double kd = std::floor((2.0 * a / us + b) * u + lambda + 0.43);
        if (us >= 0.07 && v <= vr) return static_cast<std::uint64_t>(kd);
        if (kd < 0 || (us < 0.013 && v > us)) continue;
        if (std::log(v * invAlpha / (a / (us * us) + b)) <=
            kd * logLambda - lambda - std::lgamma(kd + 1.0))
            return static_cast<std::uint64_t>(kd);
    }
}

StochasticEnsemble stochastic(const NetworkSpec& spec, const StochasticOptions& options) {
    StochasticEnsemble out;
    out.solver.method = options.tauLeap ? "tau-leap" : "gillespie";
    out.solver.seed = options.seed;
    Network net;
    std::string error;
    if (!Network::compile(spec, net, &error)) {
        out.solver.note = "network does not compile: " + error;
        return out;
    }

    // Build the elementary channels. A saturable rate law is not an elementary
    // reaction, so it is refused BY NAME rather than approximated by its
    // deterministic rate - which would put a mean-field rate inside an exact
    // stochastic algorithm and quietly destroy the variance the SSA exists to give.
    const std::size_t ns = net.speciesCount(), nr = net.reactionCount();
    const std::vector<double>& s = net.stoichiometry();
    std::vector<Channel> channels;
    for (std::size_t r = 0; r < nr; ++r) {
        const ReactionSpec& rs = spec.reactions[r];
        if (rs.law != RateLaw::MassAction && rs.law != RateLaw::ReversibleMassAction) {
            out.solver.note = "reaction '" + rs.id +
                              "' uses a saturable rate law, which is not an elementary reaction "
                              "channel; the stochastic simulation algorithm is refused rather "
                              "than approximated";
            return out;
        }
        auto build = [&](bool forward) {
            Channel ch;
            ch.k = forward ? rs.parameters[0] : rs.parameters[1];
            const auto& side = forward ? rs.reactants : rs.products;
            for (const auto& [id, stoich] : side) {
                const auto it = std::find(net.speciesIds().begin(), net.speciesIds().end(), id);
                const std::size_t idx =
                    static_cast<std::size_t>(std::distance(net.speciesIds().begin(), it));
                ch.reactants.emplace_back(idx, static_cast<int>(std::lround(stoich)));
                ch.order += static_cast<int>(std::lround(stoich));
            }
            const double sign = forward ? 1.0 : -1.0;
            for (std::size_t i = 0; i < ns; ++i) {
                const double d = sign * s[i * nr + r];
                if (d != 0 && !net.boundary(i)) ch.delta.emplace_back(i, d);
            }
            channels.push_back(std::move(ch));
        };
        build(true);
        if (rs.law == RateLaw::ReversibleMassAction) build(false);
    }

    out.times = outputGrid(options.horizon, options.outputPoints);
    const std::size_t nt = out.times.size();
    out.speciesIds = net.speciesIds();
    out.replicates = std::max(1, options.replicates);
    std::vector<std::vector<double>> sum(ns, std::vector<double>(nt, 0.0));
    std::vector<std::vector<double>> sumSq(ns, std::vector<double>(nt, 0.0));

    // Highest reaction order each species participates in: the g_i of the
    // Cao-Gillespie-Petzold tau selection.
    std::vector<int> highestOrder(ns, 1);
    for (const Channel& ch : channels)
        for (const auto& [i, order] : ch.reactants)
            highestOrder[i] = std::max(highestOrder[i], ch.order);

    const double t0 = now();
    std::int64_t events = 0;
    for (int rep = 0; rep < out.replicates; ++rep) {
        // One stream per replicate, derived from the run seed and the replicate
        // index, so a replicate is reproducible on its own and adding replicates
        // never perturbs the ones already computed.
        Pcg64Dxsm rng(options.seed, static_cast<std::uint64_t>(rep) + 1u);
        std::vector<double> n(ns);
        for (std::size_t i = 0; i < ns; ++i) n[i] = std::round(net.initialState()[i]);
        double t = 0;
        std::size_t next = 0;
        auto record = [&](const std::vector<double>& state) {
            for (std::size_t i = 0; i < ns; ++i) {
                sum[i][next] += state[i];
                sumSq[i][next] += state[i] * state[i];
            }
            ++next;
        };
        // The state is recorded for every grid point strictly BEFORE the next event
        // time, using the state that actually held over that interval. Recording
        // after applying the event instead biases every sample by one reaction -
        // which showed up as a systematic ~1-molecule deficit against the exact
        // binomial mean, well outside Monte Carlo error.
        std::vector<double> a(channels.size(), 0.0);
        while (next < nt) {
            double a0 = 0;
            for (std::size_t c = 0; c < channels.size(); ++c) {
                a[c] = propensity(channels[c], n);
                a0 += a[c];
            }
            if (a0 <= 0) break;   // absorbing state: the tail is filled below
            double tNext = t;
            std::vector<double> updated = n;
            bool leapt = false;
            if (options.tauLeap) {
                double tau = std::numeric_limits<double>::infinity();
                for (std::size_t i = 0; i < ns; ++i) {
                    double mu = 0, sigma2 = 0;
                    for (std::size_t c = 0; c < channels.size(); ++c) {
                        for (const auto& [j, d] : channels[c].delta) {
                            if (j != i) continue;
                            mu += d * a[c];
                            sigma2 += d * d * a[c];
                        }
                    }
                    if (mu == 0 && sigma2 == 0) continue;
                    const double bound =
                        std::max(options.epsilon * n[i] / highestOrder[i], 1.0);
                    if (mu != 0) tau = std::min(tau, bound / std::abs(mu));
                    if (sigma2 > 0) tau = std::min(tau, bound * bound / sigma2);
                }
                // Below 10 mean events a leap is neither cheaper nor as accurate as
                // exact steps, which is the published switch condition.
                if (tau > 10.0 / a0 && std::isfinite(tau)) {
                    // A leap must never straddle an output time: the recorded value
                    // would then be the state at the START of the leap, which biased
                    // the measured mean by ~0.7% and inflated the measured variance
                    // 13-fold against the exact binomial solution at 200000 copies.
                    double nextGrid = out.times.back();
                    for (std::size_t g = next; g < nt; ++g)
                        if (out.times[g] > t) { nextGrid = out.times[g]; break; }
                    if (nextGrid > t) tau = std::min(tau, nextGrid - t);
                    for (int attempt = 0; attempt < 20 && !leapt; ++attempt) {
                        std::vector<double> trial = n;
                        bool negative = false;
                        for (std::size_t c = 0; c < channels.size() && !negative; ++c) {
                            const std::uint64_t fired = poissonVariate(rng, a[c] * tau);
                            if (fired == 0) continue;
                            for (const auto& [j, d] : channels[c].delta) {
                                trial[j] += d * static_cast<double>(fired);
                                if (trial[j] < 0) negative = true;
                            }
                        }
                        if (negative) {
                            tau *= 0.5;
                            continue;
                        }
                        updated = std::move(trial);
                        tNext = t + tau;
                        leapt = true;
                        ++out.solver.acceptedSteps;
                    }
                }
            }
            if (!leapt) {
                const double u1 = std::max(rng.nextDouble(), 1e-300);
                tNext = t - std::log(u1) / a0;
                const double pick = rng.nextDouble() * a0;
                double acc = 0;
                std::size_t chosen = channels.size() - 1;
                for (std::size_t c = 0; c < channels.size(); ++c) {
                    acc += a[c];
                    if (acc >= pick) { chosen = c; break; }
                }
                for (const auto& [j, d] : channels[chosen].delta) updated[j] += d;
                ++out.solver.acceptedSteps;
            }
            while (next < nt && out.times[next] < tNext) record(n);
            n = std::move(updated);
            t = tNext;
            if (++events > options.maxEvents) {
                out.solver.note = "event budget exhausted";
                break;
            }
        }
        // Any grid point past the last event keeps the final state.
        while (next < nt) {
            for (std::size_t i = 0; i < ns; ++i) {
                sum[i][next] += n[i];
                sumSq[i][next] += n[i] * n[i];
            }
            ++next;
        }
    }

    const double m = static_cast<double>(out.replicates);
    out.mean.assign(ns, std::vector<double>(nt, 0.0));
    out.variance.assign(ns, std::vector<double>(nt, 0.0));
    for (std::size_t i = 0; i < ns; ++i) {
        for (std::size_t k = 0; k < nt; ++k) {
            const double mean = sum[i][k] / m;
            out.mean[i][k] = mean;
            // Unbiased sample variance: the ensemble is a sample, and dividing by m
            // would understate the spread that is the point of running it.
            out.variance[i][k] =
                m > 1 ? std::max(0.0, (sumSq[i][k] - m * mean * mean) / (m - 1.0)) : 0.0;
        }
    }
    out.solver.cpuSeconds = now() - t0;
    if (out.solver.note.empty())
        out.solver.note = options.tauLeap
                              ? "explicit tau-leaping, Cao-Gillespie-Petzold step selection"
                              : "exact stochastic simulation algorithm (Gillespie direct method)";
    return out;
}

}  // namespace biocad::sim
