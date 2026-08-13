// tests/test_network.cpp - the reaction-network core: stoichiometry, the analytic
// Jacobian, conservation and Wegscheider structure, both deterministic solvers, the
// stochastic algorithms, and the constraint-based flux path.
//
// The three checks worth reading twice:
//   * the analytic Jacobian is compared against a central difference at a state where
//     a species is EXACTLY zero, which is the state where the tempting v*a/c form
//     divides by zero and where the true derivative is finite and nonzero;
//   * the Robertson problem is integrated by both solvers, and the explicit one is
//     asserted to FAIL, so the reason the Rosenbrock method exists is in the suite
//     rather than in a comment;
//   * the Gillespie ensemble is checked against the exact binomial solution of a
//     first-order decay - mean within Monte Carlo error, variance within a few
//     percent - because a stochastic solver that gets the mean right and the variance
//     wrong has destroyed the only thing it was run for.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "sim/Control.h"
#include "sim/Flux.h"
#include "sim/Network.h"
#include "sim/Solvers.h"

using namespace biocad;
using namespace biocad::sim;

namespace {

SpeciesSpec sp(const char* id, double c0, bool boundary = false, const char* formula = "") {
    SpeciesSpec s;
    s.id = id;
    s.name = id;
    s.initialConcentration = c0;
    s.boundary = boundary;
    s.formula = formula;
    s.compartment = "cell";
    return s;
}

ReactionSpec rx(const char* id, std::vector<std::pair<std::string, double>> reactants,
                std::vector<std::pair<std::string, double>> products, RateLaw law,
                std::vector<double> parameters) {
    ReactionSpec r;
    r.id = id;
    r.reactants = std::move(reactants);
    r.products = std::move(products);
    r.law = law;
    r.parameters = std::move(parameters);
    r.reversible = law == RateLaw::ReversibleMassAction ||
                   law == RateLaw::ReversibleMichaelisMenten;
    return r;
}

double worstAbsoluteError(const TimeCourse& tc, std::size_t species,
                          const std::function<double(double)>& exact) {
    double worst = 0;
    for (std::size_t k = 0; k < tc.times.size(); ++k)
        worst = std::max(worst, std::abs(tc.trajectories[species][k] - exact(tc.times[k])));
    return worst;
}

IntegrationOptions tight(double horizon, int points, const char* method = "rosenbrock") {
    IntegrationOptions o;
    o.horizon = horizon;
    o.outputPoints = points;
    o.method = method;
    o.relativeTolerance = 1e-8;
    o.absoluteTolerance = 1e-12;
    o.fixedStep = 1e-3;
    return o;
}

// The Robertson problem, written as reactions. y3 is a catalyst in R2 and y2 in R3,
// which is exactly what makes y1 + y2 + y3 conserved in the classical right-hand side.
NetworkSpec robertson() {
    NetworkSpec n;
    n.id = "robertson";
    n.species = {sp("y1", 1.0), sp("y2", 0.0), sp("y3", 0.0)};
    n.reactions = {
        rx("R1", {{"y1", 1}}, {{"y2", 1}}, RateLaw::MassAction, {0.04}),
        rx("R2", {{"y2", 1}, {"y3", 1}}, {{"y1", 1}, {"y3", 1}}, RateLaw::MassAction, {1e4}),
        rx("R3", {{"y2", 2}}, {{"y3", 1}, {"y2", 1}}, RateLaw::MassAction, {3e7})};
    return analyze(n);
}

// A hand-solvable fermentation cartoon in which every reaction conserves C, H and O.
// Reducing equivalents are carried explicitly as H2, so the redox balance is real.
//   v1: Glc_e -> Glc            [0, 10]
//   v2: Glc -> 2 Pyr + 2 H2     C6H12O6 -> 2 C3H4O3 + 2 H2
//   v3: Pyr + H2 -> Lac         forced to at least 3 (a maintenance demand)
//   v4: Pyr -> AcAld + CO2      C3H4O3 -> C2H4O + CO2
//   v5..v8: exports
// Hand solution: v1 = 10 => Pyr = 20; v3 = 3 => v4 = v6 = 17.
NetworkSpec fermentation() {
    NetworkSpec n;
    n.id = "toy-fermentation";
    n.species = {sp("Glc_e", 0, true, "C6H12O6"), sp("Glc", 0, false, "C6H12O6"),
                 sp("Pyr", 0, false, "C3H4O3"),   sp("Lac", 0, false, "C3H6O3"),
                 sp("AcAld", 0, false, "C2H4O"),  sp("CO2", 0, false, "CO2"),
                 sp("H2", 0, false, "H2"),        sp("Lac_e", 0, true, "C3H6O3"),
                 sp("AcAld_e", 0, true, "C2H4O"), sp("CO2_e", 0, true, "CO2"),
                 sp("H2_e", 0, true, "H2")};
    const std::vector<double> k = {1.0};
    n.reactions = {rx("v1", {{"Glc_e", 1}}, {{"Glc", 1}}, RateLaw::MassAction, k),
                   rx("v2", {{"Glc", 1}}, {{"Pyr", 2}, {"H2", 2}}, RateLaw::MassAction, k),
                   rx("v3", {{"Pyr", 1}, {"H2", 1}}, {{"Lac", 1}}, RateLaw::MassAction, k),
                   rx("v4", {{"Pyr", 1}}, {{"AcAld", 1}, {"CO2", 1}}, RateLaw::MassAction, k),
                   rx("v5", {{"Lac", 1}}, {{"Lac_e", 1}}, RateLaw::MassAction, k),
                   rx("v6", {{"AcAld", 1}}, {{"AcAld_e", 1}}, RateLaw::MassAction, k),
                   rx("v7", {{"CO2", 1}}, {{"CO2_e", 1}}, RateLaw::MassAction, k),
                   rx("v8", {{"H2", 1}}, {{"H2_e", 1}}, RateLaw::MassAction, k)};
    return n;
}

const std::vector<FluxBound> kFermentationBounds = {{"v1", 0, 10}, {"v3", 3, 1000}};

}  // namespace

TEST_CASE("A -> B reproduces its analytic solution and its conservation law", "[network]") {
    const double k = 0.7;
    NetworkSpec n;
    n.id = "ab";
    n.species = {sp("A", 1.0), sp("B", 0.0)};
    n.reactions = {rx("R1", {{"A", 1}}, {{"B", 1}}, RateLaw::MassAction, {k})};
    n = analyze(n);

    REQUIRE(n.conservationLaws.size() == 1);
    REQUIRE(n.conservationLabels[0] == "A + B");
    REQUIRE(n.thermodynamicCycles.empty());

    const auto a = [k](double t) { return std::exp(-k * t); };
    const auto b = [k](double t) { return 1.0 - std::exp(-k * t); };
    for (const char* method : {"rosenbrock", "rk4"}) {
        const TimeCourse tc = integrate(n, tight(5.0, 51, method));
        REQUIRE(tc.solver.method == std::string(method));
        REQUIRE(worstAbsoluteError(tc, 0, a) < 1e-10);
        REQUIRE(worstAbsoluteError(tc, 1, b) < 1e-10);
        // Conservation is exact for both methods, not merely within the tolerance.
        REQUIRE(tc.worstConservationDrift < 1e-12);
        REQUIRE(tc.solver.nonNegativityClips == 0);
    }
}

TEST_CASE("A <=> B reproduces its analytic relaxation to equilibrium", "[network]") {
    const double kf = 0.3, kr = 0.7;
    NetworkSpec n;
    n.id = "rev";
    n.species = {sp("A", 1.0), sp("B", 0.0)};
    n.reactions = {rx("R1", {{"A", 1}}, {{"B", 1}}, RateLaw::ReversibleMassAction, {kf, kr})};
    n = analyze(n);
    REQUIRE(n.conservationLaws.size() == 1);

    const auto exact = [kf, kr](double t) {
        return (kr + kf * std::exp(-(kf + kr) * t)) / (kf + kr);
    };
    for (const char* method : {"rosenbrock", "rk4"}) {
        const TimeCourse tc = integrate(n, tight(8.0, 81, method));
        REQUIRE(worstAbsoluteError(tc, 0, exact) < 1e-10);
    }
}

TEST_CASE("A -> B -> C reproduces the consecutive first-order solution", "[network]") {
    const double k1 = 1.0, k2 = 0.3;
    NetworkSpec n;
    n.id = "chain";
    n.species = {sp("A", 1.0), sp("B", 0.0), sp("C", 0.0)};
    n.reactions = {rx("R1", {{"A", 1}}, {{"B", 1}}, RateLaw::MassAction, {k1}),
                   rx("R2", {{"B", 1}}, {{"C", 1}}, RateLaw::MassAction, {k2})};
    n = analyze(n);
    REQUIRE(n.conservationLaws.size() == 1);

    const auto b = [k1, k2](double t) {
        return k1 / (k2 - k1) * (std::exp(-k1 * t) - std::exp(-k2 * t));
    };
    const auto c = [k1, k2](double t) {
        return 1.0 + (k1 * std::exp(-k2 * t) - k2 * std::exp(-k1 * t)) / (k2 - k1);
    };
    for (const char* method : {"rosenbrock", "rk4"}) {
        IntegrationOptions o = tight(12.0, 121, method);
        o.fixedStep = 5e-4;
        const TimeCourse tc = integrate(n, o);
        REQUIRE(worstAbsoluteError(tc, 1, b) < 1e-10);
        REQUIRE(worstAbsoluteError(tc, 2, c) < 1e-10);
    }
}

TEST_CASE("the ROS3P table really is third order", "[network]") {
    // 2A -> B has the analytic solution A(t) = A0 / (1 + 2 k A0 t). A mistyped
    // coefficient still converges, just at the wrong order, so the ORDER is the check.
    NetworkSpec n;
    n.id = "order";
    n.species = {sp("A", 1.0), sp("B", 0.0)};
    n.reactions = {rx("R1", {{"A", 2}}, {{"B", 1}}, RateLaw::MassAction, {0.5})};
    Network net;
    std::string error;
    REQUIRE(Network::compile(n, net, &error));

    const double horizon = 4.0;
    const double exact = 1.0 / (1.0 + 2.0 * 0.5 * 1.0 * horizon);
    double previous = 0;
    for (int p = 0; p < 5; ++p) {
        const double h = 0.25 / std::pow(2.0, p);
        std::vector<double> y = net.initialState();
        SolverReport report;
        double t = 0;
        while (t < horizon - 1e-12) {
            const double step = std::min(h, horizon - t);
            rosenbrockStep(net, y, step, 1e-30, 1e-30, report);
            t += step;
        }
        const double err = std::abs(y[0] - exact);
        if (p >= 2) {
            const double order = std::log2(previous / err);
            REQUIRE(order > 2.7);
            REQUIRE(order < 3.4);
        }
        previous = err;
    }
}

TEST_CASE("the Robertson problem is conserved by Rosenbrock and destroyed by RK4",
          "[network]") {
    const NetworkSpec n = robertson();
    REQUIRE(n.conservationLaws.size() == 1);
    REQUIRE(n.conservationLabels[0] == "y1 + y2 + y3");

    IntegrationOptions o = tight(40.0, 41);
    o.initialStep = 1e-8;
    const TimeCourse ros = integrate(n, o);
    REQUIRE(ros.solver.method == "rosenbrock");
    REQUIRE(ros.worstConservationDrift < 1e-10);
    REQUIRE(ros.solver.nonNegativityClips == 0);
    const double sum = ros.trajectories[0].back() + ros.trajectories[1].back() +
                       ros.trajectories[2].back();
    REQUIRE(std::abs(sum - 1.0) < 1e-10);
    // One LU factorization per step, shared by all three stages.
    REQUIRE(ros.solver.jacobianEvaluations ==
            ros.solver.acceptedSteps + ros.solver.rejectedSteps);

    // Explicit RK4 at a step a person would actually choose does not merely lose
    // accuracy: it leaves the invariant manifold entirely.
    IntegrationOptions explicitOptions = tight(40.0, 41, "rk4");
    explicitOptions.fixedStep = 1e-3;
    explicitOptions.maxSteps = 100000000;
    const TimeCourse rk = integrate(n, explicitOptions);
    REQUIRE_FALSE(rk.worstConservationDrift < 1e-10);

    // The step it would need: the explicit stability limit is about 2.78/|lambda|,
    // and the Jacobian's largest diagonal entry just after the transient is ~2.3e3,
    // so a 40-second horizon costs tens of thousands of steps at best.
    Network net;
    std::string error;
    REQUIRE(Network::compile(n, net, &error));
    std::vector<double> j;
    net.jacobian({0.9, 2.2e-5, 0.1}, j);
    double biggest = 0;
    for (std::size_t i = 0; i < 3; ++i) biggest = std::max(biggest, std::abs(j[i * 3 + i]));
    REQUIRE(biggest > 1e3);
    REQUIRE(40.0 / (2.78 / biggest) > 1e4);
}

TEST_CASE("the analytic Jacobian is exact where a species is depleted", "[network]") {
    // y2 = 0 EXACTLY is the state where dv/dc_j = v * a_j / c_j divides by zero while
    // the true derivative is finite and nonzero. Every rate law is present so the whole
    // excluding-product implementation is covered at once.
    NetworkSpec n;
    n.id = "jac";
    n.species = {sp("y1", 1.0), sp("y2", 0.0), sp("y3", 0.0), sp("E", 0.5, true)};
    n.reactions = {
        rx("R1", {{"y1", 1}}, {{"y2", 1}}, RateLaw::MassAction, {0.04}),
        rx("R2", {{"y2", 1}, {"y3", 1}}, {{"y1", 1}, {"y3", 1}}, RateLaw::MassAction, {1e4}),
        rx("R3", {{"y2", 2}}, {{"y3", 1}, {"y2", 1}}, RateLaw::MassAction, {3e7}),
        rx("R4", {{"y1", 1}}, {{"y3", 1}}, RateLaw::MichaelisMenten, {2.0, 0.3}),
        rx("R5", {{"y3", 1}}, {{"y1", 1}}, RateLaw::Hill, {1.5, 0.4, 2.0}),
        rx("R6", {{"y1", 1}}, {{"y2", 1}}, RateLaw::ReversibleMichaelisMenten,
           {1.0, 0.2, 0.5, 0.4}),
        rx("R7", {{"E", 1}, {"y1", 1}}, {{"y3", 2}}, RateLaw::ReversibleMassAction, {0.9, 0.35})};
    Network net;
    std::string error;
    REQUIRE(Network::compile(n, net, &error));

    const std::size_t ns = net.speciesCount();
    const std::vector<double> c = {0.8, 0.0, 0.25, 0.5};
    std::vector<double> analytic;
    net.jacobian(c, analytic);
    double worst = 0;
    for (std::size_t k = 0; k < ns; ++k) {
        const double h = 1e-6 * std::max(1.0, std::abs(c[k]));
        std::vector<double> plus = c, minus = c, fp, fm;
        plus[k] += h;
        minus[k] -= h;
        net.derivatives(plus, fp);
        net.derivatives(minus, fm);
        for (std::size_t i = 0; i < ns; ++i) {
            const double fd = (fp[i] - fm[i]) / (2 * h);
            worst = std::max(worst, std::abs(fd - analytic[i * ns + k]) /
                                        std::max(1.0, std::abs(fd)));
        }
    }
    REQUIRE(worst < 1e-7);
    // The derivative with respect to the depleted species is finite and NONZERO.
    REQUIRE(std::abs(analytic[0 * ns + 1]) > 1.0);
    REQUIRE(std::isfinite(analytic[1 * ns + 1]));
    // A boundary species is clamped, so its whole row is exactly zero.
    for (std::size_t k = 0; k < ns; ++k) REQUIRE(analytic[3 * ns + k] == 0.0);
}

TEST_CASE("a violated Wegscheider cycle condition is an error, not a warning", "[network]") {
    NetworkSpec n;
    n.id = "cycle";
    n.species = {sp("A", 1.0), sp("B", 0.0), sp("C", 0.0)};
    n.reactions = {rx("R1", {{"A", 1}}, {{"B", 1}}, RateLaw::ReversibleMassAction, {2.0, 1.0}),
                   rx("R2", {{"B", 1}}, {{"C", 1}}, RateLaw::ReversibleMassAction, {3.0, 1.0}),
                   rx("R3", {{"C", 1}}, {{"A", 1}}, RateLaw::ReversibleMassAction, {5.0, 1.0})};
    const NetworkSpec bad = analyze(n);
    REQUIRE(bad.thermodynamicCycles.size() == 1);
    REQUIRE(bad.wegscheiderViolations.size() == 1);
    // Integration is refused outright: a cycle that does net work at equilibrium
    // produces a smooth, plausible, thermodynamically impossible steady state.
    const TimeCourse refused = integrate(bad, tight(10.0, 11));
    REQUIRE(refused.times.empty());
    REQUIRE_FALSE(refused.warnings.empty());

    // 2 * 3 * (1/6) = 1 is consistent.
    n.reactions[2].parameters = {1.0, 6.0};
    const NetworkSpec good = analyze(n);
    REQUIRE(good.wegscheiderViolations.empty());
    REQUIRE_FALSE(integrate(good, tight(10.0, 11)).times.empty());
}

TEST_CASE("the Gillespie ensemble matches the exact binomial solution", "[network]") {
    // A first-order decay from n0 molecules leaves EXACTLY Binomial(n0, exp(-k t))
    // molecules at time t: mean n0.p, variance n0.p.(1-p).
    const double k = 0.5, n0 = 100;
    NetworkSpec n;
    n.id = "birth-death";
    n.species = {sp("A", n0), sp("B", 0.0)};
    n.reactions = {rx("R1", {{"A", 1}}, {{"B", 1}}, RateLaw::MassAction, {k})};

    StochasticOptions o;
    o.horizon = 4.0;
    o.outputPoints = 5;
    o.replicates = 4000;
    o.seed = 12345;
    const StochasticEnsemble e = stochastic(n, o);
    REQUIRE(e.solver.method == "gillespie");
    REQUIRE(e.replicates == 4000);
    for (std::size_t i = 0; i < e.times.size(); ++i) {
        const double p = std::exp(-k * e.times[i]);
        const double mean = n0 * p, variance = n0 * p * (1 - p);
        if (variance <= 0) continue;
        const double standardError = std::sqrt(variance / e.replicates);
        REQUIRE(std::abs(e.mean[0][i] - mean) < 4 * standardError);
        // The variance is the reason to run an SSA at all; 8% over 4000 replicates is
        // a bound the run must clear, not a threshold tuned to the observed value.
        REQUIRE(std::abs(e.variance[0][i] / variance - 1.0) < 0.08);
        // A + B is conserved exactly by every event.
        REQUIRE(std::abs(e.mean[0][i] + e.mean[1][i] - n0) < 1e-12);
    }

    // The same seed reproduces byte-identically; a different seed does not.
    const StochasticEnsemble again = stochastic(n, o);
    for (std::size_t i = 0; i < e.mean.size(); ++i)
        for (std::size_t j = 0; j < e.mean[i].size(); ++j)
            REQUIRE(e.mean[i][j] == again.mean[i][j]);
    o.seed = 12346;
    const StochasticEnsemble other = stochastic(n, o);
    bool differs = false;
    for (std::size_t j = 0; j < e.mean[0].size(); ++j)
        if (e.mean[0][j] != other.mean[0][j]) differs = true;
    REQUIRE(differs);
}

TEST_CASE("tau-leaping is cheaper, and its bias tracks its own error parameter",
          "[network]") {
    const double k = 0.5, n0 = 200000;
    NetworkSpec n;
    n.id = "big-decay";
    n.species = {sp("A", n0), sp("B", 0.0)};
    n.reactions = {rx("R1", {{"A", 1}}, {{"B", 1}}, RateLaw::MassAction, {k})};

    StochasticOptions o;
    o.horizon = 4.0;
    o.outputPoints = 5;
    o.replicates = 40;
    o.seed = 99;
    const StochasticEnsemble exact = stochastic(n, o);
    // The EXACT algorithm has to agree with the analytic mean to Monte Carlo error.
    for (std::size_t i = 1; i < exact.times.size(); ++i) {
        const double p = std::exp(-k * exact.times[i]);
        const double mean = n0 * p, variance = n0 * p * (1 - p);
        REQUIRE(std::abs(exact.mean[0][i] - mean) < 4 * std::sqrt(variance / o.replicates));
    }

    double previousBias = 1.0;
    for (double epsilon : {0.03, 0.005}) {
        o.tauLeap = true;
        o.epsilon = epsilon;
        const StochasticEnsemble leap = stochastic(n, o);
        REQUIRE(leap.solver.method == "tau-leap");
        // At 200000 copies a leap is legitimate and enormously cheaper. (At 100 copies
        // the published switch condition correctly refuses to leap at all.)
        REQUIRE(leap.solver.acceptedSteps * 5 < exact.solver.acceptedSteps);
        double worstBias = 0;
        for (std::size_t i = 1; i < leap.times.size(); ++i) {
            const double p = std::exp(-k * leap.times[i]);
            const double mean = n0 * p, variance = n0 * p * (1 - p);
            worstBias = std::max(worstBias, std::abs((leap.mean[0][i] - mean) / mean));
            // 40 replicates give a variance estimate with a ~23% relative standard
            // error, so 35% is the sampling bound, not a tuned tolerance.
            REQUIRE(std::abs(leap.variance[0][i] / variance - 1.0) < 0.35);
        }
        REQUIRE(worstBias < 2.0 * epsilon);
        REQUIRE(worstBias < previousBias);
        previousBias = worstBias;
    }
}

TEST_CASE("a saturable rate law is refused by name by the stochastic algorithm",
          "[network]") {
    NetworkSpec n;
    n.id = "not-elementary";
    n.species = {sp("S", 100), sp("P", 0)};
    n.reactions = {rx("E1", {{"S", 1}}, {{"P", 1}}, RateLaw::MichaelisMenten, {2.0, 0.5})};
    const StochasticEnsemble e = stochastic(n, {});
    REQUIRE(e.mean.empty());
    REQUIRE(e.solver.note.find("E1") != std::string::npos);
    REQUIRE(e.solver.note.find("elementary") != std::string::npos);
}

TEST_CASE("the LP solves hand-computed problems and names infeasible and unbounded ones",
          "[network][flux]") {
    // max 3x + 5y subject to x + 2y <= 14, 3x - y >= 0, x - y <= 2. Vertex
    // enumeration puts the optimum at x = 6, y = 4, value 38.
    LpProblem lp;
    lp.variables = 5;
    lp.constraints = 3;
    lp.a = {1, 2, 1, 0, 0, 3, -1, 0, -1, 0, 1, -1, 0, 0, 1};
    lp.b = {14, 0, 2};
    lp.c = {-3, -5, 0, 0, 0};
    lp.lower.assign(5, 0.0);
    lp.upper.assign(5, kInfiniteBound);
    const LpResult r = solveLp(lp);
    REQUIRE(r.feasible);
    REQUIRE(r.bounded);
    REQUIRE(std::abs(r.x[0] - 6.0) < 1e-9);
    REQUIRE(std::abs(r.x[1] - 4.0) < 1e-9);
    REQUIRE(std::abs(r.objective + 38.0) < 1e-9);

    LpProblem infeasible;
    infeasible.variables = 2;
    infeasible.constraints = 2;
    infeasible.a = {1, -1, 1, 1};
    infeasible.b = {3, 1};
    infeasible.c = {1, 0};
    infeasible.lower.assign(2, 0.0);
    infeasible.upper.assign(2, kInfiniteBound);
    REQUIRE(solveLp(infeasible).status == "infeasible");

    LpProblem unbounded;
    unbounded.variables = 2;
    unbounded.constraints = 1;
    unbounded.a = {0, 1};
    unbounded.b = {0};
    unbounded.c = {-1, 0};
    unbounded.lower.assign(2, 0.0);
    unbounded.upper.assign(2, kInfiniteBound);
    REQUIRE(solveLp(unbounded).status == "unbounded");
}

TEST_CASE("FBA reproduces the hand-solved toy model", "[network][flux]") {
    const NetworkSpec n = fermentation();
    REQUIRE(balance(n).massBalanced);

    const FluxSolution f = fba(n, "v6", kFermentationBounds);
    REQUIRE(f.feasible);
    REQUIRE(f.massBalanced);
    REQUIRE(std::abs(f.objectiveValue - 17.0) < 1e-9);
    REQUIRE(std::abs(f.fluxes[0] - 10.0) < 1e-9);   // v1 uptake at its ceiling
    REQUIRE(std::abs(f.fluxes[2] - 3.0) < 1e-9);    // v3 at its forced floor
    REQUIRE(std::abs(f.fluxes[3] - 17.0) < 1e-9);   // v4
    // Every bound in force is reported, not only the two that were passed in: a flux
    // without the medium it was allowed is not interpretable.
    REQUIRE(f.exchangeBounds.size() == n.reactions.size());
    REQUIRE(f.objectiveReactionId == "v6");

    // At 100% of the optimum the pathway is fully determined; at 50% the lactate
    // branch absorbs the slack and its range must open up.
    const auto pinned = fva(n, "v6", kFermentationBounds, 1.0);
    REQUIRE(pinned.size() == n.reactions.size());
    for (const FluxRange& r : pinned) REQUIRE(std::abs(r.maximum - r.minimum) < 1e-8);
    const auto loose = fva(n, "v6", kFermentationBounds, 0.5);
    REQUIRE(loose[2].maximum - loose[2].minimum > 1.0);

    const FluxSolution p = parsimonious(n, "v6", kFermentationBounds);
    REQUIRE(std::abs(p.objectiveValue - 17.0) < 1e-9);
    double parsimoniousTotal = 0, plainTotal = 0;
    for (double v : p.fluxes) parsimoniousTotal += std::abs(v);
    for (double v : f.fluxes) plainTotal += std::abs(v);
    REQUIRE(parsimoniousTotal <= plainTotal + 1e-9);

    const auto single = deletions(n, "v6", kFermentationBounds, 1);
    const auto objectiveAfter = [&](const std::string& id) {
        for (const FluxRange& r : single)
            if (r.reactionId == id) return r.minimum;
        return -1.0;
    };
    REQUIRE(std::abs(objectiveAfter("v1")) < 1e-9);   // uptake is essential
    REQUIRE(std::abs(objectiveAfter("v4")) < 1e-9);
    // Deleting v3 REMOVES the forced lactate demand, so acetaldehyde rises to 20.
    REQUIRE(std::abs(objectiveAfter("v3") - 20.0) < 1e-9);
    REQUIRE(deletions(n, "v6", kFermentationBounds, 2).size() == 21);
}

TEST_CASE("an unbalanced or unchecked reaction stops FBA before it runs",
          "[network][flux]") {
    NetworkSpec bad = fermentation();
    // AcAld -> Lac is C2H4O -> C3H6O3: a carbon, two hydrogens and two oxygens from
    // nowhere.
    bad.reactions.push_back(
        rx("vbad", {{"AcAld", 1}}, {{"Lac", 1}}, RateLaw::MassAction, {1.0}));
    const FluxSolution check = balance(bad);
    REQUIRE_FALSE(check.massBalanced);
    REQUIRE(check.unbalancedReactions.size() == 1);
    REQUIRE(check.unbalancedReactions[0].find("vbad") != std::string::npos);

    const FluxSolution refused = fba(bad, "v6", kFermentationBounds);
    REQUIRE_FALSE(refused.feasible);
    REQUIRE(refused.fluxes.empty());
    REQUIRE(refused.solverStatus.find("refused") != std::string::npos);

    // A species with no formula is UNCHECKED, and unchecked counts as unbalanced
    // rather than as "probably fine".
    NetworkSpec unchecked = fermentation();
    unchecked.species[4].formula.clear();
    const FluxSolution u = balance(unchecked);
    REQUIRE_FALSE(u.massBalanced);
    REQUIRE(u.unbalancedReactions[0].find("UNCHECKED") != std::string::npos);
}

TEST_CASE("control analysis reproduces hand-computed control coefficients", "[network]") {
    // S <=> X (kf = kr = 1, S clamped at 1), X -> P (k2 = 2, P clamped).
    // Steady state X = kf.S/(kr + k2) = 1/3, J = k2.X = 2/3. By hand:
    //   C^J_1 = k2/(kr+k2) = 2/3, C^J_2 = kr/(kr+k2) = 1/3   (sum to 1)
    //   C^S_1 = 2/3, C^S_2 = -2/3                            (sum to 0)
    //   eps_1,X = -kr.X/v1 = -0.5, eps_2,X = 1
    //   connectivity: (2/3)(-0.5) + (1/3)(1) = 0
    NetworkSpec n;
    n.id = "mca";
    n.species = {sp("S", 1.0, true), sp("X", 0.2), sp("P", 0.0, true)};
    n.reactions = {rx("v1", {{"S", 1}}, {{"X", 1}}, RateLaw::ReversibleMassAction, {1.0, 1.0}),
                   rx("v2", {{"X", 1}}, {{"P", 1}}, RateLaw::MassAction, {2.0})};

    const SteadyState ss = steadyState(n, 40.0);
    REQUIRE(ss.converged);
    REQUIRE(std::abs(ss.concentrations[1] - 1.0 / 3.0) < 1e-10);
    REQUIRE(std::abs(ss.fluxes[1] - 2.0 / 3.0) < 1e-10);

    const ControlAnalysis ca = controlAnalysis(n, 40.0);
    REQUIRE(std::abs(ca.fluxControlCoefficients[1][0] - 2.0 / 3.0) < 1e-6);
    REQUIRE(std::abs(ca.fluxControlCoefficients[1][1] - 1.0 / 3.0) < 1e-6);
    REQUIRE(std::abs(ca.concentrationControlCoefficients[1][0] - 2.0 / 3.0) < 1e-6);
    REQUIRE(std::abs(ca.concentrationControlCoefficients[1][1] + 2.0 / 3.0) < 1e-6);
    REQUIRE(std::abs(ca.elasticities[0][1] + 0.5) < 1e-9);
    REQUIRE(std::abs(ca.elasticities[1][1] - 1.0) < 1e-9);
    // The theorems, which are how a wrong Jacobian is caught.
    REQUIRE(ca.summationResidual < 1e-6);
    REQUIRE(ca.connectivityResidual < 1e-6);

    // A parameter sweep against the exact steady state X = 1/(1 + 2f), J = 2f.X.
    for (const SweepPoint& point : sweep(n, 1, {0.5, 1.0, 2.0, 4.0}, 40.0)) {
        const double x = 1.0 / (1.0 + 2.0 * point.factor);
        REQUIRE(std::abs(point.state.concentrations[1] - x) < 1e-9);
        REQUIRE(std::abs(point.state.fluxes[1] - 2.0 * point.factor * x) < 1e-9);
    }
}
