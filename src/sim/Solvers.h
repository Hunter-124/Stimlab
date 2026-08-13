// sim/Solvers.h - deterministic and stochastic integration of a sim::Network.
//
// THREE SOLVERS, ONE NETWORK, AND A REPORT THAT SAYS WHICH RAN.
//
// * rk4  - numeric::rk4Integrate, fixed step. Correct and cheap for the smooth,
//          well-scaled teaching networks (A -> B, A <=> B, A -> B -> C) whose
//          analytic solutions the tests check against. It is EXPLICIT, so on a
//          stiff problem its stability limit collapses: the Robertson problem needs
//          a step around 1e-6 for a 40-second horizon, i.e. tens of millions of
//          steps, and taking a sane step instead produces a smooth, plausible,
//          catastrophically wrong curve. The test suite demonstrates that rather
//          than asserting it in a comment.
// * rosenbrock - ROS3P (Lang & Verwer), a 3-stage third-order L-stable Rosenbrock
//          method with an embedded second-order solution for step control. ONE LU
//          factorization of (I - gamma.h.J) per step is shared by all three stages;
//          that is the whole reason a Rosenbrock method is used here rather than a
//          Newton-iterated BDF. L-stability (not merely A-stability) is what keeps
//          the fast Robertson transient from ringing.
// * gillespie / tau-leap - the exact stochastic simulation algorithm and Cao-
//          Gillespie-Petzold explicit tau-leaping, over the SAME network and the
//          same seeded PCG64-DXSM stream the population layer uses. At low copy
//          number the variance IS the result, and a deterministic mean hides it.
//
// EVERY result fills a SolverReport: method, tolerances, accepted and rejected
// steps, Jacobian evaluations, non-negativity clips and CPU seconds. A nonzero clip
// count is a WARNING surfaced in TimeCourse::warnings, never a silent repair - a
// solver that quietly clamps a negative concentration is hiding a step-size failure.
//
// STOCHASTIC STATE IS COPY NUMBER, NOT CONCENTRATION. The SSA counts molecules, so
// SpeciesSpec::initialConcentration is read as a copy number in a unit volume and
// mass-action rate constants are read as stochastic propensity constants with the
// combinatorial factor prod(n_i choose a_i) * a_i! applied. Saturable rate laws are
// NOT elementary reaction channels, so a network containing one is refused by name
// rather than approximated.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "data/Systems.h"
#include "sim/Network.h"
#include "sim/Random.h"

namespace biocad::sim {

struct IntegrationOptions {
    std::string method = "rosenbrock";   // "rk4" | "rosenbrock"
    double      horizon = 10.0;
    double      relativeTolerance = 1e-8;
    double      absoluteTolerance = 1e-10;
    // Output grid: `outputPoints` equally spaced samples including both ends.
    int         outputPoints = 201;
    // rk4 only: the fixed step. Zero means horizon / (outputPoints - 1) / 10.
    double      fixedStep = 0.0;
    double      initialStep = 1e-6;      // rosenbrock only
    double      maxStep = 0.0;           // 0 => horizon
    std::int64_t maxSteps = 5'000'000;
};

// Integrates `spec` and fills the conserved-quantity audit from the conservation
// laws already in the spec (run sim::analyze first, or they are simply absent).
TimeCourse integrate(const NetworkSpec& spec, const IntegrationOptions& options);

struct StochasticOptions {
    double        horizon = 10.0;
    int           replicates = 100;
    int           outputPoints = 101;
    std::uint64_t seed = 1;
    bool          tauLeap = false;
    double        epsilon = 0.03;        // tau-leap error control (Cao et al. 2006)
    std::int64_t  maxEvents = 20'000'000;
};

// Refuses (through SolverReport::note and an empty ensemble) a network whose rate
// laws are not elementary mass action, naming the offending reaction.
StochasticEnsemble stochastic(const NetworkSpec& spec, const StochasticOptions& options);

// ------------------------------------------------------------------ internals
// Exposed because they are separately testable: a Poisson sampler whose mean and
// variance are wrong would corrupt every tau-leap result silently.

// One adaptive ROS3P step from (t, y) with step h. Returns the error norm scaled by
// the tolerances (accept when <= 1). `report` accumulates Jacobian evaluations.
double rosenbrockStep(const Network& net, std::vector<double>& y, double h, double relTol,
                      double absTol, SolverReport& report);

// Poisson variate with mean `lambda` from the shared stream: Knuth's product method
// below 30 and Hormann's PTRS transformed rejection above it, which is where the
// product method's expected iteration count stops being acceptable.
std::uint64_t poissonVariate(Pcg64Dxsm& rng, double lambda);

}  // namespace biocad::sim
