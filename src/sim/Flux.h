// sim/Flux.h - constraint-based metabolic flux: mass/charge balance, FBA, FVA,
// parsimonious FBA and deletion scans.
//
// WHICH LP SOLVER SHIPPED, AND WHY. The plan offered HiGHS from vcpkg. HiGHS is not
// present in this development environment and this project is developed on Linux
// while it builds only on Windows, so adding a vcpkg dependency here would mean
// shipping an LP path that has never once been executed by the person who wrote it.
// Instead this file contains a dense two-phase primal simplex, written here and
// exercised by the test suite against hand-solved fluxes. It is Bland's-rule
// pivoted, so it provably terminates on the degenerate problems FBA produces, and it
// is exact on the small networks BioCAD works with. It is NOT a genome-scale solver:
// a dense tableau over a 2000-reaction reconstruction is out of reach, and
// FluxSolution::warnings says so when the problem is large.
//
// THE BALANCE GATE IS NOT ADVISORY. balance() must pass before fba() will return a
// flux: a reaction whose two sides do not have the same elemental composition and
// the same charge can carry flux out of nowhere, and an objective optimised over
// such a reaction is a number about the arithmetic, not about metabolism. Species
// without a stated formula are reported as unchecked, and an unchecked reaction
// counts as unbalanced.
#pragma once

#include <string>
#include <vector>

#include "data/Systems.h"

namespace biocad::sim {

// -------------------------------------------------------------------- the LP
// minimize c'x subject to Ax = b, l <= x <= u. Infinite bounds are represented by
// +/- kInfiniteBound.
inline constexpr double kInfiniteBound = 1e30;

struct LpProblem {
    std::size_t         variables = 0;
    std::size_t         constraints = 0;
    std::vector<double> a;   // row-major constraints x variables
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> lower;
    std::vector<double> upper;
};

struct LpResult {
    bool                feasible = false;
    bool                bounded = true;
    double              objective = 0;
    std::vector<double> x;
    int                 iterations = 0;
    std::string         status;
};

// Two-phase primal simplex. Bounded variables are handled by shifting to the lower
// bound and adding an explicit slack row per finite upper bound, which keeps the
// inner loop a plain standard-form simplex.
LpResult solveLp(const LpProblem& problem);

// ------------------------------------------------------------------ the module
// Element-by-element and charge balance for every reaction, through
// chem::parseFormula. `massBalanced` is the AND over all reactions.
FluxSolution balance(const NetworkSpec& network);

FluxSolution fba(const NetworkSpec& network, const std::string& objectiveReactionId,
                 const std::vector<FluxBound>& bounds);

// Flux variability: the min and max of every reaction while the objective is held at
// `objectiveFraction` of its optimum.
std::vector<FluxRange> fva(const NetworkSpec& network, const std::string& objectiveReactionId,
                           const std::vector<FluxBound>& bounds, double objectiveFraction);

// Parsimonious FBA: among the optima, the one minimising the total absolute flux.
// |v| is linearised with t_j >= v_j, t_j >= -v_j, which is exact for a minimisation.
FluxSolution parsimonious(const NetworkSpec& network, const std::string& objectiveReactionId,
                          const std::vector<FluxBound>& bounds);

// Single (order 1) or double (order 2) reaction deletions. Each FluxRange reports the
// objective WITH the deletion as `minimum` and the undeleted optimum as `maximum`, so
// a reader can see the loss without a second call.
std::vector<FluxRange> deletions(const NetworkSpec& network,
                                 const std::string& objectiveReactionId,
                                 const std::vector<FluxBound>& bounds, int order);

}  // namespace biocad::sim
