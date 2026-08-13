// sim/Network.h - one reaction-network representation and its analytic derivatives.
//
// WHY ONE REPRESENTATION. A deterministic time course, a stochastic trajectory, a
// control analysis and a flux balance are four questions about the SAME
// stoichiometric matrix S. Writing S once (species x reaction, dense, row-major)
// and deriving dc/dt = S.v(c) from it means a reaction cannot exist in the ODE and
// be missing from the FBA, which is the failure mode of every codebase that keeps
// a second copy of the network for the "other" solver.
//
// WHY THE ANALYTIC JACOBIAN IS WRITTEN AS AN EXCLUDING PRODUCT. The tempting
// shortcut for a mass-action rate v = k.prod(c_i^a_i) is
//
//     dv/dc_j = v * a_j / c_j
//
// which is algebraically correct and numerically indefensible: the instant a
// species is depleted (c_j = 0) it divides by zero, and the derivative it should
// return there is generally FINITE and NONZERO. This file therefore always
// recomputes the product over the OTHER species and multiplies by
// a_j * c_j^(a_j - 1), which is exact at c_j = 0. The test suite checks the
// analytic Jacobian against a central difference on a network with a depleted
// species precisely to keep that shortcut from creeping back in.
//
// Structural facts come out of the same matrix: conservation laws are the LEFT
// null space (y with y^T.S = 0, so y.c is invariant), and thermodynamic
// (Wegscheider) cycle conditions are the RIGHT null space (k with S.k = 0, so the
// equilibrium constants around the cycle must multiply to 1). A violated
// Wegscheider condition is an error in the PARAMETERS - it means the rate
// constants describe a perpetual motion machine - and is reported as such, never
// silently integrated.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "data/Systems.h"

namespace biocad::sim {

// A compiled network: identifiers resolved to indices, S materialised once.
class Network {
public:
    // Resolves species references, validates rate-law parameter counts and builds S.
    // Returns false and fills `error` (by name) on an unresolvable reference or a
    // rate law given the wrong number of parameters.
    static bool compile(const NetworkSpec& spec, Network& out, std::string* error);

    [[nodiscard]] std::size_t speciesCount() const { return species_.size(); }
    [[nodiscard]] std::size_t reactionCount() const { return reactions_.size(); }
    [[nodiscard]] const std::vector<std::string>& speciesIds() const { return speciesIds_; }
    [[nodiscard]] const std::vector<std::string>& reactionIds() const { return reactionIds_; }
    // Row-major speciesCount() x reactionCount().
    [[nodiscard]] const std::vector<double>& stoichiometry() const { return s_; }
    [[nodiscard]] bool boundary(std::size_t i) const { return species_[i].boundary; }
    [[nodiscard]] std::vector<double> initialState() const;

    // v(c), one entry per reaction.
    void rates(const std::vector<double>& c, std::vector<double>& v) const;
    // dc/dt = S.v(c); boundary species get an exact zero, because a clamped pool is
    // held constant by definition and must not drift by integration error.
    void derivatives(const std::vector<double>& c, std::vector<double>& dcdt) const;
    // dv/dc, row-major reactionCount() x speciesCount(). Excluding-product form.
    void rateJacobian(const std::vector<double>& c, std::vector<double>& dvdc) const;
    // J = S . dv/dc, row-major speciesCount() x speciesCount().
    void jacobian(const std::vector<double>& c, std::vector<double>& j) const;

    // Mutable parameter access, for sweeps, sensitivity and control analysis: the
    // network is recompiled once and only the numbers move.
    [[nodiscard]] const std::vector<double>& parameters(std::size_t reaction) const {
        return reactions_[reaction].parameters;
    }
    void setParameter(std::size_t reaction, std::size_t index, double value) {
        reactions_[reaction].parameters[index] = value;
    }

private:
    struct Term { std::size_t species = 0; double stoichiometry = 0; };
    struct Reaction {
        RateLaw           law = RateLaw::MassAction;
        std::vector<Term> reactants;
        std::vector<Term> products;
        std::vector<double> parameters;
    };
    // Only the fields the maths needs; the DTO keeps the presentation ones.
    struct Species { bool boundary = false; double initial = 0; };

    std::vector<Species>     species_;
    std::vector<Reaction>    reactions_;
    std::vector<std::string> speciesIds_;
    std::vector<std::string> reactionIds_;
    std::vector<double>      s_;
};

// How many parameters a rate law takes, and what they are called. Used by the
// compiler, by the SBML reader and by the panel, so the three cannot disagree.
std::size_t rateLawParameterCount(RateLaw law);
const char* const* rateLawParameterNames(RateLaw law);

// Fills conservationLaws / conservationLabels (left null space of S),
// thermodynamicCycles (right null space) and wegscheiderViolations. Everything
// else in the returned copy is untouched.
NetworkSpec analyze(const NetworkSpec& spec);

}  // namespace biocad::sim
