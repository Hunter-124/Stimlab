// sim/Control.h - steady states, parameter sweeps and metabolic control analysis.
//
// WHY THE SUMMATION AND CONNECTIVITY RESIDUALS ARE PART OF THE RESULT. Flux control
// coefficients must sum to 1 across reactions and concentration control coefficients
// to 0; and for every species, sum_k C^J_k * eps_k,j must be 0 (the connectivity
// theorem). Those are theorems, not fit statistics, so a nonzero residual means the
// steady state was not reached, the Jacobian is wrong, or the perturbation was too
// large - and reporting the residual is how that gets caught instead of shipping a
// confident table of wrong control coefficients.
//
// The coefficients are computed by the OPERATIONAL definition: each reaction's rate is
// scaled by a small factor (which for a reversible step scales BOTH directions, so the
// equilibrium constant is untouched and the perturbation is an enzyme-amount change
// rather than a thermodynamic one), the steady state is recomputed, and the log-log
// derivative is taken by central difference.
#pragma once

#include <string>
#include <vector>

#include "data/Systems.h"
#include "sim/Network.h"

namespace biocad::sim {

struct SteadyState {
    std::vector<double> concentrations;
    std::vector<double> fluxes;
    double              residual = 0;   // max |dc/dt| at the returned point
    bool                converged = false;
    std::string         note;
};

// Integrates to `horizon` with the Rosenbrock solver, then polishes with a
// least-squares Newton step (the Jacobian is singular by exactly the number of
// conservation laws, so a plain solve would fail).
SteadyState steadyState(const NetworkSpec& spec, double horizon);

// Scales reaction `reaction`'s rate by each factor in `factors` and reports the
// resulting steady state. This is the parameter sweep and the raw material of the
// control analysis.
struct SweepPoint {
    double              factor = 1.0;
    SteadyState         state;
};
std::vector<SweepPoint> sweep(const NetworkSpec& spec, std::size_t reaction,
                              const std::vector<double>& factors, double horizon);

ControlAnalysis controlAnalysis(const NetworkSpec& spec, double horizon,
                               double perturbation = 1e-4);

}  // namespace biocad::sim
