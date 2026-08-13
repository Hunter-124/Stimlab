// pkpd/PkEngine.h - PK structural models integrated with numeric::rk4Integrate.
//
// The integrator, not a closed form, is the primary path: the Bateman function is
// singular at ka == ke and no closed form exists for Michaelis-Menten elimination.
// The closed forms below are still public because the summary readout and the unit
// tests must share ONE implementation - a test oracle that is a second, private copy
// of the maths proves only that the copy agrees with itself.
//
// SAFETY SCOPE: everything here is an exposure scenario under stated assumptions.
// Nothing in this header returns, or may be turned into, a dose recommendation.
#pragma once

#include "data/Domain.h"

namespace biocad::pkpd {

// One-compartment oral concentration (mg/L) at time t, the Bateman function
//   C = F*D*ka / (V*(ka-ke)) * (exp(-ke*t) - exp(-ka*t)).
// When ka == ke that expression is 0/0; the analytic limit
//   C = F*D*ke*t*exp(-ke*t) / V
// is finite and is what this function returns, so a caller never has to special-case
// the flip-flop boundary.
double batemanConcentration(double dose, double F, double ka, double ke, double V,
                            double t);

// Accumulation ratio for evenly spaced doses: Rac = 1 / (1 - exp(-ke*tau)).
// Returns 0 for a non-positive ke or tau rather than dividing by zero.
double accumulationRatio(double ke, double tauH);

// Average steady-state concentration, Cav,ss = F*D / (CL*tau), in mg/L.
double steadyStateAverage(double F, double doseMg, double clearance, double tauH);

// Integrate the regimen. Dose events are applied BETWEEN integration segments, so an
// instantaneous dose is never smeared across an RK4 step; a zero-order infusion is a
// rate term that is switched on and off exactly at the window boundaries.
PkProfile simulate(const PkModelSpec& spec, const DoseRegimen& regimen);

// Fractional occupancy theta = [A]free / (Kd + [A]free) over the profile's unbound
// series. A missing or non-positive Kd yields an EMPTY curve naming the missing Kd:
// an occupancy curve computed from a guessed Kd is the single most misleading thing
// this subsystem could draw.
OccupancyCurve occupancy(const PkProfile& profile, const Quantity& kd);

}  // namespace biocad::pkpd
