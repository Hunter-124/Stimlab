// pkpd/Fits.h - pure pharmacodynamic fitting functions. Part of target biocad_pkpd.
//
// These are free functions with no state and no I/O so they can be tested against
// closed forms directly; RealPharmacodynamics only forwards to them.
//
// SAFETY SCOPE: nothing here returns or formats a dose or a regimen. The outputs are
// potency, affinity, and antagonism parameters under explicitly stated assumptions.
#pragma once

#include <vector>

#include "data/Domain.h"

namespace biocad::pkpd {

// Four-parameter logistic in log-concentration form:
//   E(x) = Bottom + (Top - Bottom) / (1 + 10^(nH * (log10 EC50 - x))),  x = log10([A])
// log10 EC50 is the fitted parameter, not EC50: the parameter is then scale-free and
// the Jacobian column stays well-conditioned across the ~9 decades an assay spans.
//
// `inverseSquareWeighting` applies 1/y^2 weights (constant relative error), skipping
// non-positive responses because their weight is undefined.
//
// Top and Bottom come back Heuristic with an EMPTY unit: the assay response scale is
// arbitrary (RLU, %, mOD), so attaching a physical unit to them would be a lie and
// makeQuantity() would throw. EC50 is mol/L and the slope is dimensionless.
//
// Fewer than 4 points, or any non-positive concentration, yields converged = false
// with every Quantity NotComputed and a note naming the failed precondition.
CurveFit fitFourParameterLogistic(const std::vector<DoseResponsePoint>& points,
                                        bool inverseSquareWeighting = false);

// Cheng-Prusoff Ki from an IC50. The modality decides the equation AND the required
// inputs; a missing input returns NotComputed naming that field. It never falls back
// to the competitive form, because competitive and uncompetitive differ by 10x at
// [S] = 10*Km (and 100x at [S] = 100*Km), and ChEMBL carries no UNCOMPETITIVE action
// type to disambiguate with, so a guess would fabricate orders of magnitude.
//
// When enzymeConc >= 0 and Ki is within ~10x of [E]t the classic Ki is no longer valid
// (it assumes [I] >> [E]t), so the depletion-corrected Ki - [E]t/2 is reported in the
// returned Quantity's `source` beside the classic value.
Quantity kiFromIc50(const ChengPrusoffInput& in);

// Schild regression of log10(DR - 1) on log10([B]). Reports pA2, the slope, and a 95%
// confidence interval for the slope. KB = 10^(-pA2) is only computed when that CI
// includes 1; otherwise the antagonism is not simple competitive and a KB from a
// non-unit-slope Schild plot has no meaning, so kb is NotComputed.
SchildResult schild(const std::vector<SchildPoint>& points);

}  // namespace biocad::pkpd
