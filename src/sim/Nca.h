// sim/Nca.h - noncompartmental analysis of an observed concentration series.
//
// NCA is the one PK analysis that assumes no structural model, which is exactly
// why it has to be strict about the two places it can still go wrong:
//
//   1. The terminal slope. lambda_z is fitted ONLY to points strictly after Tmax,
//      over at least three of them, and the window is chosen by ADJUSTED R-squared
//      so adding a point has to earn its degree of freedom. A half-life from two
//      points, or from a window that includes the distribution phase, is the most
//      common wrong number in a PK report.
//   2. The extrapolation. AUCinf adds Clast_pred/lambda_z beyond the last sample.
//      Above 20% extrapolated, the number is mostly the fit and not the data, so
//      extrapolationUnreliable is set and every quantity that inherits the
//      extrapolation carries a warning naming it - AUCinf, CL, Vz, Vss, MRT.
//
// The AUC rule is linear-up/log-down: the trapezoid is linear while concentrations
// rise and logarithmic while they fall, because a linear trapezoid over a decaying
// exponential systematically OVERestimates the area, and NCA is not permitted a
// systematic bias in its headline number.
#pragma once

#include "data/Population.h"

namespace biocad::sim {

// Analyse one series. Returns a result whose every field is NotComputed rather
// than zero when its prerequisite is missing: fewer than three post-Tmax points
// means no lambda_z, and no lambda_z means no AUCinf, no half-life, no Vz.
NcaResult noncompartmental(const ConcentrationSeries& observed);

}  // namespace biocad::sim
