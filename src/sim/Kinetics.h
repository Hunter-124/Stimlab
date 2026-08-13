// sim/Kinetics.h - chemical kinetics fitted to data the user actually has.
//
// WHY THIS FILE EXISTS AT ALL. BioCAD used to answer "how long does this keep?" by
// counting perceived functional groups: an ester subtracted a fixed number from a
// 0-100 score, the score picked one of four hard-coded strings ("~24 months @
// 25C/60%RH"), and a pH window came from adding invented interval bounds per group.
// Those three functions - predictPhWindow, predictThermalWindow and the shelf-life
// bucket - are deleted, not softened. A shelf life is an extrapolation of a MEASURED
// degradation rate; there is no structure-only predictor of one, and a number that
// looks like months but is really a flag count is exactly the fabrication this phase
// exists to remove.
//
// What replaces them: Arrhenius and Eyring fits to user-entered rate constants at
// several temperatures, with a covariance matrix, a joint confidence ellipse for the
// activation parameters, and a prediction interval on the extrapolated rate. Fewer
// than three distinct temperatures cannot support an extrapolation - two points fit
// a straight line exactly and say nothing about how wrong it is - so the fit comes
// back with extrapolationSupported = false and the dependent quantities read
// notComputed.
//
// All fitting goes through numeric::levenbergMarquardt with an analytic Jacobian.
// There is one fitter in this application; adding a second here would produce two
// slightly different convergence behaviours for the same job.
#pragma once

#include <string>
#include <vector>

#include "data/Systems.h"

namespace biocad::sim {

// Molar gas constant in kJ/(mol*K), so activation energies come out in kJ/mol
// (CODATA 2018: R = 8.314462618 J/(mol*K), exact by the 2019 SI redefinition).
inline constexpr double kGasConstantKjPerMolK = 8.314462618e-3;
inline constexpr double kBoltzmannJPerK = 1.380649e-23;      // exact, SI 2019
inline constexpr double kPlanckJs = 6.62607015e-34;          // exact, SI 2019

// k = A * exp(-Ea / (R T)). Fitted in (ln A, Ea) so that A, which spans many orders
// of magnitude, is a well-scaled parameter. `rateConstants` must be positive.
KineticsFit arrhenius(const std::vector<double>& temperaturesK,
                      const std::vector<double>& rateConstants);

// k = kappa * (kB T / h) * exp(dS/R) * exp(-dH/(RT)), fitted in (dH kJ/mol,
// dS J/(mol*K)). `transmission` is kappa; it is 1 unless the user has a reason,
// and the value used is recorded in KineticsFit::assumptions because kappa scales
// dS directly.
//
// The dH/dS confidence ellipse is filled from the 2x2 covariance block: the two
// parameters are strongly correlated (the classic isokinetic artefact), so separate
// error bars overstate what the experiment determined.
KineticsFit eyring(const std::vector<double>& temperaturesK,
                   const std::vector<double>& rateConstants, double transmission = 1.0);

// Time to lose `fractionLost` of a first-order-degrading substance at
// `storageTemperatureK`, in months, extrapolated from an Arrhenius fit:
// t = -ln(1 - fractionLost) / k(T). Returns notComputed when the fit does not
// support extrapolation, and carries the fit's own prediction interval as its error.
Quantity shelfLife(const KineticsFit& fit, double storageTemperatureK, double fractionLost);

// k_obs = kH[H+] + k0 + kOH[OH-], fitted in (kH, k0, kOH). The minimum is reported
// from its closed form pH = 0.5*(pKw + log10(kH/kOH)) with minimum rate
// k0 + 2*sqrt(kH*kOH*Kw), not by scanning a grid - a scanned minimum inherits the
// grid spacing as a fake precision.
//
// `pKw` defaults to 14.0 (water at 25 C, ionic strength zero); it is an ASSUMPTION
// recorded in the result, because at 37 C pKw is 13.62 and the reported minimum
// would shift by 0.19 pH units.
PhRateProfile phRate(const std::vector<double>& pHValues,
                     const std::vector<double>& rateConstants, double pKw = 14.0);

}  // namespace biocad::sim
