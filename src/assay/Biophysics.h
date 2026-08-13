// assay/Biophysics.h - biophysical binding, melting and calorimetry fits.
//
// SPR/BLI kinetics, DSF/nanoDSF melts and ITC isotherms are the three places where
// a chemist's raw instrument trace becomes an affinity, a Tm or a binding enthalpy.
// Every fit here goes through numeric::levenbergMarquardt with an analytic
// Jacobian, and every ODE model through numeric::rk4Integrate - there is one fitter
// and one integrator in BioCAD, not one per module.
//
// All three functions return the shared biocad::FitResult, so a kinetic KD, a Tm and
// a binding enthalpy are rendered by the same panel code and carry the same
// provenance rule: the trace is Measured, every fitted parameter is Model.
//
// The refusals below are the point of the file:
//   - a steady-state KD is withheld when the association phase never approached
//     equilibrium, because the plateau it would be read off does not exist;
//   - a SYPRO Orange trace is truncated at its fluorescence maximum, because the
//     post-peak decay is dye behaviour and fitting it moves Tm;
//   - an ITC fit without a blank heat-of-dilution run does not happen at all.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "data/Assay.h"
#include "data/Domain.h"

namespace biocad::assay {

// ---------------------------------------------------------------------------
// SPR / BLI kinetics
// ---------------------------------------------------------------------------

// One sensorgram: association from t = 0 to `dissociationStartS`, dissociation
// afterwards. Both phases live in one curve because they share ka, kd and Rmax -
// fitting the dissociation alone is what produces the classic "kd from the tail"
// artefact when there is rebinding.
struct KineticCurve {
    double              concentrationM = 0;      // analyte concentration during association
    double              dissociationStartS = 0;  // buffer switch; C = 0 after it
    std::vector<double> timeS;                   // strictly increasing
    std::vector<double> responseRu;
};

// A concentration series analysed as ONE global fit. ka, kd and Rmax are shared
// across curves; only a baseline offset is per-curve, because a global fit is the
// only way a 1:1 model can be falsified - per-curve ka values always "fit".
struct KineticExperiment {
    std::string                seriesId;
    std::vector<KineticCurve>  curves;
    // Rmax expected from the immobilisation level and the analyte/ligand mass ratio,
    // in RU. 0 means the user did not supply it. Reported beside the fitted Rmax
    // because a fitted Rmax far above the theoretical one means the surface is not
    // behaving as a 1:1 monolayer, whatever the residuals look like.
    double                     theoreticalRmaxRu = 0;
};

struct KineticFitOptions {
    // Fraction of the fitted equilibrium response the association phase must reach
    // before a steady-state KD is reported. 0.90 is stated, not tuned: below it the
    // "plateau" is still rising and the Req read off it is biased low.
    double equilibriumFraction = 0.90;
    double maxStepS = 0.5;          // RK4 substep cap between sample times
    bool   fitOffsets = true;       // per-curve baseline offset
    // Starting values. 0 means "estimate from the trace": ka from the initial
    // association slope, kd from the dissociation tail, Rmax from the maximum
    // response, kt from ka*Rmax. An override is for a user who knows better, not a
    // constant the fitter depends on.
    double kaInitial = 0;           // 1/(M s)
    double kdInitial = 0;           // 1/s
    double rmaxInitial = 0;         // RU
    double ktInitial = 0;           // RU/(M s), mass transport only
};

// Global 1:1 Langmuir: dR/dt = ka*C*(Rmax - R) - kd*R, with C = 0 after the buffer
// switch. FitResult::derivedKd is the kinetic KD = kd/ka with its error propagated
// from the fit covariance in log space. The steady-state KD appears as the
// FittedParameter "KD (steady state)" and is NotComputed when the equilibrium test
// above fails, naming the fraction reached.
FitResult fitLangmuirKinetics(const KineticExperiment& experiment,
                              const KineticFitOptions& options = {});

// Two-compartment mass-transport-limited 1:1. The surface compartment is at quasi
// steady state, which is what makes kt an instrument-scaled transport coefficient in
// RU/(M s) - the form the instruments themselves report:
//   dR/dt = kt*(ka*C*(Rmax - R) - kd*R) / (kt + ka*(Rmax - R))
// As kt -> infinity this collapses to the plain Langmuir model, so fitting both and
// comparing is a real test for transport limitation rather than a stylistic choice.
FitResult fitMassTransportKinetics(const KineticExperiment& experiment,
                                   const KineticFitOptions& options = {});

// ---------------------------------------------------------------------------
// DSF / nanoDSF
// ---------------------------------------------------------------------------

struct MeltCurve {
    std::string         seriesId;
    std::vector<double> temperatureC;   // increasing; near-uniform spacing for the derivative
    std::vector<double> signal;         // arbitrary units (RFU, 350/330 ratio)
    // SYPRO Orange (and other extrinsic dyes) lose signal above the unfolding
    // transition as the dye leaves aggregating protein. That decay is dye behaviour,
    // not unfolding, so the trace is truncated at its maximum before fitting.
    bool                syproOrange = false;
};

struct MeltFitOptions {
    // Heat capacity change of unfolding, kcal/(mol K). It is an INPUT, not a fitted
    // parameter: dCp is not identifiable from a single transition, and letting the
    // fitter absorb it into dHm is how a two-state fit produces a confident wrong
    // enthalpy. 0 reduces the two-state model to the van't Hoff form.
    double deltaCpKcalPerMolK = 0.0;
    // Disagreement above this, in degrees C, between the model Tm and the
    // Savitzky-Golay derivative Tm raises a warning: the two agree for a clean
    // two-state transition and diverge when the trace has a second event.
    double derivativeToleranceC = 1.0;
};

// Boltzmann sigmoid with independent sloping native and unfolded baselines:
//   y(T) = (Fn + Sn*T) + [(Fu + Su*T) - (Fn + Sn*T)] / (1 + exp((Tm - T)/a))
// Six parameters. `a` is the transition width in degrees C and is reported as such,
// never converted into an enthalpy - the Boltzmann form has no thermodynamics in it.
FitResult fitBoltzmannMelt(const MeltCurve& curve, const MeltFitOptions& options = {});

// Two-state thermodynamic model with the Gibbs-Helmholtz temperature dependence:
//   dG(T) = dHm*(1 - T/Tm) - dCp*[(Tm - T) + T*ln(T/Tm)]      (T, Tm in kelvin)
//   fu(T) = 1 / (1 + exp(dG/(R*T)))
// with the same sloping baselines. Fits Fn, Sn, Fu, Su, Tm and dHm; dCp is fixed by
// MeltFitOptions and echoed as an assumption.
FitResult fitTwoStateMelt(const MeltCurve& curve, const MeltFitOptions& options = {});

// Tm from the maximum of the Savitzky-Golay first derivative (window 9, order 2),
// computed independently of any model so that it can contradict one. Requires
// near-uniform temperature spacing; a spacing CV above 5% is refused rather than
// smoothed over, because the convolution coefficients assume a uniform grid.
Quantity derivativeTm(const MeltCurve& curve);

// ---------------------------------------------------------------------------
// ITC
// ---------------------------------------------------------------------------

struct ItcInjection {
    double volumeL = 0;    // injected volume
    double heatUcal = 0;   // integrated peak area, microcalories (exothermic is negative)
};

struct ItcExperiment {
    std::string               seriesId;
    double                    cellVolumeL = 0;     // V0, the active cell volume
    double                    macromoleculeM = 0;  // Mt(0), cell
    double                    titrantM = 0;        // X(0), syringe
    double                    temperatureK = 298.15;
    std::vector<ItcInjection> injections;
    // Heat of dilution per injection, microcalories, from a titration of titrant into
    // buffer. REQUIRED: the dilution heat is the same order as the last few binding
    // injections, which is exactly where n and K are determined from, so fitting
    // without it biases stoichiometry. Absent -> converged = false.
    std::vector<double>       blankHeatUcal;
};

// Wiseman one-set-of-sites isotherm with the displaced-volume correction, fitting
// n, K and dH, then reporting dG = -RT ln K and -T dS = dG - dH.
//
// FitResult::parameters[0] is the Wiseman c = n*K*Mt(0), placed first on purpose: c
// decides whether the experiment could have determined K at all, and reading it after
// the parameters is reading it too late. Outside c ~ 1-1000 a warning is raised.
FitResult fitWisemanIsotherm(const ItcExperiment& experiment);

}  // namespace biocad::assay
