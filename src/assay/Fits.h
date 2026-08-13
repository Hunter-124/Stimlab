// assay/Fits.h - real curve fitting for experimental assay data. Part of target
// biocad_assay.
//
// Every fit here goes through numeric::levenbergMarquardt (or its Tukey-biweight
// IRLS wrapper) with an ANALYTIC Jacobian. There is exactly one fitter in the
// tree; nothing in this file differentiates by finite differences, and nothing
// here implements a second optimiser.
//
// Concentration-response models are fitted with the curve's horizontal location
// as log10, matching src/pkpd/Fits.cpp: the parameter is then scale-free and the
// Jacobian column stays conditioned across the nine decades an assay spans.
//
// Response asymptotes come back Heuristic with an EMPTY unit, because a plate
// readout scale (RLU, mOD, %) is arbitrary; makeQuantity() throws on any unit
// there. Concentrations are mol/L and slopes are dimensionless.
//
// SAFETY SCOPE: potency, affinity and kinetic parameters from measured data. No
// dose, no regimen, no synthesis.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "data/Assay.h"
#include "data/Domain.h"

namespace biocad::assay {

// How a fit is run. Weighting, robustness and profile intervals compose: prior
// (replicate) weights multiply the robust weights inside numeric::tukeyBiweight.
struct FitOptions {
    // 1/y^2, the constant-relative-error weighting. Non-positive responses have
    // no defined relative error and are given weight 0 and counted in the note.
    bool   inverseSquareWeighting = false;
    // 1/sd^2 from a per-point replicate SD. Points with sd <= 0 get weight 0.
    // Takes precedence over inverseSquareWeighting when both are set.
    bool   sdWeighting = false;
    // Tukey-biweight IRLS instead of plain least squares. An outlier is
    // downweighted, never deleted: the point stays on the plot and in the record.
    bool   robust = false;
    // Profile-likelihood intervals for the fitted parameters. Off by default
    // because each one costs a full refit per trial value.
    bool   profileLikelihood = false;
    double confidence = 0.95;   // 0.68 / 0.90 / 0.95 / 0.99 only; see numeric::profileLikelihood
};

// One concentration-response observation. Replicates are separate points; that is
// what makes replicate weighting and the residual count honest. `sd` < 0 means
// "no replicate SD available".
struct DosePoint {
    double concentration = 0;   // mol/L, must be > 0 (it is log-transformed)
    double response = 0;        // arbitrary assay units
    double sd = -1;             // replicate SD in response units; < 0 = absent
};

// One initial-rate observation for an enzyme model.
struct KineticPoint {
    double substrate = 0;   // [S], mol/L
    double velocity = 0;    // rate, arbitrary units per time
    double sd = -1;
};

// One cell of the [S] x [I] matrix. The global modality fit needs the full matrix:
// a single [S] cannot distinguish competitive from uncompetitive inhibition at all
// (the two forms coincide at [S] = Km), which is why the fit reports Unknown
// rather than picking one.
struct InhibitionPoint {
    double substrate = 0;    // [S], mol/L
    double inhibitor = 0;    // [I], mol/L; 0 is the uninhibited control row
    double velocity = 0;
    double sd = -1;
};

// ---------------------------------------------------------------------------
// Concentration-response
// ---------------------------------------------------------------------------

// Four-parameter logistic in log10-concentration space, identical in
// parameterisation to pkpd::fitFourParameterLogistic:
//   y = D + (A - D) / (1 + 10^(B * (log10 C - x))),   x = log10([conc])
// Fitted parameters, in order and by FittedParameter::name:
//   "responseAtHighConc" (A), "responseAtLowConc" (D), "log10EC50" (log10 C),
//   "hillSlope" (B).
// The asymptotes are named by their position on the CONCENTRATION axis, not on the
// plot: "top"/"bottom" read inverted on a descending inhibition curve, where the
// high-concentration plateau is the low response.
// B is an EMPIRICAL SLOPE: amplification and receptor reserve bend it, so it is
// not evidence about binding stoichiometry.
//
// FitResult::derivedEc50 is the linear mol/L EC50. profileLower/Upper on
// "log10EC50" are in LOG10 units.
//
// Fewer than 4 points, or any non-positive concentration, returns
// converged = false with a note naming the failed precondition.
FitResult fitFourParameterLogistic(const std::vector<DosePoint>& points,
                                   const FitOptions& options = {});

// Five-parameter (asymmetric) logistic in LINEAR concentration:
//   y = D + (A - D) / [1 + (x / C)^B]^G
// Fitted parameters: "responseAtLowConc" (A), "responseAtHighConc" (D),
// "log10C" (log10 C), "slopeB" (B), "asymmetryG" (G). NOTE the asymptote order is
// the OPPOSITE of the 4PL above, because in this published parameterisation A is
// the response as concentration goes to zero - which is exactly why the names say
// which end they belong to.
//
// WHY derivedEc50 is not C: with G != 1 the curve is asymmetric and C is no
// longer the half-maximal concentration. The half-maximal point solves
// [1 + (x/C)^B]^G = 2, i.e.
//   EC50 = C * (2^(1/G) - 1)^(1/B)
// and that - never C - is what lands in FitResult::derivedEc50. Reporting C as an
// EC50 for an asymmetric curve is a real, common, order-of-magnitude error.
//
// REFUSED with converged = false and an explanatory note below eight DISTINCT
// concentrations: the fifth parameter is only identifiable from a well-sampled
// asymmetry, and a 5PL on seven points fits the noise.
FitResult fitFiveParameterLogistic(const std::vector<DosePoint>& points,
                                   const FitOptions& options = {});

// ---------------------------------------------------------------------------
// Enzyme kinetics
// ---------------------------------------------------------------------------

// v = Vmax * S / (Km + S). Parameters "Vmax", "Km". Fitted NONLINEARLY on the
// untransformed rates; see lineweaverBurk() for why the double-reciprocal form is
// a plot and not a fit.
FitResult fitMichaelisMenten(const std::vector<KineticPoint>& points,
                             const FitOptions& options = {});

// v = Vmax * S^h / (K^h + S^h). Parameters "Vmax", "log10K", "hillCoefficient".
FitResult fitHill(const std::vector<KineticPoint>& points, const FitOptions& options = {});

// v = Vmax * S / (Km + S * (1 + S / Ki)). Parameters "Vmax", "Km", "Ki".
// FitResult::derivedKd carries the substrate concentration of maximal rate,
// sqrt(Km * Ki), which is the diagnostic a chemist reads off the curve.
FitResult fitSubstrateInhibition(const std::vector<KineticPoint>& points,
                                 const FitOptions& options = {});

// Morrison tight-binding inhibition, as the QUADRATIC solution and not the
// classic 1/(1 + I/Ki_app) approximation:
//   v/v0 = 1 - ([E]t + [I] + Ki_app - sqrt(([E]t + [I] + Ki_app)^2 - 4[E]t[I])) / (2[E]t)
// WHY: the classic form assumes [I] >> [E]t, i.e. that free inhibitor equals
// total inhibitor. When Ki_app approaches [E]t that assumption fails badly - at
// [E]t = [I] = Ki_app the quadratic gives 0.618 of control where the
// approximation gives 0.500, a 24% error in the measured effect, and the fitted
// Ki inherits it.
//
// `points` are (inhibitor concentration, velocity); `enzymeTotal` is [E]t in the
// same molar units and is a MEASURED input, not a fitted one. Parameters "v0" and
// "KiApp"; Ki_app = Ki * (1 + [S]/Km) for a competitive tight binder, so
// converting to Ki needs the substrate conditions the caller holds.
FitResult fitMorrisonTightBinding(const std::vector<DosePoint>& points, double enzymeTotal,
                                  const FitOptions& options = {});

// The two fractional-velocity forms, exposed so a panel or a test can show the
// disagreement instead of asserting it.
double morrisonFraction(double enzymeTotal, double inhibitor, double kiApp);
double classicInhibitionFraction(double inhibitor, double kiApp);

// ---------------------------------------------------------------------------
// Global inhibition modality over the full [S] x [I] matrix
// ---------------------------------------------------------------------------

// Fits FOUR candidate models to the SAME data and ranks them by numeric::aicc:
//   competitive     v = Vmax S / (Km (1 + I/Ki) + S)
//   uncompetitive   v = Vmax S / (Km + S (1 + I/Ki))
//   noncompetitive  v = Vmax S / ((Km + S)(1 + I/Ki))
//   mixed           v = Vmax S / (Km (1 + I/Kic) + S (1 + I/Kiu))
//
// This is the ONE producer of InhibitionModality in BioCAD; Phase 4's
// Cheng-Prusoff consumes it. When the top-two AICc difference is < 2 the data do
// not choose, so ModelComparison::decisive = false and the winning
// FitResult::modality is InhibitionModality::Unknown - the candidates and their
// AICc values are still returned, ascending, so the reader can see the tie.
//
// Two models that both fit to double-precision are ranked by their parameter
// count, not by rounding noise: the sum of squares entering aicc() is floored at
// the relative machine-precision level of the data. Otherwise a residual of
// 1e-31 versus 1e-33 would "decide" a modality.
ModelComparison fitInhibitionModality(const std::vector<InhibitionPoint>& matrix,
                                      const FitOptions& options = {});

// ---------------------------------------------------------------------------
// Inverse prediction and diagnostics
// ---------------------------------------------------------------------------

// Concentration that produces `response` on an already-fitted 4PL or 5PL curve.
// The interval comes from the profile-likelihood interval of the curve's location
// parameter (log10EC50 / log10C), shifted onto the requested response: the
// horizontal uncertainty of a logistic is dominated by that parameter. It is
// therefore an approximate interval and says so in `note`. `fit` must carry a
// computed profile interval (FitOptions::profileLikelihood) or the interval is
// reported as not found.
struct InversePrediction {
    double      concentration = 0;   // mol/L
    double      lower = 0;
    double      upper = 0;
    bool        intervalFound = false;
    bool        extrapolated = false;   // outside the tested concentration range
    std::string note;
};
InversePrediction inversePredict(const FitResult& fit, double response,
                                 const std::vector<DosePoint>& points);

// Lineweaver-Burk transform, FOR DIAGNOSTIC PLOTTING ONLY. It is not a fit entry
// point and no function in this file will fit a line to it: the double-reciprocal
// transform inflates the error of the smallest rates without bound, so a
// regression on it returns a biased Vmax and Km. The nonlinear fitters above are
// the only way to get parameters out of BioCAD.
struct LineweaverBurkPoint {
    double inverseSubstrate = 0;   // 1/[S]
    double inverseVelocity = 0;    // 1/v
};
std::vector<LineweaverBurkPoint> lineweaverBurk(const std::vector<KineticPoint>& points);

// ---------------------------------------------------------------------------
// Well adapters (used by the module and the panels)
// ---------------------------------------------------------------------------

// Included wells only; an excluded well keeps its exclusion rule in the dataset
// and is simply not fitted. Replicate SDs are computed per distinct
// concentration and copied onto every point of that concentration, so
// FitOptions::sdWeighting has something to weight with.
std::vector<DosePoint>       doseSeriesFromWells(const std::vector<Well>& wells);
std::vector<KineticPoint>    kineticSeriesFromWells(const std::vector<Well>& wells);
// Matrix convention: Well::concentration is [S] and Well::seriesId parses as the
// numeric [I] for that column (e.g. seriesId "1e-6"). A seriesId that does not
// parse is skipped and reported by the caller-visible fit note.
std::vector<InhibitionPoint> inhibitionMatrixFromWells(const std::vector<Well>& wells);

FitResult       fitSeries(const std::vector<Well>& series, AssayModel model,
                          const FitOptions& options = {});
ModelComparison compareModels(const std::vector<Well>& series,
                              const std::vector<AssayModel>& models,
                              const FitOptions& options = {});
ModelComparison fitInhibitionModalityFromWells(const std::vector<Well>& matrix,
                                               const FitOptions& options = {});

}  // namespace biocad::assay
