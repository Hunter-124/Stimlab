// data/Population.h - JSON-serializable DTOs for population PK, uncertainty
// propagation, noncompartmental analysis and mechanistic drug interactions.
//
// SAFETY SCOPE, permanent: everything here is an EXPOSURE SCENARIO. There is no
// dose, no dose adjustment, no dosing-interval recommendation and no per-patient
// prediction anywhere in these types, and no field in which to hide one. A
// percentile band describes the variability that was ENTERED; it is not a
// prediction about an individual, which is why every result carries the Omega,
// the sampler and the seed that produced it.
//
// Whole-body PBPK is deliberately absent. Its required fu, blood-to-plasma
// ratio, tissue partition, Papp and transporter inputs are not derivable from
// anything BioCAD has, so a PBPK profile here would be precise fiction.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "data/Domain.h"

namespace biocad {

// ---------------------------------------------------------------------------
// 13.1 - The reproducible population layer.
// ---------------------------------------------------------------------------

// The three variability layers are separately toggled because they answer
// different questions: BSV asks "how much do people differ", parameter
// uncertainty asks "how well do I know the typical value", and residual error
// asks "how noisy is the assay". Collapsing them into one band makes all three
// unanswerable.
struct VariabilitySpec {
    bool                             betweenSubject = false;
    bool                             parameterUncertainty = false;
    bool                             residualError = false;
    // Omega: the BSV covariance of the log-scale random effects, row-major and
    // square over the varying parameters named in `parameters`.
    std::vector<std::string>         parameters;
    std::vector<double>              omega;
    // The fit covariance for the uncertainty layer, row-major over `parameters`.
    std::vector<double>              parameterCovariance;
    double                           proportionalResidualCv = 0;
    double                           additiveResidualSd = 0;   // mg/L
    std::uint64_t                    seed = 0;
    int                              subjects = 500;
    std::string                      sampler;   // "monte-carlo" | "latin-hypercube"
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(VariabilitySpec, betweenSubject, parameterUncertainty,
                                   residualError, parameters, omega, parameterCovariance,
                                   proportionalResidualCv, additiveResidualSd, seed, subjects,
                                   sampler)

// One time point's percentile band across the simulated population.
struct PercentileBand {
    double timeH = 0;
    double p5 = 0;
    double p50 = 0;
    double p95 = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PercentileBand, timeH, p5, p50, p95)

// The result of a population simulation. `provenanceStatement` is a required
// sentence, not decoration: it is what stops a band from being read as a
// prediction about the reader.
struct PopulationProfile {
    VariabilitySpec              spec;
    std::vector<PercentileBand>  bands;
    // Up to a stated maximum of individual trajectories, for the faint-line
    // overlay. Storing all of them would be 32 MB for 1000 x 4000 points.
    std::vector<std::vector<double>> sampleTrajectories;
    std::vector<double>          times;
    Quantity                     medianAuc;      // mg*h/L
    Quantity                     medianCmax;     // mg/L
    Quantity                     aucCv;          // percent, across subjects
    std::string                  provenanceStatement;
    std::vector<std::string>     assumptions;
    std::vector<std::string>     warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PopulationProfile, spec, bands, sampleTrajectories, times,
                                   medianAuc, medianCmax, aucCv, provenanceStatement,
                                   assumptions, warnings)

// ---------------------------------------------------------------------------
// 13.2 - Noncompartmental analysis.
// ---------------------------------------------------------------------------

// One observed concentration-time series. This is DATA, so every value is
// Measured; NCA turns it into parameters without ever assuming a compartment.
struct ConcentrationSeries {
    std::string         subjectId;
    std::vector<double> timeH;
    std::vector<double> concentration;   // mg/L
    double              dose = 0;        // mg
    bool                intravenous = false;
    double              tauH = 0;        // dosing interval; 0 = single dose
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ConcentrationSeries, subjectId, timeH, concentration, dose,
                                   intravenous, tauH)

// `lambdaZPointCount` and `adjustedRSquared` are reported because a terminal
// slope fitted through two points, or chosen to maximise R-squared over an
// arbitrary window, is the single most common way an NCA half-life goes wrong.
// Anything derived from an extrapolation above 20% is flagged unreliable.
struct NcaResult {
    std::string subjectId;
    Quantity    cmax;                 // mg/L
    Quantity    tmax;                 // h
    Quantity    aucLast;              // mg*h/L, linear-up/log-down
    Quantity    aucInfinity;
    Quantity    percentExtrapolated;  // percent of AUCinf beyond the last point
    Quantity    lambdaZ;              // 1/h
    Quantity    halfLife;             // h
    Quantity    clearance;            // L/h; CL for IV, CL/F otherwise
    Quantity    volumeZ;              // L; Vz or Vz/F
    Quantity    volumeSteadyState;    // L; IV only, NotComputed otherwise
    Quantity    aumc;                 // mg*h^2/L
    Quantity    meanResidenceTime;    // h
    Quantity    aucTau;               // mg*h/L, steady state only
    Quantity    cAverage;             // mg/L
    Quantity    swing;                // dimensionless
    int         lambdaZPointCount = 0;
    double      adjustedRSquared = 0;
    bool        extrapolationUnreliable = false;   // percentExtrapolated > 20
    std::vector<std::string> warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(NcaResult, subjectId, cmax, tmax, aucLast, aucInfinity,
                                   percentExtrapolated, lambdaZ, halfLife, clearance, volumeZ,
                                   volumeSteadyState, aumc, meanResidenceTime, aucTau, cAverage,
                                   swing, lambdaZPointCount, adjustedRSquared,
                                   extrapolationUnreliable, warnings)

// ---------------------------------------------------------------------------
// 13.3 - Mechanistic drug interactions.
// ---------------------------------------------------------------------------

// One perpetrator's in vitro parameters against one enzyme. Every field is an
// INPUT with a citation; nothing here is predicted from structure.
struct PerpetratorSpec {
    std::string label;
    std::string enzyme;          // "CYP3A4"
    double      ki = -1;         // uM, reversible inhibition; < 0 = absent
    double      kinact = -1;     // 1/h, time-dependent inactivation
    double      kI = -1;         // uM, TDI concentration for half-maximal kinact
    double      indEmax = -1;    // fold induction maximum
    double      indEc50 = -1;    // uM
    double      indD = 1.0;      // in vitro-to-in vivo scaling factor for induction
    double      unboundHepaticInletUM = -1;   // [I]h,inlet,u
    double      unboundSystemicUM = -1;       // [I]sys,u
    double      enterocyteUM = -1;            // [I]g
    std::string source;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PerpetratorSpec, label, enzyme, ki, kinact, kI, indEmax,
                                   indEc50, indD, unboundHepaticInletUM, unboundSystemicUM,
                                   enterocyteUM, source)

// The victim's disposition. fm and Fg are inputs and have no default: assuming
// fm = 1 is what turns a 1.3-fold interaction into a 5-fold one.
struct VictimSpec {
    std::string label;
    double      fractionMetabolizedByEnzyme = -1;   // fm; < 0 = absent
    double      intestinalAvailability = -1;        // Fg; < 0 = absent
    double      fractionExcretedUnchanged = -1;     // fe, for the renal scenario
    std::string source;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(VictimSpec, label, fractionMetabolizedByEnzyme,
                                   intestinalAvailability, fractionExcretedUnchanged, source)

// FDA basic-model screening R-values plus the mechanistic static AUCR. The
// theoretical ceiling (1/(1-fm)) is displayed beside the AUCR because it is the
// number that says whether a predicted ratio is even attainable.
struct InteractionReport {
    PerpetratorSpec perpetrator;
    VictimSpec      victim;
    Quantity        r1;              // reversible, hepatic: 1 + [I]u/Ki
    Quantity        r1Gut;           // 1 + [I]g/Ki
    Quantity        r2;              // TDI
    Quantity        rInduction;      // induction
    Quantity        aucRatio;        // mechanistic static, all mechanisms combined
    Quantity        aucRatioHepaticOnly;
    Quantity        theoreticalCeiling;  // 1/(1 - fm)
    std::string     dominantMechanism;
    bool            gutIncluded = false;
    std::vector<std::string> assumptions;
    std::vector<std::string> warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(InteractionReport, perpetrator, victim, r1, r1Gut, r2,
                                   rInduction, aucRatio, aucRatioHepaticOnly,
                                   theoreticalCeiling, dominantMechanism, gutIncluded,
                                   assumptions, warnings)

// A dynamic enzyme-activity simulation:
//   dE/dt = kdeg*(1 + d*Emax*I/(EC50+I)) - kdeg*E - kinact*I/(KI+I)*E
//   CLint(t) = CLint0 * E / (1 + I/Ki)
// At constant I the steady state of this must EQUAL the static model exactly,
// which is the assertion that keeps the two implementations honest.
struct EnzymeTimeCourse {
    std::vector<double> timeH;
    std::vector<double> relativeActivity;   // E, 1.0 = uninhibited baseline
    std::vector<double> clearanceRatio;     // CLint(t)/CLint0
    double              steadyStateActivity = 0;
    double              staticModelActivity = 0;   // the static equivalent
    double              agreement = 0;              // |dynamic - static|
    double              kdegUsed = 0;
    std::string         kdegSource;
    std::vector<std::string> assumptions;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(EnzymeTimeCourse, timeH, relativeActivity, clearanceRatio,
                                   steadyStateActivity, staticModelActivity, agreement,
                                   kdegUsed, kdegSource, assumptions)

// An organ-impairment EXPOSURE SCENARIO. Renal is an explicit function ratio from
// fe; hepatic is an editable well-stirred scenario. Neither is a Child-Pugh or
// creatinine-clearance-to-dose formula, because that is a dosing decision and
// BioCAD does not make those.
struct ImpairmentScenario {
    std::string label;
    double      renalFunctionRatio = 1.0;    // 1.0 = normal
    double      hepaticClintRatio = 1.0;
    Quantity    exposureRatio;               // AUC impaired / AUC normal
    std::vector<std::string> assumptions;
    std::string boundaryStatement;           // the explicit "this is not a dose" sentence
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ImpairmentScenario, label, renalFunctionRatio,
                                   hepaticClintRatio, exposureRatio, assumptions,
                                   boundaryStatement)

}  // namespace biocad
