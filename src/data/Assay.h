// data/Assay.h - JSON-serializable DTOs for real experimental data: plates,
// wells, QC, curve fits, biophysics traces and prospective assay design.
//
// WHY these types exist: every other module in BioCAD turns a structure into a
// prediction. This one turns a plate reader export into numbers with error bars.
// That makes the provenance rule sharper here than anywhere else - a raw well
// value is Measured, a fitted parameter is Model, and the two must never be
// rendered in the same colour.
//
// SAFETY SCOPE: this is experimental design and analysis. There is no synthesis
// procedure, no reagent sourcing and no dose: a concentration ladder for a plate
// is not a regimen for a person.
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "data/Domain.h"

namespace biocad {

// ---------------------------------------------------------------------------
// 10.1 - One auditable representation.
// ---------------------------------------------------------------------------

// What a well is FOR. The role drives QC: a plate with no positive and negative
// control cannot have a Z-prime, and saying so is better than computing one from
// whatever the extreme wells happened to be.
enum class WellRole { Unknown, Sample, PositiveControl, NegativeControl, Blank, Reference, Empty };

NLOHMANN_JSON_SERIALIZE_ENUM(WellRole, {
    {WellRole::Unknown, "unknown"},
    {WellRole::Sample, "sample"},
    {WellRole::PositiveControl, "positive"},
    {WellRole::NegativeControl, "negative"},
    {WellRole::Blank, "blank"},
    {WellRole::Reference, "reference"},
    {WellRole::Empty, "empty"},
})

// One measurement. `excluded` never deletes the point: it hollows the marker and
// records which rule hollowed it, so an exclusion is auditable after the fact.
struct Well {
    std::string plateId;
    std::string well;          // "A1", "AF48"
    int         row = 0;       // 0-based, parsed from the letter block
    int         column = 0;    // 0-based
    WellRole    role = WellRole::Unknown;
    std::string sampleId;
    std::string seriesId;      // groups the wells of one concentration series
    double      concentration = 0;
    std::string concUnit;
    int         replicate = 0;
    double      readout = 0;
    std::string readoutUnit;
    double      timeS = 0;
    double      temperatureC = 0;
    bool        excluded = false;
    std::string exclusionRule;  // rule id that excluded it; empty when included
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Well, plateId, well, row, column, role, sampleId, seriesId,
                                   concentration, concUnit, replicate, readout, readoutUnit,
                                   timeS, temperatureC, excluded, exclusionRule)

struct Plate {
    std::string       id;
    int               rows = 0;      // 8, 16, 32
    int               columns = 0;   // 12, 24, 48
    std::vector<Well> wells;
    std::string       readoutUnit;
    std::string       barcode;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Plate, id, rows, columns, wells, readoutUnit, barcode)

// An imported data set. Unknown input columns survive in `metadata` rather than
// being dropped: the instrument knows things about the run that BioCAD does not,
// and losing them makes the import irreversible.
struct AssayDataset {
    std::string                                      id;
    std::vector<Plate>                               plates;
    std::vector<std::pair<std::string, std::string>> metadata;
    std::vector<std::string>                         warnings;
    std::string                                      sourceFile;
    std::string                                      detectedLayout;  // "long" or "96/384/1536 grid"
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(AssayDataset, id, plates, metadata, warnings, sourceFile,
                                   detectedLayout)

// ---------------------------------------------------------------------------
// 10.2 - Plate QC and normalization.
// ---------------------------------------------------------------------------

// A per-plate QC verdict. Bands are the published Z-prime ones (>= 0.5
// excellent, 0 to 0.5 marginal, <= 0 unusable) and are stated in `interpretation`
// rather than being collapsed into a colour.
struct QcReport {
    std::string plateId;
    Quantity    zPrime;         // Zhang et al. 1999; NotComputed without both controls
    Quantity    robustZPrime;   // median/MAD variant
    Quantity    ssmd;
    Quantity    signalToBackground;
    Quantity    signalToNoise;
    Quantity    positiveMean;
    Quantity    positiveSd;
    Quantity    negativeMean;
    Quantity    negativeSd;
    Quantity    cvPositivePct;  // only for positive ratio-scale data
    Quantity    cvNegativePct;
    // Edge and pattern effects are REPORTED, never auto-corrected: silently
    // median-polishing a plate hides the pipetting problem that caused it.
    Quantity    edgeEffectP;    // Mann-Whitney outer ring vs interior
    Quantity    rowEffectP;     // Kruskal-Wallis across rows
    Quantity    columnEffectP;  // Kruskal-Wallis across columns
    std::string interpretation;
    std::vector<std::string> warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(QcReport, plateId, zPrime, robustZPrime, ssmd,
                                   signalToBackground, signalToNoise, positiveMean, positiveSd,
                                   negativeMean, negativeSd, cvPositivePct, cvNegativePct,
                                   edgeEffectP, rowEffectP, columnEffectP, interpretation, warnings)

enum class Normalization { None, PercentOfControl, NormalizedPercentInhibition, ZScore,
                           RobustZScore, BScore };

NLOHMANN_JSON_SERIALIZE_ENUM(Normalization, {
    {Normalization::None, "none"},
    {Normalization::PercentOfControl, "percent of control"},
    {Normalization::NormalizedPercentInhibition, "npi"},
    {Normalization::ZScore, "z"},
    {Normalization::RobustZScore, "robust z"},
    {Normalization::BScore, "b-score"},
})

// ---------------------------------------------------------------------------
// 10.3 - Fitting.
// ---------------------------------------------------------------------------

enum class AssayModel {
    FourParameterLogistic,
    FiveParameterLogistic,
    MichaelisMenten,
    Hill,
    SubstrateInhibition,
    MorrisonTightBinding,
    LangmuirKinetics,      // SPR/BLI global 1:1
    MassTransportKinetics, // two-compartment SPR
    BoltzmannMelt,         // DSF
    TwoStateThermodynamic, // DSF with dCp
    WisemanIsotherm        // ITC one set of sites
};

NLOHMANN_JSON_SERIALIZE_ENUM(AssayModel, {
    {AssayModel::FourParameterLogistic, "4pl"},
    {AssayModel::FiveParameterLogistic, "5pl"},
    {AssayModel::MichaelisMenten, "michaelis-menten"},
    {AssayModel::Hill, "hill"},
    {AssayModel::SubstrateInhibition, "substrate inhibition"},
    {AssayModel::MorrisonTightBinding, "morrison"},
    {AssayModel::LangmuirKinetics, "langmuir 1:1"},
    {AssayModel::MassTransportKinetics, "mass transport"},
    {AssayModel::BoltzmannMelt, "boltzmann"},
    {AssayModel::TwoStateThermodynamic, "two-state"},
    {AssayModel::WisemanIsotherm, "wiseman"},
})

// Inhibition modality lives in data/Domain.h because Phase 4's Cheng-Prusoff
// consumes it; the global fit over the full [S] x [I] matrix below is its ONE
// producer, and answers InhibitionModality::Unknown when the AICc difference
// between modalities is under 2.

struct FittedParameter {
    std::string name;
    Quantity    value;          // Model provenance; error from the covariance
    double      profileLower = 0;   // profile-likelihood CI, when computed
    double      profileUpper = 0;
    bool        profileComputed = false;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(FittedParameter, name, value, profileLower, profileUpper,
                                   profileComputed)

// One fit. `extrapolated` is the flag that stops a curve from lying: an EC50
// outside the tested concentration range is greyed, not reported as a result.
struct FitResult {
    std::string                  seriesId;
    AssayModel                   model = AssayModel::FourParameterLogistic;
    std::vector<FittedParameter> parameters;
    Quantity                     derivedEc50;   // for 5PL this is C*(2^(1/G)-1)^(1/B), never C
    Quantity                     derivedKd;     // kd/ka for kinetics
    double                       rSquared = 0;
    double                       aicc = 0;
    double                       conditionNumber = 0;
    std::size_t                  rank = 0;
    std::size_t                  observations = 0;
    bool                         converged = false;
    bool                         extrapolated = false;   // parameter outside the tested range
    bool                         robust = false;         // Tukey-biweight IRLS was used
    InhibitionModality           modality = InhibitionModality::Unknown;
    std::vector<double>          fittedX;
    std::vector<double>          fittedY;
    std::vector<double>          residuals;
    std::vector<std::string>     assumptions;
    std::vector<std::string>     warnings;
    std::string                  note;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(FitResult, seriesId, model, parameters, derivedEc50, derivedKd,
                                   rSquared, aicc, conditionNumber, rank, observations, converged,
                                   extrapolated, robust, modality, fittedX, fittedY, residuals,
                                   assumptions, warnings, note)

// A model-selection comparison. Reporting the runner-up and the delta is what
// makes "this is competitive inhibition" checkable.
struct ModelComparison {
    std::vector<FitResult> candidates;   // ascending AICc
    double                 deltaAicc = 0;
    bool                   decisive = false;   // false when deltaAicc < 2
    std::string             conclusion;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelComparison, candidates, deltaAicc, decisive, conclusion)

// ---------------------------------------------------------------------------
// 10.4 - Assay design simulation.
// ---------------------------------------------------------------------------

// The truth model and error structure a simulated plate is generated from. Every
// field is the user's stated belief, which is why the design report echoes them
// back: a power calculation is only as good as the noise you assumed.
struct AssayDesignSpec {
    AssayModel          truthModel = AssayModel::FourParameterLogistic;
    std::vector<double> truthParameters;
    std::vector<double> concentrations;      // mol/L, the ladder to be simulated
    int                 replicates = 3;
    int                 rows = 8;
    int                 columns = 12;
    double              additiveNoiseSd = 0;      // readout units
    double              proportionalNoiseCv = 0;  // fraction of signal
    double              pipettingCv = 0;          // compounding lognormal, per transfer
    double              plateGradientPct = 0;     // linear row/column gradient, percent of signal
    double              dmsoTolerancePct = 0;     // signal loss at the top DMSO concentration
    std::uint64_t       seed = 0;
    int                 replicateRuns = 1000;     // seeded Monte Carlo repetitions
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(AssayDesignSpec, truthModel, truthParameters, concentrations,
                                   replicates, rows, columns, additiveNoiseSd,
                                   proportionalNoiseCv, pipettingCv, plateGradientPct,
                                   dmsoTolerancePct, seed, replicateRuns)

// What the simulation actually measured over `replicateRuns` seeded repetitions.
// Empirical CI coverage is the number that matters: a fitter whose 95% interval
// covers the truth 70% of the time is not reporting a 95% interval.
struct AssayDesignReport {
    AssayDesignSpec          spec;
    Quantity                 medianZPrime;
    Quantity                 medianEc50;
    Quantity                 ec50CiWidthLog10;
    Quantity                 empiricalCoveragePct;
    Quantity                 convergenceRatePct;
    std::vector<double>      recoveredEc50;       // one per replicate run
    std::vector<double>      optimalConcentrations;  // Fedorov D-optimal, ladder-constrained
    double                   dOptimalityGain = 0;    // relative to the entered ladder
    std::vector<std::string> assumptions;
    std::vector<std::string> warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(AssayDesignReport, spec, medianZPrime, medianEc50,
                                   ec50CiWidthLog10, empiricalCoveragePct, convergenceRatePct,
                                   recoveredEc50, optimalConcentrations, dOptimalityGain,
                                   assumptions, warnings)

}  // namespace biocad
