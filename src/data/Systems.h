// data/Systems.h - JSON-serializable DTOs for reaction-network simulation,
// chemical kinetics, metabolic flux and pathway enrichment.
//
// WHY the solver settings are embedded in every result: a time course whose
// tolerance, step count and rejected-step count are hidden cannot be checked, and
// a stiff network integrated with the wrong solver produces a smooth, plausible,
// wrong curve. TimeCourse therefore carries its own SolverReport, and no DTO here
// has a defaulted objective or background - those are required fields, because a
// flux without its bounds and an enrichment without its background are both
// meaningless.
//
// SAFETY SCOPE: reaction networks describe kinetics the user entered. A computed
// growth rate is a property of the model, never a measurement of an organism, and
// no docking score may ever propagate into a flux, a pathway or a network number.
#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "data/Domain.h"

namespace biocad {

enum class RateLaw { MassAction, ReversibleMassAction, MichaelisMenten, ReversibleMichaelisMenten,
                     Hill };

NLOHMANN_JSON_SERIALIZE_ENUM(RateLaw, {
    {RateLaw::MassAction, "mass action"},
    {RateLaw::ReversibleMassAction, "reversible mass action"},
    {RateLaw::MichaelisMenten, "michaelis-menten"},
    {RateLaw::ReversibleMichaelisMenten, "reversible michaelis-menten"},
    {RateLaw::Hill, "hill"},
})

struct SpeciesSpec {
    std::string id;
    std::string name;
    double      initialConcentration = 0;   // mol/L or copy number for the SSA
    std::string compartment;
    bool        boundary = false;   // held constant (a clamped pool), not integrated
    // Elemental composition and net charge, for the flux module's mass and charge
    // balance. Empty means "not stated", and IFluxModule::balance() reports that as
    // an unbalanced reaction rather than assuming the reaction balances - an FBA over
    // reactions nobody checked is the classic way a model produces mass from nothing.
    std::string formula;
    int         charge = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SpeciesSpec, id, name, initialConcentration, compartment,
                                   boundary, formula, charge)

struct ReactionSpec {
    std::string                      id;
    std::vector<std::pair<std::string, double>> reactants;   // species id, stoichiometry
    std::vector<std::pair<std::string, double>> products;
    std::vector<std::string>         modifiers;
    RateLaw                          law = RateLaw::MassAction;
    std::vector<double>              parameters;   // k, or Vmax/Km, or Vmax/Km/n
    std::vector<std::string>         parameterNames;
    bool                             reversible = false;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ReactionSpec, id, reactants, products, modifiers, law,
                                   parameters, parameterNames, reversible)

// A network plus the structural facts derived from its stoichiometric matrix.
// Conservation laws come from the LEFT null space and thermodynamic (Wegscheider)
// cycle constraints from the RIGHT null space; both are displayed, and a violated
// Wegscheider condition is an error in the parameters, not a rounding issue.
struct NetworkSpec {
    std::string                       id;
    std::vector<SpeciesSpec>          species;
    std::vector<ReactionSpec>         reactions;
    std::vector<std::vector<double>>  conservationLaws;      // left null space basis
    std::vector<std::string>          conservationLabels;
    std::vector<std::vector<double>>  thermodynamicCycles;   // right null space basis
    std::vector<std::string>          wegscheiderViolations;
    std::vector<std::string>          warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(NetworkSpec, id, species, reactions, conservationLaws,
                                   conservationLabels, thermodynamicCycles,
                                   wegscheiderViolations, warnings)

// What the solver actually did. A hidden solver is not a trustworthy solver.
struct SolverReport {
    std::string   method;          // "rk4" | "rosenbrock" | "gillespie" | "tau-leap"
    double        relativeTolerance = 0;
    double        absoluteTolerance = 0;
    std::int64_t  acceptedSteps = 0;
    std::int64_t  rejectedSteps = 0;
    std::int64_t  jacobianEvaluations = 0;
    std::int64_t  nonNegativityClips = 0;   // a nonzero count is a warning, not a fix
    double        cpuSeconds = 0;
    std::uint64_t seed = 0;                 // stochastic methods only
    std::string   note;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SolverReport, method, relativeTolerance, absoluteTolerance,
                                   acceptedSteps, rejectedSteps, jacobianEvaluations,
                                   nonNegativityClips, cpuSeconds, seed, note)

struct TimeCourse {
    std::vector<double>              times;
    // Column-major by species: `trajectories[i]` is one species over all times.
    std::vector<std::vector<double>> trajectories;
    std::vector<std::string>         speciesIds;
    // Each conservation law's value over time. It must stay constant; the worst
    // drift is the number that says whether the integration can be trusted.
    std::vector<std::vector<double>> conservedQuantities;
    double                           worstConservationDrift = 0;
    SolverReport                     solver;
    std::vector<std::string>         warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(TimeCourse, times, trajectories, speciesIds,
                                   conservedQuantities, worstConservationDrift, solver, warnings)

// An ensemble of stochastic trajectories, summarised. Mean and variance are
// reported because for the models where the SSA matters, the variance IS the
// result - a deterministic mean would hide it.
struct StochasticEnsemble {
    std::vector<double>              times;
    std::vector<std::vector<double>> mean;
    std::vector<std::vector<double>> variance;
    std::vector<std::string>         speciesIds;
    int                              replicates = 0;
    SolverReport                     solver;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(StochasticEnsemble, times, mean, variance, speciesIds,
                                   replicates, solver)

// ---------------------------------------------------------------------------
// Chemical kinetics: Arrhenius, Eyring, and stability fits to REAL data.
// ---------------------------------------------------------------------------

// An Arrhenius/Eyring fit to user-entered rate data. Fewer than three
// temperatures cannot support an extrapolation, and the DTO says so rather than
// extrapolating from two.
struct KineticsFit {
    std::vector<double> temperaturesK;
    std::vector<double> rateConstants;
    Quantity            preExponential;      // A, same units as k
    Quantity            activationEnergy;    // Ea, kJ/mol
    Quantity            enthalpyOfActivation;  // dH*, kJ/mol
    Quantity            entropyOfActivation;   // dS*, J/(mol*K)
    Quantity            predictedRateAt25C;
    double              predictionIntervalLow = 0;
    double              predictionIntervalHigh = 0;
    // The dH*/dS* joint confidence ellipse, as (dH, dS) boundary points. The two
    // are strongly correlated, so separate error bars overstate what is known.
    std::vector<std::pair<double, double>> confidenceEllipse;
    double              rSquared = 0;
    bool                extrapolationSupported = false;   // >= 3 temperatures
    std::vector<std::string> assumptions;
    std::vector<std::string> warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(KineticsFit, temperaturesK, rateConstants, preExponential,
                                   activationEnergy, enthalpyOfActivation, entropyOfActivation,
                                   predictedRateAt25C, predictionIntervalLow,
                                   predictionIntervalHigh, confidenceEllipse, rSquared,
                                   extrapolationSupported, assumptions, warnings)

// pH-rate profile k_obs = kH[H+] + k0 + kOH[OH-]. The minimum has the closed form
// pH = 0.5*(pKw + log10(kH/kOH)), which is reported rather than found by scanning.
struct PhRateProfile {
    std::vector<double> pHValues;
    std::vector<double> rateConstants;
    Quantity            kAcid;
    Quantity            kNeutral;
    Quantity            kBase;
    Quantity            minimumPh;
    Quantity            minimumRate;
    double              rSquared = 0;
    std::vector<std::string> warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PhRateProfile, pHValues, rateConstants, kAcid, kNeutral,
                                   kBase, minimumPh, minimumRate, rSquared, warnings)

// ---------------------------------------------------------------------------
// Metabolic flux.
// ---------------------------------------------------------------------------

struct FluxBound {
    std::string reactionId;
    double      lower = 0;
    double      upper = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(FluxBound, reactionId, lower, upper)

// A flux solution. `objective` and `exchangeBounds` are REQUIRED: a flux
// distribution without the objective it optimised and the medium it was allowed
// is not interpretable, so no renderer may display one without them.
struct FluxSolution {
    std::string              objectiveReactionId;
    double                   objectiveValue = 0;
    std::vector<double>      fluxes;
    std::vector<std::string> reactionIds;
    std::vector<FluxBound>   exchangeBounds;
    bool                     massBalanced = false;   // formula/charge balance passed
    std::vector<std::string> unbalancedReactions;
    bool                     feasible = false;
    std::string              solverStatus;
    std::vector<std::string> warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(FluxSolution, objectiveReactionId, objectiveValue, fluxes,
                                   reactionIds, exchangeBounds, massBalanced,
                                   unbalancedReactions, feasible, solverStatus, warnings)

struct FluxRange {
    std::string reactionId;
    double      minimum = 0;
    double      maximum = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(FluxRange, reactionId, minimum, maximum)

// ---------------------------------------------------------------------------
// Sensitivity, control analysis, enrichment and networks.
// ---------------------------------------------------------------------------

// Metabolic control analysis with its summation and connectivity checks attached:
// flux control coefficients must sum to 1 and concentration control coefficients
// to 0, and reporting the residual is how a wrong Jacobian is caught.
struct ControlAnalysis {
    std::vector<std::string>         reactionIds;
    std::vector<std::string>         speciesIds;
    std::vector<std::vector<double>> fluxControlCoefficients;
    std::vector<std::vector<double>> concentrationControlCoefficients;
    std::vector<std::vector<double>> elasticities;
    double                           summationResidual = 0;
    double                           connectivityResidual = 0;
    std::vector<std::string>         warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ControlAnalysis, reactionIds, speciesIds,
                                   fluxControlCoefficients, concentrationControlCoefficients,
                                   elasticities, summationResidual, connectivityResidual,
                                   warnings)

// One over-represented pathway. The q-value is Benjamini-Hochberg over the whole
// tested set; a raw p-value alone from a few hundred pathway tests is not a
// result.
struct EnrichmentHit {
    std::string pathwayId;
    std::string pathwayName;
    int         inSetAndPathway = 0;
    int         pathwaySize = 0;
    double      pValue = 0;
    double      qValue = 0;
    double      foldEnrichment = 0;
    std::string source;        // database name with its release date
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(EnrichmentHit, pathwayId, pathwayName, inSetAndPathway,
                                   pathwaySize, pValue, qValue, foldEnrichment, source)

// `background` is a required field with no default: the hypergeometric test's
// answer is a function of the background, and silently using "all annotated
// genes" is the most common way an enrichment result becomes wrong.
struct EnrichmentReport {
    std::vector<std::string>   querySet;
    std::vector<std::string>   background;
    std::vector<EnrichmentHit> hits;        // ascending q
    std::string                databaseRelease;
    std::string                method;      // "hypergeometric + BH"
    std::vector<std::string>   warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(EnrichmentReport, querySet, background, hits,
                                   databaseRelease, method, warnings)

struct NetworkEdge {
    std::string source;
    std::string target;
    double      weight = 0;
    std::string evidence;      // required: how this edge is known
    Provenance  provenance = Provenance::NotComputed;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(NetworkEdge, source, target, weight, evidence, provenance)

struct GraphMetrics {
    std::vector<std::string>      nodes;
    std::vector<NetworkEdge>      edges;
    std::vector<int>              degree;
    std::vector<double>           betweenness;      // Brandes
    std::vector<int>              community;        // Louvain assignment
    int                           componentCount = 0;
    double                        modularity = 0;
    std::vector<std::string>      warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(GraphMetrics, nodes, edges, degree, betweenness, community,
                                   componentCount, modularity, warnings)

}  // namespace biocad
