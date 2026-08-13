// data/Domain.h - JSON-serializable domain types shared across the whole app.
// These are the "lingua franca": modules produce them, the UI renders them,
// storage persists them. Pure data; no behavior beyond (de)serialization.
//
// SAFETY SCOPE: every type here describes WHAT A COMPOUND IS AND DOES
// (identity, properties, pharmacology, stability, absorption/PK, similarity,
// legal status). There is deliberately NO synthesis/route/manufacturability
// type anywhere in the domain - that is out of scope by design.
#pragma once

#include <string>
#include <vector>

#include <nlohmann/json.hpp>

namespace biocad {

// ---------------------------------------------------------------------------
// Shared verdict scale used by ADMET / absorption / safety flags.
// ---------------------------------------------------------------------------
enum class Verdict { Info, Good, Warn, Danger };

NLOHMANN_JSON_SERIALIZE_ENUM(Verdict, {
    {Verdict::Info, "info"},
    {Verdict::Good, "good"},
    {Verdict::Warn, "warn"},
    {Verdict::Danger, "danger"},
})

const char* verdictLabel(Verdict v);

// ---------------------------------------------------------------------------
// Provenance - why a number may be trusted. This is the single most important
// rule in the domain: no derived number is rendered without one.
// The enumerators are ordered strongest-first, so the *weakest* of a set of
// inputs (which is what a derived quantity inherits) is std::max, not std::min;
// use weakest() rather than open-coding the comparison.
// ---------------------------------------------------------------------------
enum class Provenance {
    Measured,     // exact geometry/statistics, or an experimental value retrieved with a citation
    Predicted,    // a published model actually ran; physical units; benchmark error is mandatory
    Model,        // a constructed artefact (built structure, docked pose) with no energy claim
    Heuristic,    // rank-ordering only; arbitrary units; physical units are forbidden
    NotComputed   // a prerequisite was missing
};

NLOHMANN_JSON_SERIALIZE_ENUM(Provenance, {
    {Provenance::Measured,    "measured"},
    {Provenance::Predicted,   "predicted"},
    {Provenance::Model,       "model"},
    {Provenance::Heuristic,   "heuristic"},
    {Provenance::NotComputed, "not computed"},
})

// "measured" | "predicted" | "model" | "heuristic" | "not computed"
const char* provenanceLabel(Provenance p);

// The weaker (less trustworthy) of two tiers. A derived quantity is never more
// trustworthy than its worst input.
constexpr Provenance weakest(Provenance a, Provenance b) {
    return static_cast<int>(a) >= static_cast<int>(b) ? a : b;
}

// A number that knows how well it is known. `unit` MUST be empty when
// provenance == Heuristic: a rank-ordering score has no physical units, and
// makeQuantity() enforces that rather than merely documenting it.
struct Quantity {
    double      value = 0.0;
    std::string unit;          // empty string is REQUIRED when provenance == Heuristic
    double      error = 0.0;   // 0 means "no error bar available"
    Provenance  provenance = Provenance::NotComputed;
    std::string source;        // citation, engine name, or model+benchmark
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Quantity, value, unit, error, provenance, source)

// The only sanctioned way to build a Quantity. Throws std::invalid_argument when
// a Heuristic carries a unit - that combination is a programming error, and making
// it unrepresentable is the point of this type.
Quantity makeQuantity(double value, std::string unit, double error, Provenance p,
                      std::string source);

// A quantity whose prerequisite was missing. `source` names what is missing.
Quantity notComputed(std::string missingPrerequisite);

// ---------------------------------------------------------------------------
// Molecule - the central identity + computed-property record.
// ---------------------------------------------------------------------------
struct Molecule {
    std::string id;            // stable slug, e.g. "amphetamine"
    std::string name;          // display name
    std::string smiles;        // canonical SMILES
    std::string formula;       // molecular formula, e.g. "C9H13N"
    double      molWeight = 0; // g/mol
    double      logP = 0;      // octanol-water partition (lipophilicity)
    double      tpsa = 0;      // topological polar surface area (A^2)
    int         hbd = 0;       // H-bond donors
    int         hba = 0;       // H-bond acceptors
    int         rotatableBonds = 0;
    std::string drugClass;     // e.g. "Phenethylamine stimulant"
    std::string legalStatus;   // e.g. "Schedule II (US)"
    std::string notes;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Molecule, id, name, smiles, formula, molWeight,
                                   logP, tpsa, hbd, hba, rotatableBonds, drugClass,
                                   legalStatus, notes)

// ---------------------------------------------------------------------------
// Stability - replaces the out-of-scope "manufacturability" score.
// Higher overallScore = more chemically stable (slower degradation).
// ---------------------------------------------------------------------------
struct StabilityFactor {
    std::string name;        // e.g. "Hydrolysis", "Oxidation", "Photolysis"
    double      score = 0;   // 0..100, higher = more resistant
    std::string rationale;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(StabilityFactor, name, score, rationale)

struct Degradant {
    std::string name;        // likely degradation product
    std::string pathway;     // e.g. "N-oxidation", "Ester hydrolysis"
    std::string note;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Degradant, name, pathway, note)

struct StabilityReport {
    std::string moleculeId;
    double      overallScore = 0;          // 0..100
    std::string shelfLifeEstimate;          // human-readable, e.g. "~24 months @ 25C/60%RH"
    std::vector<StabilityFactor> factors;
    std::vector<Degradant>       degradants;
    std::string summary;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(StabilityReport, moleculeId, overallScore,
                                   shelfLifeEstimate, factors, degradants, summary)

// ---------------------------------------------------------------------------
// ADMET / metabolism - "what the body does to the drug" (D, M, E, T).
// (Absorption gets its own richer report below.)
// ---------------------------------------------------------------------------
struct AdmetEndpoint {
    std::string name;        // e.g. "CYP2D6 substrate", "hERG liability"
    Verdict     verdict = Verdict::Info;
    std::string detail;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(AdmetEndpoint, name, verdict, detail)

struct AdmetReport {
    std::string moleculeId;
    Verdict     overall = Verdict::Info;
    std::vector<AdmetEndpoint> endpoints;
    std::string summary;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(AdmetReport, moleculeId, overall, endpoints, summary)

// ---------------------------------------------------------------------------
// Structural alerts - LIABILITY FLAGS, never toxicity verdicts.
//
// A flag means: this substructure has been associated with reactive-metabolite
// formation in the medicinal-chemistry literature. It does NOT mean the compound
// is toxic, and it is not a prediction that bioactivation occurs - that depends
// on the enzymes present, the competing clearance routes, the dose and the
// detoxication capacity, none of which a substructure knows. Many marketed drugs
// match several alerts.
//
// `severity` is therefore Verdict::Info or Verdict::Warn ONLY. Verdict::Danger is
// not reachable from this type by construction: a substructure match cannot
// support a danger verdict, so the producer never emits one and no consumer needs
// to render one.
// ---------------------------------------------------------------------------
struct AlertFlag {
    std::string label;      // e.g. "para-Aminophenol / anilide (quinone-imine former)"
    std::string mechanism;  // the metabolic route, stated as a route
    std::string citation;   // never empty: an uncitable alert is not an alert
    int         atomCount = 0;  // atoms of this molecule that raised the flag
    Verdict     severity = Verdict::Info;  // Info or Warn only
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(AlertFlag, label, mechanism, citation, atomCount, severity)

struct AlertReport {
    std::string            moleculeId;
    std::vector<AlertFlag> flags;
    std::string            summary;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(AlertReport, moleculeId, flags, summary)

// ---------------------------------------------------------------------------
// Known metabolites - CURATED FACTS, not enumerated hypotheses.
//
// Every entry here is a transformation that has actually been characterised in a
// citable source, which is why it can be rendered Provenance::Measured. The
// alternative - enumerating a rule-based metabolite tree - was measured by an
// independent EPA cross-tool benchmark at 1.1-29% precision and 14.7-28.3%
// sensitivity (Boyce et al. 2022, Comput Toxicol 21:100208), so the two must never
// share a surface: a hypothesis rendered next to a fact inherits its credibility.
// ---------------------------------------------------------------------------
struct MetaboliteFact {
    std::string parentId;          // library compound id this metabolite comes from
    std::string metaboliteName;    // display name of the product
    std::string metaboliteSmiles;  // empty when the structure is not certain - never guessed
    std::string enzyme;            // e.g. "CYP2D6", "UGT2B7", "Butyrylcholinesterase"
    std::string reaction;          // e.g. "O-demethylation at the 3-position"
    std::string significance;      // why this transformation matters
    std::string citation;          // real, checkable reference; never fabricated
    bool        polymorphic = false;  // the enzyme has clinically relevant genetic variation
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MetaboliteFact, parentId, metaboliteName,
                                   metaboliteSmiles, enzyme, reaction, significance,
                                   citation, polymorphic)

struct MetabolismReport {
    std::string moleculeId;
    std::vector<MetaboliteFact> known;
    std::string summary;
    // ALWAYS populated, and rendered whether or not `known` is empty. An empty
    // list means BioCAD has no curated fact for this compound - it does NOT mean
    // the compound has no metabolites. Implying the latter is the single easiest
    // lie this surface could tell, so the note is a required field rather than a
    // conditional string the renderer might skip.
    std::string coverageNote;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MetabolismReport, moleculeId, known, summary,
                                   coverageNote)

// ---------------------------------------------------------------------------
// Absorption / Pharmacokinetics - the "A" of ADMET, modeled in depth so it can
// discriminate between candidate analogs (per user request).
// ---------------------------------------------------------------------------
struct PkMetric {
    std::string name;        // e.g. "Caco-2 Papp", "Oral bioavailability"
    double      value = 0;
    std::string unit;        // e.g. "%", "log(cm/s)", "logBB"
    Verdict     band = Verdict::Info;  // qualitative band for coloring
    std::string rationale;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PkMetric, name, value, unit, band, rationale)

struct AbsorptionReport {
    std::string moleculeId;
    double      bioavailabilityPct = 0;   // predicted oral F (%)
    double      hiaPct = 0;               // human intestinal absorption (%)
    double      caco2LogPapp = 0;         // log(cm/s) apparent permeability
    double      logBB = 0;                // blood-brain barrier partition (logBB)
    // Aqueous solubility as log10 mol/L. This was a fabricated regression on MW and
    // logP; it is now a Quantity so an absent melting point reads notComputed rather
    // than a plausible-looking number. The JSON key is unchanged on purpose.
    Quantity    logS;
    bool        pgpSubstrate = false;     // P-glycoprotein efflux substrate
    bool        cnsPenetrant = false;     // crosses BBB meaningfully
    std::vector<PkMetric> metrics;        // for the PK dashboard / charts
    std::string summary;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(AbsorptionReport, moleculeId, bioavailabilityPct,
                                   hiaPct, caco2LogPapp, logBB, logS, pgpSubstrate,
                                   cnsPenetrant, metrics, summary)

// ---------------------------------------------------------------------------
// Similarity vs a curated known-substance reference library.
// ---------------------------------------------------------------------------
struct SimilarityHit {
    std::string referenceName;
    std::string referenceClass;
    std::string legalStatus;
    double      tanimoto = 0;          // 0..1 fingerprint similarity
    double      pharmacophore = 0;     // 0..1 pharmacophore overlap
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SimilarityHit, referenceName, referenceClass,
                                   legalStatus, tanimoto, pharmacophore)

struct SimilarityReport {
    std::string moleculeId;
    std::string nearestName;
    double      nearestScore = 0;
    std::vector<SimilarityHit> hits;     // ranked descending
    std::string summary;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SimilarityReport, moleculeId, nearestName,
                                   nearestScore, hits, summary)

// ---------------------------------------------------------------------------
// Legal-analog "substantially similar" scorecard (policy/forensic use).
// ---------------------------------------------------------------------------
struct LegalScorecard {
    std::string moleculeId;
    std::string jurisdiction;            // e.g. "US Federal Analogue Act"
    double      substantialSimilarity = 0;  // 0..100
    std::string classification;          // e.g. "Likely controlled analog"
    std::vector<std::string> rationale;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(LegalScorecard, moleculeId, jurisdiction,
                                   substantialSimilarity, classification, rationale)

// ---------------------------------------------------------------------------
// Docking = ligand->protein BINDING AFFINITY (pharmacology/activity).
// NOTE: binding energy is a target-engagement signal, never a make-it signal.
// ---------------------------------------------------------------------------
struct Pose {
    int    rank = 0;
    double affinityKcalPerMol = 0;   // more negative = stronger predicted binding
    double rmsd = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Pose, rank, affinityKcalPerMol, rmsd)

struct DockingResult {
    std::string moleculeId;
    std::string targetName;          // e.g. "DAT (dopamine transporter)"
    double      bestAffinity = 0;
    std::vector<Pose> poses;
    std::string summary;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(DockingResult, moleculeId, targetName,
                                   bestAffinity, poses, summary)

// ---------------------------------------------------------------------------
// Run history (for the Runs/Results page).
// ---------------------------------------------------------------------------
struct RunRecord {
    std::string id;
    std::string kind;        // "Docking", "ADMET", "Stability", "Absorption", ...
    std::string subject;     // molecule / target description
    std::string status;      // "complete", "running", "failed"
    std::string createdAt;   // ISO-ish string
    std::string summary;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RunRecord, id, kind, subject, status, createdAt, summary)

// ---------------------------------------------------------------------------
// Pharmacodynamics and pharmacokinetics.
//
// SAFETY SCOPE: everything here is an EXPOSURE SCENARIO under stated assumptions.
// None of it is a dose recommendation, and the UI and the agent are both forbidden
// from turning it into one.
// ---------------------------------------------------------------------------

// One concentration/effect observation. `concentration` is molar; a non-positive
// value cannot be log-transformed and is rejected by the fitter.
struct DoseResponsePoint {
    double concentration = 0.0;  // mol/L
    double effect = 0.0;         // assay response, arbitrary units
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(DoseResponsePoint, concentration, effect)

// A four-parameter logistic fit. `hillSlope` is deliberately named an EMPIRICAL
// SLOPE, never "cooperativity": in a functional assay amplification and receptor
// reserve bend the slope, so a value above 1 is not evidence of cooperative binding.
struct CurveFit {
    Quantity top;         // upper asymptote, assay units
    Quantity bottom;      // lower asymptote, assay units
    Quantity ec50;        // mol/L
    Quantity hillSlope;   // empirical slope, dimensionless
    double   rSquared = 0.0;
    int      iterations = 0;
    bool     converged = false;
    std::string note;     // why a fit failed, or what was assumed
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(CurveFit, top, bottom, ec50, hillSlope, rSquared,
                                   iterations, converged, note)

// How an inhibitor acts. There is no default: competitive and uncompetitive differ
// by 10x at [S] = 10*Km and by 100x at [S] = 100*Km, so guessing would silently
// magnitude. ChEMBL has no UNCOMPETITIVE action type to disambiguate with.
// Unknown and Mixed exist for the one producer that can legitimately return them:
// the Phase 10 global fit over the full [S] x [I] matrix, which selects a modality
// by AICc and answers Unknown when the difference is under 2. Cheng-Prusoff
// refuses to convert an IC50 under Unknown rather than falling back to
// competitive.
enum class InhibitionModality { Competitive, Noncompetitive, Uncompetitive, Mixed,
                               RadioligandBinding, Unknown };

NLOHMANN_JSON_SERIALIZE_ENUM(InhibitionModality, {
    {InhibitionModality::Competitive,        "competitive"},
    {InhibitionModality::Noncompetitive,     "noncompetitive"},
    {InhibitionModality::Uncompetitive,      "uncompetitive"},
    {InhibitionModality::Mixed,              "mixed"},
    {InhibitionModality::RadioligandBinding, "radioligand-binding"},
    {InhibitionModality::Unknown,            "unknown"},
})

// Cheng-Prusoff input. The modality decides which fields are REQUIRED; a missing
// one yields a NotComputed Quantity naming it rather than an assumed competitive fit.
struct ChengPrusoffInput {
    InhibitionModality modality = InhibitionModality::Competitive;
    double ic50 = 0.0;              // mol/L
    double substrate = -1.0;        // [S], mol/L; < 0 = absent
    double km = -1.0;               // Km, mol/L; < 0 = absent
    double radioligand = -1.0;      // [L*], mol/L; < 0 = absent
    double kdRadioligand = -1.0;    // Kd of the radioligand, mol/L; < 0 = absent
    double enzymeConc = -1.0;       // [E]t, mol/L; >= 0 enables the tight-binding view
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ChengPrusoffInput, modality, ic50, substrate, km,
                                   radioligand, kdRadioligand, enzymeConc)

// One Schild point: an antagonist concentration and the agonist dose ratio it caused.
struct SchildPoint {
    double antagonist = 0.0;  // [B], mol/L
    double doseRatio = 1.0;   // DR, dimensionless
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SchildPoint, antagonist, doseRatio)

// A Schild regression. KB is only meaningful when the slope is indistinguishable
// from 1; otherwise the antagonism is not simple competitive and `kbUsable` is false.
struct SchildResult {
    Quantity pA2;          // -log10 of the [B] giving DR = 2
    Quantity slope;        // Schild slope, dimensionless
    double   slopeCiLow = 0.0;
    double   slopeCiHigh = 0.0;
    Quantity kb;           // NotComputed when the slope CI excludes 1
    bool     kbUsable = false;
    std::string note;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SchildResult, pA2, slope, slopeCiLow, slopeCiHigh,
                                   kb, kbUsable, note)

// Which PK structural model to integrate.
enum class PkModel { IvBolus, IvInfusion, OralOneCompartment, OralTwoCompartment };

NLOHMANN_JSON_SERIALIZE_ENUM(PkModel, {
    {PkModel::IvBolus,            "iv-bolus"},
    {PkModel::IvInfusion,         "iv-infusion"},
    {PkModel::OralOneCompartment, "oral-1c"},
    {PkModel::OralTwoCompartment, "oral-2c"},
})

// One dosing event in the simulated regimen.
struct DoseEvent {
    double timeH = 0.0;        // hours from t0
    double amountMg = 0.0;
    double durationH = 0.0;    // 0 = bolus / instantaneous oral input
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(DoseEvent, timeH, amountMg, durationH)

// A PK parameter set. Every field carries its own Provenance because F, ka and fu
// have NO credible structure-only predictor: they default to assumed, and the panel
// says so rather than implying they were computed.
struct PkModelSpec {
    PkModel  model = PkModel::OralOneCompartment;
    Quantity bioavailability;      // F, dimensionless 0..1
    Quantity absorptionRate;       // ka, 1/h
    Quantity clearance;            // CL, L/h
    Quantity volume;               // V (central), L
    Quantity volumePeripheral;     // V2, L (two-compartment only)
    Quantity intercompartmental;   // Q, L/h (two-compartment only)
    Quantity unboundFraction;      // fu, dimensionless
    // Optional Michaelis-Menten elimination; used when vmax > 0.
    Quantity vmax;                 // mg/h
    Quantity km;                   // mg/L
    double   stepH = 0.01;         // RK4 fixed step, hours
    double   horizonH = 24.0;      // simulated span, hours
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PkModelSpec, model, bioavailability, absorptionRate,
                                   clearance, volume, volumePeripheral, intercompartmental,
                                   unboundFraction, vmax, km, stepH, horizonH)

struct DoseRegimen {
    std::vector<DoseEvent> doses;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(DoseRegimen, doses)

// A simulated exposure profile. `assumptions` is rendered verbatim under the plot:
// a curve whose assumptions are not visible is a curve that misleads.
struct PkProfile {
    std::vector<double> timeH;
    std::vector<double> concentrationMgPerL;   // total plasma concentration
    std::vector<double> unboundMgPerL;         // fu * total
    Quantity cmax;         // mg/L
    Quantity tmax;         // h
    Quantity auc;          // mg*h/L, trapezoidal over the simulated horizon
    Quantity halfLife;     // h
    Quantity accumulation; // Rac, dimensionless (multiple dosing only)
    bool flipFlop = false; // ka < ke: the terminal phase reflects absorption, not elimination
    std::vector<std::string> assumptions;
    std::string note;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PkProfile, timeH, concentrationMgPerL, unboundMgPerL,
                                   cmax, tmax, auc, halfLife, accumulation, flipFlop,
                                   assumptions, note)

// Fractional target occupancy over time: theta = [A] / (Kd + [A]).
// This is the honest synthesis of PK and PD - it needs no Emax, no tau and no
// tissue assumption, only a free concentration and a Kd.
struct OccupancyCurve {
    std::vector<double> timeH;
    std::vector<double> occupancy;   // 0..1
    Quantity peakOccupancy;
    Quantity timeAbove50Pct;         // hours above 50% occupancy
    std::string note;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(OccupancyCurve, timeH, occupancy, peakOccupancy,
                                   timeAbove50Pct, note)

// ---------------------------------------------------------------------------
// Protein core (Phase 5) results. These are the UI/agent-facing shapes; the
// geometry and alignment maths live in bio:: and never leak their own types
// through a contract.
// ---------------------------------------------------------------------------

// One pairwise sequence alignment, already rendered to strings so the panel and
// the agent see the same characters. `midline` is the BLAST-style match line
// ('|' identity, '+' positive substitution, ' ' otherwise) and is exactly as
// long as the two aligned strings.
//
// eValue is notComputed("local alignment required") for a global alignment:
// Karlin-Altschul statistics are defined for local alignments only, so a global
// E-value would be a number with no meaning. The UI omits the row entirely.
struct SequenceAlignment {
    std::string aligned1;
    std::string aligned2;
    std::string midline;
    Quantity    score;            // alignment score in the matrix' half-bit units
    Quantity    identityPct;
    Quantity    similarityPct;    // positive-scoring pairs; a DIFFERENT number to identity
    Quantity    eValue;
    int         gapOpens = 0;
    int         alignedLength = 0;
    std::string note;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SequenceAlignment, aligned1, aligned2, midline, score,
                                   identityPct, similarityPct, eValue, gapOpens,
                                   alignedLength, note)

// One structure-vs-structure comparison. rmsd and lddt are exact geometry
// (Provenance::Measured). tmScore stays notComputed until Phase 5.8 vendors
// TM-align: an approximated TM-score is worse than none, because the number is
// read as if it came from the reference implementation.
struct StructureComparison {
    Quantity    rmsd;             // Angstrom, over the matched atom set
    Quantity    lddt;             // 0..1, superposition-free
    Quantity    tmScore;
    Quantity    alignedResidues;  // count of residues actually matched
    int         unmatchedResidues = 0;
    std::string note;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(StructureComparison, rmsd, lddt, tmScore, alignedResidues,
                                   unmatchedResidues, note)

}  // namespace biocad
