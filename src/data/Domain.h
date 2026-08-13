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
    double      logS = 0;                 // aqueous solubility (logS, mol/L)
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

}  // namespace biocad
