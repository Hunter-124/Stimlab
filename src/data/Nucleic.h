// data/Nucleic.h - JSON-serializable DTOs for the DNA/RNA workbench: sequences
// and features, restriction maps, translation, oligo thermodynamics and codon
// optimization.
//
// BIOSECURITY BOUNDARY, permanent: there is no synthesis-vendor integration, no
// order-sheet export, no pathogen-driven batch design, no therapeutic or germline
// CRISPR framing, and no "synthesise this" affordance anywhere in these types.
// Export is FASTA and GenBank only. Codon optimization is labelled constraint
// satisfaction and never an expression prediction - the DTO carries no predicted
// yield or titre field to put one in.
//
// PROVENANCE RULE: a guide-search result carries the exact reference it searched
// and the number of bases actually examined, because "no off-target found" in
// 2686 bp of plasmid is not a genome-wide specificity claim.
#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "data/Domain.h"

namespace biocad {

enum class NucKind { Dna, Rna };

NLOHMANN_JSON_SERIALIZE_ENUM(NucKind, {
    {NucKind::Dna, "dna"},
    {NucKind::Rna, "rna"},
})

enum class Strand { Forward, Reverse };

NLOHMANN_JSON_SERIALIZE_ENUM(Strand, {
    {Strand::Forward, "+"},
    {Strand::Reverse, "-"},
})

// One GenBank-style feature interval. `parts` holds a join()/order() location as
// half-open [begin, end) 0-based pairs so a spliced CDS is one feature rather than
// several; a single-span feature has exactly one part.
struct NucFeature {
    std::string                                  type;      // CDS, gene, promoter...
    std::vector<std::pair<int, int>>             parts;     // [begin, end), 0-based
    Strand                                       strand = Strand::Forward;
    std::vector<std::pair<std::string, std::string>> qualifiers;  // /gene=, /translation=
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(NucFeature, type, parts, strand, qualifiers)

struct NucRecord {
    std::string              id;
    std::string              description;
    NucKind                  kind = NucKind::Dna;
    std::string              sequence;     // IUPAC, uppercase
    bool                     circular = false;
    std::vector<NucFeature>  features;
    std::vector<std::string> warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(NucRecord, id, description, kind, sequence, circular,
                                   features, warnings)

// ---------------------------------------------------------------------------
// Restriction mapping.
// ---------------------------------------------------------------------------

struct RestrictionSite {
    std::string enzyme;
    std::string recognition;   // IUPAC pattern as written in the pack
    int         position = 0;  // 0-based index of the cut on the top strand
    Strand      strand = Strand::Forward;
    int         overhang = 0;  // >0 5' overhang, <0 3' overhang, 0 blunt
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictionSite, enzyme, recognition, position, strand,
                                   overhang)

// A digest. `fragmentLengths` sums to the sequence length exactly - for a circular
// template as well, which is the arithmetic that catches an off-by-one in the
// wrap-around cut.
struct RestrictionDigest {
    std::string                  recordId;
    std::vector<RestrictionSite> sites;
    std::vector<int>             fragmentLengths;  // descending, as a gel would show
    bool                         circular = false;
    std::vector<std::string>     warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RestrictionDigest, recordId, sites, fragmentLengths, circular,
                                   warnings)

// ---------------------------------------------------------------------------
// Translation.
// ---------------------------------------------------------------------------

struct OpenReadingFrame {
    int         begin = 0;      // 0-based, inclusive of the start codon
    int         end = 0;        // exclusive, inclusive of the stop codon
    Strand      strand = Strand::Forward;
    int         frame = 0;      // 0, 1, 2 relative to the strand's 5' end
    std::string protein;        // one-letter, no trailing stop
    bool        stopped = false;  // false when the ORF ran off the end
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(OpenReadingFrame, begin, end, strand, frame, protein, stopped)

struct TranslationResult {
    std::string                   recordId;
    int                           geneticCodeId = 1;   // NCBI transl_table
    std::vector<OpenReadingFrame> orfs;                // descending length
    std::vector<std::string>      frames;              // six-frame translation
    std::vector<std::string>      warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(TranslationResult, recordId, geneticCodeId, orfs, frames,
                                   warnings)

// ---------------------------------------------------------------------------
// Oligo thermodynamics.
// ---------------------------------------------------------------------------

// Nearest-neighbour thermodynamics for one oligo. dH/dS/dG carry their units and
// the SantaLucia parameter set as the Quantity source: a Tm quoted without its
// salt, oligo concentration and parameter set is not reproducible.
struct OligoThermo {
    std::string sequence;
    Quantity    deltaH;       // kcal/mol
    Quantity    deltaS;       // cal/(mol*K)
    Quantity    deltaG37;     // kcal/mol
    Quantity    tm;           // degrees C
    double      gcPercent = 0;
    double      naMolar = 0;      // monovalent salt actually used
    double      mgMolar = 0;
    double      oligoMolar = 0;
    double      dntpMolar = 0;
    std::vector<std::string> assumptions;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(OligoThermo, sequence, deltaH, deltaS, deltaG37, tm,
                                   gcPercent, naMolar, mgMolar, oligoMolar, dntpMolar,
                                   assumptions)

struct SecondaryStructure {
    std::string kind;         // "hairpin", "self-dimer", "hetero-dimer"
    Quantity    deltaG37;     // kcal/mol; more negative is more stable
    int         position = 0;
    std::string alignment;    // two-line ASCII rendering of the duplex
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SecondaryStructure, kind, deltaG37, position, alignment)

struct PrimerPair {
    OligoThermo                     forwardOligo;
    OligoThermo                     reverseOligo;
    int                             productBegin = 0;
    int                             productEnd = 0;
    int                             productLength = 0;
    double                          tmDifference = 0;
    std::vector<SecondaryStructure> liabilities;
    std::vector<std::string>        warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PrimerPair, forwardOligo, reverseOligo, productBegin,
                                   productEnd, productLength, tmDifference, liabilities, warnings)

// ---------------------------------------------------------------------------
// Codon usage and constraint-based optimization.
// ---------------------------------------------------------------------------

struct CodonUsageEntry {
    std::string codon;
    std::string aminoAcid;
    double      relativeAdaptiveness = 0;   // w_ij, 0..1 within the synonymous family
    double      frequencyPerThousand = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(CodonUsageEntry, codon, aminoAcid, relativeAdaptiveness,
                                   frequencyPerThousand)

struct CodonMetrics {
    Quantity    cai;              // geometric mean of w_ij; Measured arithmetic on a cited table
    Quantity    gcPercent;
    Quantity    gc3Percent;
    std::string usageTableName;   // the cited table the metrics are relative to
    std::vector<std::string> warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(CodonMetrics, cai, gcPercent, gc3Percent, usageTableName,
                                   warnings)

// The result of constraint satisfaction, NOT a prediction. There is deliberately
// no expression, yield or titre field: the guarantee this type makes is that
// `optimized` translates to exactly the input protein and contains none of the
// forbidden sites, and that is the only guarantee it makes.
struct CodonOptimizationResult {
    std::string              optimized;
    std::string              protein;          // back-translation check, must equal the input
    bool                     translationPreserved = false;
    std::vector<std::string> forbiddenSites;   // patterns that were excluded
    std::vector<std::string> remainingViolations;  // empty on success
    CodonMetrics             before;
    CodonMetrics             after;
    std::vector<std::string> assumptions;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(CodonOptimizationResult, optimized, protein,
                                   translationPreserved, forbiddenSites, remainingViolations,
                                   before, after, assumptions)

// ---------------------------------------------------------------------------
// Guide search - with its reference scope attached, always.
// ---------------------------------------------------------------------------

struct GuideCandidate {
    std::string protospacer;     // 20 nt
    std::string pam;
    int         position = 0;
    Strand      strand = Strand::Forward;
    double      gcPercent = 0;
    int         exactOffTargets = 0;      // within the reference ACTUALLY searched
    int         oneMismatchOffTargets = 0;
    int         twoMismatchOffTargets = 0;
    std::vector<std::string> warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(GuideCandidate, protospacer, pam, position, strand, gcPercent,
                                   exactOffTargets, oneMismatchOffTargets, twoMismatchOffTargets,
                                   warnings)

// `referenceName` and `basesSearched` are required fields, not decoration: an
// off-target count is meaningless without the scope it was counted in, and
// `genomeWideClaimPossible` is false unless a complete genome was supplied.
struct GuideSearchResult {
    std::string                 referenceName;
    std::int64_t                basesSearched = 0;
    bool                        genomeWideClaimPossible = false;
    std::string                 scopeStatement;
    std::vector<GuideCandidate> guides;
    std::vector<std::string>    warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(GuideSearchResult, referenceName, basesSearched,
                                   genomeWideClaimPossible, scopeStatement, guides, warnings)

}  // namespace biocad
