// data/Variants.h - JSON-serializable DTOs for protein variant analysis and
// point-mutation modelling.
//
// HONESTY SCOPE, permanent. This is the least certain area BioCAD touches, so the
// types are shaped to make the uncertainty unavoidable:
//  - A conservation score carries the homolog COUNT and the median identity of the
//    alignment it came from, and is NotComputed below a minimum homolog count. A
//    conservation score from five ad-hoc sequences is an arbitrary number wearing
//    a threshold.
//  - A rebuilt side chain is Provenance::Model: a constructed artefact with NO
//    energy claim attached. There is no field here for a rebuild score.
//  - Any predicted stability change carries its benchmark error IN the value and
//    the model's positive predictive value beside it, because a coin flip
//    presented as a number is the failure mode of this whole area.
//  - There is no pathogenicity verdict, no clinical significance field, and no
//    binder-design output. Those absences are deliberate.
#pragma once

#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "data/Domain.h"

namespace biocad {

// The alignment a conservation score was computed FROM. Without these numbers the
// score is not interpretable, which is why they are part of every result rather
// than a separate query.
struct HomologSetSummary {
    int         sequenceCount = 0;
    double      medianIdentityPct = 0;
    double      minIdentityPct = 0;
    double      maxIdentityPct = 0;
    int         alignmentColumns = 0;
    int         effectiveSequenceCount = 0;   // after identity-weighting
    std::string source;                        // where the homologs came from
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(HomologSetSummary, sequenceCount, medianIdentityPct,
                                   minIdentityPct, maxIdentityPct, alignmentColumns,
                                   effectiveSequenceCount, source)

// Per-column conservation. `shannonEntropy` is in bits over the 20 amino acids;
// `pssm` is the log-odds column against the background frequencies actually used.
struct ConservationColumn {
    int                 position = 0;      // 1-based in the query numbering
    char                queryResidue = 'X';
    double              shannonEntropy = 0;
    double              gapFraction = 0;
    std::vector<double> frequencies;       // 20 entries, alphabetical one-letter order
    std::vector<double> pssm;              // 20 entries, log-odds
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ConservationColumn, position, queryResidue, shannonEntropy,
                                   gapFraction, frequencies, pssm)

struct ConservationProfile {
    std::string                     queryId;
    HomologSetSummary               homologs;
    std::vector<ConservationColumn> columns;
    // Below the minimum homolog count nothing here is usable, and `usable` is
    // false rather than the caller being trusted to check the count.
    bool                            usable = false;
    int                             minimumHomologsRequired = 0;
    std::string                     backgroundFrequencySource;
    // The query sequence itself, kept even when the profile is unusable. Without
    // it a BLOSUM62 delta - which does not depend on the homolog set at all -
    // would be unavailable purely because the conservation half was refused, and
    // the cheap honest number would be lost with the expensive uncertain one.
    std::string                     query;
    std::vector<std::string>        warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ConservationProfile, queryId, homologs, columns, usable,
                                   minimumHomologsRequired, backgroundFrequencySource, query,
                                   warnings)

// One substitution's conservation-based scores. The two published thresholds are
// carried WITH the scores so the reader is not left to remember them: SIFT is
// deleterious below 0.05 and PROVEAN below -2.282 (the latter at a balanced
// accuracy of 79.05% on its 58,684-variant human validation set).
struct VariantScore {
    int         position = 0;
    char        wildType = 'X';
    char        mutant = 'X';
    Quantity    blosum62Delta;      // Measured: an exact table lookup difference
    Quantity    siftScore;          // p(mutant)/max_y p(y)
    Quantity    proveanScore;
    Quantity    columnEntropy;
    double      siftDeleteriousBelow = 0.05;
    double      proveanDeleteriousBelow = -2.282;
    HomologSetSummary homologs;
    // Deliberately NOT a pathogenicity call. This string states what the numbers
    // are and what they are not.
    std::string interpretation;
    std::vector<std::string> warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(VariantScore, position, wildType, mutant, blosum62Delta,
                                   siftScore, proveanScore, columnEntropy, siftDeleteriousBelow,
                                   proveanDeleteriousBelow, homologs, interpretation, warnings)

// A rebuilt side chain. `rotamerLibrarySource` is required: the Dunbrack
// backbone-dependent library is CC BY 4.0 and redistributable with attribution,
// and a rebuild with no named library is not reproducible.
//
// There is no energy field. The rebuild is a geometric construction chosen by
// dead-end elimination over the library, and attaching a score to it would imply
// a force field that BioCAD does not have.
struct RotamerRebuild {
    int                      position = 0;
    char                     wildType = 'X';
    char                     mutant = 'X';
    std::vector<double>      chiAngles;          // degrees, chi1..chi4 as applicable
    double                   rotamerProbability = 0;   // from the library, not an energy
    int                      clashCount = 0;     // heavy-atom overlaps after repack
    std::vector<std::string> repackedNeighbours;
    Provenance               provenance = Provenance::Model;
    std::string              rotamerLibrarySource;
    std::vector<std::string> assumptions;
    std::vector<std::string> warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RotamerRebuild, position, wildType, mutant, chiAngles,
                                   rotamerProbability, clashCount, repackedNeighbours,
                                   provenance, rotamerLibrarySource, assumptions, warnings)

// A predicted stability change. `positivePredictiveValuePct` is a REQUIRED field
// because it is the number that stops a stabilising prediction from being read as
// a result: for the current best CPU-feasible models it is 46-56%, which is
// roughly a coin flip, and the panel says so.
struct StabilityPrediction {
    int         position = 0;
    char        wildType = 'X';
    char        mutant = 'X';
    Quantity    deltaDeltaG;                  // kcal/mol, error = the benchmark RMSE
    std::string modelName;
    std::string benchmarkName;
    double      benchmarkRmse = 0;
    double      benchmarkPearson = 0;
    double      positivePredictiveValuePct = 0;
    std::vector<std::string> assumptions;
    std::vector<std::string> warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(StabilityPrediction, position, wildType, mutant, deltaDeltaG,
                                   modelName, benchmarkName, benchmarkRmse, benchmarkPearson,
                                   positivePredictiveValuePct, assumptions, warnings)

}  // namespace biocad
