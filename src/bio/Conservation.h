// bio/Conservation.h - homolog-set profiling, per-column entropy, a PSSM, and the
// two conservation-based substitution scores (SIFT-style and PROVEAN-style).
//
// WHY THERE IS A MINIMUM HOMOLOG COUNT. A conservation score is an estimate of a
// column's amino-acid distribution. With five sequences the estimate has at most
// five distinct observations, the pseudocount dominates, and the resulting SIFT
// ratio is essentially a function of the pseudocount scheme rather than of the
// protein - yet it is read against a published 0.05 threshold as if it were the
// same quantity SIFT computes over a PSI-BLAST hit list of dozens to hundreds of
// sequences. The two published methods are explicit about this: SIFT's own
// pipeline warns when the median sequence information of the alignment is too
// high (too few, too similar sequences), and PROVEAN's threshold of -2.282 was
// calibrated on clusters of up to 45 supporting sequences per query. So this file
// sets kMinimumHomologs = 15 and refuses below it: `ConservationProfile::usable`
// is false, and every VariantScore quantity derived from such a profile is
// notComputed() naming the shortfall. 15 is a judgement, not a published
// constant, and it is stated as one in the profile's warnings.
//
// The alignment is the tree's ONE aligner (bio::alignGlobal, Gotoh affine gaps,
// BLOSUM62); nothing here re-implements alignment.
#pragma once

#include <array>
#include <string>
#include <vector>

#include "bio/Align.h"
#include "data/Variants.h"

namespace biocad::bio {

// The 20 standard amino acids in alphabetical one-letter order. This order is the
// index order of ConservationColumn::frequencies and ::pssm, so it is defined once
// here and never re-spelled.
inline constexpr const char* kAminoAcids = "ACDEFGHIKLMNPQRSTVWY";
inline constexpr int kAminoAcidCount = 20;

// Below this many homologs a profile is not usable. See the file header for why.
inline constexpr int kMinimumHomologs = 15;

// -1 for anything that is not one of the 20 standard one-letter codes.
[[nodiscard]] int aminoAcidIndex(char c);

// Background amino-acid frequencies. The set actually used travels with the
// profile in `backgroundFrequencySource`, because a log-odds column against an
// unnamed background is not reproducible.
struct BackgroundFrequencies {
    std::array<double, kAminoAcidCount> p{};
    std::string source;
};

// Robinson & Robinson 1991 (PNAS 88:8880-8884) composition, which is the
// background NCBI BLAST uses for BLOSUM62 statistics. Sums to 1 to within 1e-12.
[[nodiscard]] const BackgroundFrequencies& robinsonBackground();

struct ProfileOptions {
    // Dirichlet pseudocount added to every observed weighted count. 1/20 per
    // residue keeps a fully conserved column from producing an infinite log-odds
    // and keeps SIFT's ratio finite for an unobserved residue.
    double pseudocount = 0.05;
    // Two sequences more identical than this contribute as a fraction of one
    // sequence each (position-independent Henikoff-style weighting is NOT used
    // here; this is the simpler and more legible 1/n_cluster scheme used by
    // PSI-BLAST's -use_sw_tback clustering and by PROVEAN's 0.75 clustering).
    double clusterIdentity = 0.80;
    int    minimumHomologs = kMinimumHomologs;
    std::string source = "user-supplied homolog set";
};

// Builds the profile. `homologs` are unaligned sequences; each is aligned to the
// query with the tree's Gotoh aligner and projected onto QUERY columns, so the
// profile has exactly query.size() columns and column i is query position i+1.
//
// Throws nothing: a missing matrix, an empty query or too few homologs all come
// back as usable == false with the reason in `warnings`.
[[nodiscard]] ConservationProfile buildProfile(const std::string& queryId,
                                               const std::string& query,
                                               const std::vector<std::string>& homologs,
                                               const SubstitutionMatrix& matrix,
                                               const ProfileOptions& opts = {});

// Shannon entropy in bits over a 20-vector of frequencies. A perfectly conserved
// column is exactly 0; a uniform column is exactly log2(20) = 4.321928...
// ConservationColumn::shannonEntropy is the entropy of the OBSERVED weighted
// distribution with no pseudocount, so those two exact values survive into the
// profile; ConservationColumn::frequencies is the pseudocounted posterior, which
// is what a log-odds column and a SIFT ratio must be computed from.
[[nodiscard]] double shannonEntropyBits(const std::vector<double>& freq);

// SIFT-style tolerance index: p(mutant) / max_y p(y) over the column's weighted,
// pseudocounted frequencies. By construction it is exactly 1.0 when the mutant IS
// the column's most frequent residue, which for a conserved column is the wild
// type. Deleterious below 0.05 (Ng & Henikoff 2003).
[[nodiscard]] double siftScore(const ConservationColumn& col, char mutant);

// PROVEAN-style delta alignment score. PROVEAN's delta is the change in the
// summed alignment score of the query against every supporting sequence when the
// one position is substituted, averaged over the sequences. For a single point
// substitution under a FIXED alignment that sum telescopes exactly: every column
// except this one is unchanged, and a gap in this column contributes zero because
// the gap penalty does not depend on which residue the query has. So the delta is
// mean over non-gap homologs of score(mutant, h) - score(wildType, h), which is
// computed here rather than re-running |homologs| alignments for one substitution.
// Negative means the variant is less compatible with the homolog set; PROVEAN's
// default deleterious threshold is -2.282 (Choi et al. 2012, PLoS ONE 7:e46688;
// balanced accuracy 79.05% on their 58,684-variant human validation set).
[[nodiscard]] double proveanDelta(const std::vector<char>& columnResidues, char wildType,
                                  char mutant, const SubstitutionMatrix& matrix);

// The same delta taken as the expectation over a profile column's distribution
// rather than over a raw residue list: sum_y f_y * (score(mutant, y) -
// score(wildType, y)). It differs from the list form only in that near-duplicate
// homologs have already been down-weighted and the 0.05 pseudocount is included,
// both of which are stated wherever the number is rendered.
[[nodiscard]] double proveanDelta(const ConservationColumn& col, char wildType, char mutant,
                                  const SubstitutionMatrix& matrix);

// The per-homolog residues aligned to a query column, in homolog order. Empty
// when the position is out of range. Kept beside the profile rather than inside
// ConservationColumn because the DTO is frozen.
struct ColumnObservations {
    std::vector<std::vector<char>> perColumn;   // [column][homolog]
};

// Same alignment work as buildProfile, exposed so a caller that needs the raw
// column residues (proveanDelta does) does not re-align.
[[nodiscard]] ColumnObservations observeColumns(const std::string& query,
                                                const std::vector<std::string>& homologs,
                                                const SubstitutionMatrix& matrix,
                                                GapCost gaps = {});

// BLOSUM62 delta for a substitution: score(wt, mut) - score(wt, wt). An exact
// table lookup difference, hence Provenance::Measured at the call site.
[[nodiscard]] int blosum62Delta(const SubstitutionMatrix& matrix, char wildType, char mutant);

}  // namespace biocad::bio
