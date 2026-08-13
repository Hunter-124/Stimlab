// bio/Align.h - Gotoh affine-gap sequence alignment and Karlin-Altschul statistics.
//
// Two entry points, deliberately with DIFFERENT return types:
//   alignGlobal() - Needleman-Wunsch, end-to-end, returns GlobalAlignment
//   alignLocal()  - Smith-Waterman, best subsegment, returns LocalAlignment
// Only LocalAlignment can carry an E-value. An E-value is the expected number of
// distinct local hits scoring at least S by chance; a global alignment has no
// such null distribution, so an E-value on one is meaningless. That is enforced
// here by structure, not by a runtime warning: GlobalAlignment has no evalue
// field and evalueOf() only overloads on LocalAlignment.
//
// Substitution matrices are DATA, loaded from assets/packs/matrices/*.json
// (see blosum62.json), so BLOSUM45/80 or a custom matrix need no rebuild.
#pragma once

#include <cstddef>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include <nlohmann/json.hpp>

namespace biocad::bio {

// --------------------------------------------------------------------------
// Substitution matrix

// Karlin-Altschul parameters for one gap-cost pair. These are NOT derivable at
// runtime for gapped alignment - they are fitted constants published per matrix
// per gap cost (NCBI blast_stat.c), so they travel with the matrix file.
struct KarlinAltschul {
    int    gapOpen = 11;
    int    gapExtend = 1;
    double lambda = 0.267;
    double K = 0.041;
    double H = 0.14;
    double alpha = 1.9;
    double beta = -30.0;
    std::string source;
};

// Scores are indexed by ASCII character so lookup is a table hit, not a search.
// Unknown residues fall back to the matrix's own 'X' row when it has one, which
// is what BLAST does; otherwise to `unknownScore`.
class SubstitutionMatrix {
public:
    [[nodiscard]] int score(char a, char b) const;
    [[nodiscard]] bool knows(char c) const { return known_[index(c)]; }

    [[nodiscard]] const std::string& id() const { return id_; }
    [[nodiscard]] const std::string& source() const { return source_; }
    [[nodiscard]] const std::string& alphabet() const { return alphabet_; }
    [[nodiscard]] int defaultGapOpen() const { return defaultGapOpen_; }
    [[nodiscard]] int defaultGapExtend() const { return defaultGapExtend_; }

    // The statistics row matching these gap costs, or nullptr when the matrix
    // file does not publish one. Missing statistics means no E-value, never a
    // guessed lambda.
    [[nodiscard]] const KarlinAltschul* statisticsFor(int gapOpen, int gapExtend) const;

    friend SubstitutionMatrix parseSubstitutionMatrix(const nlohmann::json& j);

private:
    static constexpr int kSlots = 128;
    static std::size_t index(char c);

    std::string id_, source_, alphabet_;
    int defaultGapOpen_ = 11;
    int defaultGapExtend_ = 1;
    int unknownScore_ = -4;
    bool hasX_ = false;
    std::vector<int> table_ = std::vector<int>(kSlots * kSlots, -4);
    std::vector<char> known_ = std::vector<char>(kSlots, 0);
    std::vector<KarlinAltschul> stats_;
};

// Both throw core::Error on a malformed or wrong-schema document.
SubstitutionMatrix parseSubstitutionMatrix(const nlohmann::json& j);
SubstitutionMatrix loadSubstitutionMatrix(const std::filesystem::path& file);

// --------------------------------------------------------------------------
// Alignment

// NCBI gap-cost convention: a gap of length k costs `open + k * extend`, so with
// the BLASTP BLOSUM62 defaults below a length-1 gap costs 12, not 11.
struct GapCost {
    int open = 11;     // BLASTP default for BLOSUM62
    int extend = 1;    // BLASTP default for BLOSUM62
};

// Two numbers that are routinely conflated, and are defined here once:
//   percentIdentity   = 100 * identical / alignedColumns
//   percentSimilarity = 100 * positive  / alignedColumns
// where `alignedColumns` counts only columns with a residue on BOTH sides (gap
// columns are excluded from both denominators), `identical` counts columns whose
// two residues are the same character, and `positive` counts columns whose
// substitution score is > 0 - which includes every identical column plus the
// conservative substitutions. Similarity is therefore always >= identity, and
// quoting similarity as if it were identity oversells the alignment.
struct AlignmentStats {
    std::size_t alignedColumns = 0;   // columns with residues on both sides
    std::size_t identical = 0;
    std::size_t positive = 0;
    std::size_t gapColumns = 0;       // columns containing a gap character
    std::size_t gapOpens = 0;         // number of gap runs, both sequences
    double percentIdentity = 0.0;
    double percentSimilarity = 0.0;
};

// Gapped rows, equal length, '-' for gaps, plus the conventional midline
// ('|' identical, '+' positive-but-different, ' ' otherwise).
struct AlignedRows {
    std::string a;
    std::string b;
    std::string midline;
};

struct GlobalAlignment {
    int            score = 0;
    AlignedRows    rows;
    AlignmentStats stats;
    // No E-value here, by design. See the file header.
};

struct LocalAlignment {
    int            score = 0;
    AlignedRows    rows;
    AlignmentStats stats;
    // Half-open [begin, end) spans in the ORIGINAL sequences, 0-based.
    std::size_t    aBegin = 0, aEnd = 0;
    std::size_t    bBegin = 0, bEnd = 0;
};

GlobalAlignment alignGlobal(std::string_view a, std::string_view b,
                            const SubstitutionMatrix& m, GapCost gaps = {});
LocalAlignment  alignLocal(std::string_view a, std::string_view b,
                           const SubstitutionMatrix& m, GapCost gaps = {});

// --------------------------------------------------------------------------
// Karlin-Altschul statistics - local alignments only.

struct Significance {
    double bitScore = 0.0;
    double evalue = 0.0;
    // Effective lengths after the edge-effect correction, so a reader can see
    // what the E-value was actually computed against.
    double effectiveQueryLength = 0.0;
    double effectiveDbLength = 0.0;
    long   lengthAdjustment = 0;
};

// BLAST_ComputeLengthAdjustment (NCBI blast_stat.c): the fixed-point solution of
//   ell = alpha/lambda * log(K * (m - ell) * (n - N*ell)) + beta
// which removes the edge effect - a local alignment cannot start arbitrarily
// close to the end of either sequence, so the raw m*n search space is too large.
long computeLengthAdjustment(const KarlinAltschul& ka, std::size_t queryLength,
                             double dbLength, long dbNumSeqs);

// E = K * m' * n' * exp(-lambda * S). Only overloaded on LocalAlignment.
Significance evalueOf(const LocalAlignment& hit, const KarlinAltschul& ka,
                      std::size_t queryLength, double dbLength, long dbNumSeqs = 1);

}  // namespace biocad::bio
