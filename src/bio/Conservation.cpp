#include "bio/Conservation.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace biocad::bio {
namespace {

// Robinson & Robinson 1991 (PNAS 88:8880-8884), the composition NCBI BLAST uses
// with BLOSUM62. Counts in their Table 1 order, renormalised here so the vector
// sums to 1 exactly in double arithmetic rather than to 0.99997.
constexpr double kRobinsonCounts[kAminoAcidCount] = {
    35155,   // A
    8669,    // C
    24161,   // D
    28354,   // E
    17367,   // F
    33229,   // G
    9906,    // H
    23161,   // I
    25872,   // K
    40625,   // L
    10101,   // M
    20212,   // N
    23435,   // P
    19208,   // Q
    23105,   // R
    32070,   // S
    26311,   // T
    29012,   // V
    5990,    // W
    14488,   // Y
};

double clampProbability(double p) { return p > 0.0 ? p : 0.0; }

}  // namespace

int aminoAcidIndex(char c) {
    if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
    const char* p = std::strchr(kAminoAcids, c);
    if (!p || c == '\0') return -1;
    return static_cast<int>(p - kAminoAcids);
}

const BackgroundFrequencies& robinsonBackground() {
    static const BackgroundFrequencies bg = [] {
        BackgroundFrequencies b;
        double total = 0;
        for (double c : kRobinsonCounts) total += c;
        for (int i = 0; i < kAminoAcidCount; ++i) b.p[i] = kRobinsonCounts[i] / total;
        b.source = "Robinson & Robinson 1991 (PNAS 88:8880-8884) amino-acid composition, "
                   "the background NCBI BLAST pairs with BLOSUM62";
        return b;
    }();
    return bg;
}

double shannonEntropyBits(const std::vector<double>& freq) {
    double h = 0.0;
    for (double p : freq) {
        if (p <= 0.0) continue;   // 0 log 0 == 0, and log(0) is not evaluated
        h -= p * std::log2(p);
    }
    return h;
}

double siftScore(const ConservationColumn& col, char mutant) {
    const int m = aminoAcidIndex(mutant);
    if (m < 0 || static_cast<int>(col.frequencies.size()) != kAminoAcidCount) return 0.0;
    const double best = *std::max_element(col.frequencies.begin(), col.frequencies.end());
    if (best <= 0.0) return 0.0;
    return clampProbability(col.frequencies[static_cast<std::size_t>(m)]) / best;
}

double proveanDelta(const std::vector<char>& columnResidues, char wildType, char mutant,
                    const SubstitutionMatrix& matrix) {
    double sum = 0.0;
    int n = 0;
    for (char h : columnResidues) {
        if (h == '-' || aminoAcidIndex(h) < 0) continue;
        sum += static_cast<double>(matrix.score(mutant, h) - matrix.score(wildType, h));
        ++n;
    }
    if (n == 0) return 0.0;
    return sum / static_cast<double>(n);
}

double proveanDelta(const ConservationColumn& col, char wildType, char mutant,
                    const SubstitutionMatrix& matrix) {
    if (static_cast<int>(col.frequencies.size()) != kAminoAcidCount) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < kAminoAcidCount; ++i) {
        const char y = kAminoAcids[i];
        sum += col.frequencies[static_cast<std::size_t>(i)] *
               static_cast<double>(matrix.score(mutant, y) - matrix.score(wildType, y));
    }
    return sum;
}

int blosum62Delta(const SubstitutionMatrix& matrix, char wildType, char mutant) {
    return matrix.score(wildType, mutant) - matrix.score(wildType, wildType);
}

ColumnObservations observeColumns(const std::string& query,
                                  const std::vector<std::string>& homologs,
                                  const SubstitutionMatrix& matrix, GapCost gaps) {
    ColumnObservations obs;
    obs.perColumn.assign(query.size(), {});
    for (const std::string& h : homologs) {
        const GlobalAlignment g = alignGlobal(query, h, matrix, gaps);
        std::size_t qcol = 0;
        for (std::size_t i = 0; i < g.rows.a.size() && qcol < query.size(); ++i) {
            const char qa = g.rows.a[i];
            if (qa == '-') continue;   // an insertion in the homolog: no query column
            obs.perColumn[qcol].push_back(g.rows.b[i]);
            ++qcol;
        }
        // A truncated homolog leaves trailing query columns unobserved; record the
        // gap explicitly so every column vector has one entry per homolog and the
        // gap fraction is right.
        for (; qcol < query.size(); ++qcol) obs.perColumn[qcol].push_back('-');
    }
    return obs;
}

ConservationProfile buildProfile(const std::string& queryId, const std::string& query,
                                 const std::vector<std::string>& homologs,
                                 const SubstitutionMatrix& matrix,
                                 const ProfileOptions& opts) {
    ConservationProfile out;
    out.queryId = queryId;
    out.query = query;
    out.minimumHomologsRequired = opts.minimumHomologs;
    out.backgroundFrequencySource = robinsonBackground().source;
    out.homologs.sequenceCount = static_cast<int>(homologs.size());
    out.homologs.source = opts.source;
    out.usable = false;

    if (query.empty()) {
        out.warnings.push_back("The query sequence is empty, so there is nothing to profile.");
        return out;
    }

    const ColumnObservations obs = observeColumns(query, homologs, matrix, GapCost{});
    out.homologs.alignmentColumns = static_cast<int>(query.size());

    // Per-homolog identity to the query, over columns where both have a residue.
    // These are the numbers a reader needs to judge the set, so they are reported
    // even when the set is too small to score with.
    std::vector<double> identities;
    identities.reserve(homologs.size());
    for (std::size_t h = 0; h < homologs.size(); ++h) {
        std::size_t both = 0, same = 0;
        for (std::size_t c = 0; c < query.size(); ++c) {
            const char hr = obs.perColumn[c][h];
            if (hr == '-') continue;
            ++both;
            if (hr == query[c]) ++same;
        }
        identities.push_back(both ? 100.0 * static_cast<double>(same) /
                                        static_cast<double>(both)
                                  : 0.0);
    }
    if (!identities.empty()) {
        std::vector<double> sorted = identities;
        std::sort(sorted.begin(), sorted.end());
        const std::size_t n = sorted.size();
        out.homologs.medianIdentityPct =
            (n % 2) ? sorted[n / 2] : 0.5 * (sorted[n / 2 - 1] + sorted[n / 2]);
        out.homologs.minIdentityPct = sorted.front();
        out.homologs.maxIdentityPct = sorted.back();
    }

    // Sequence weighting. Homologs whose pairwise identity exceeds the threshold
    // are single-linkage clustered and each member carries 1/clusterSize, so a
    // set padded with near-duplicates does not look like independent evidence.
    const std::size_t n = homologs.size();
    std::vector<int> cluster(n, -1);
    int clusters = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (cluster[i] >= 0) continue;
        cluster[i] = clusters++;
        for (std::size_t j = i + 1; j < n; ++j) {
            if (cluster[j] >= 0) continue;
            std::size_t both = 0, same = 0;
            for (std::size_t c = 0; c < query.size(); ++c) {
                const char a = obs.perColumn[c][i], b = obs.perColumn[c][j];
                if (a == '-' || b == '-') continue;
                ++both;
                if (a == b) ++same;
            }
            const double id = both ? static_cast<double>(same) / static_cast<double>(both) : 0.0;
            if (id >= opts.clusterIdentity) cluster[j] = cluster[i];
        }
    }
    std::vector<double> clusterSize(static_cast<std::size_t>(clusters), 0.0);
    for (int c : cluster) clusterSize[static_cast<std::size_t>(c)] += 1.0;
    std::vector<double> weight(n, 1.0);
    double effective = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        weight[i] = 1.0 / clusterSize[static_cast<std::size_t>(cluster[i])];
        effective += weight[i];
    }
    out.homologs.effectiveSequenceCount = static_cast<int>(std::lround(effective));

    if (static_cast<int>(n) < opts.minimumHomologs) {
        out.warnings.push_back(
            "Only " + std::to_string(n) + " homolog(s) were supplied; " +
            std::to_string(opts.minimumHomologs) +
            " is the minimum this build will score from. Below it a column "
            "distribution is dominated by the pseudocount rather than by the "
            "protein, so no conservation score is produced. The minimum is a "
            "BioCAD judgement, not a published constant.");
        return out;
    }

    const BackgroundFrequencies& bg = robinsonBackground();
    out.columns.reserve(query.size());
    for (std::size_t c = 0; c < query.size(); ++c) {
        ConservationColumn col;
        col.position = static_cast<int>(c) + 1;
        col.queryResidue = query[c];

        std::vector<double> counts(kAminoAcidCount, 0.0);
        double observed = 0.0, gapWeight = 0.0, totalWeight = 0.0;
        for (std::size_t h = 0; h < n; ++h) {
            const char r = obs.perColumn[c][h];
            totalWeight += weight[h];
            const int idx = aminoAcidIndex(r);
            if (idx < 0) {
                gapWeight += weight[h];
                continue;
            }
            counts[static_cast<std::size_t>(idx)] += weight[h];
            observed += weight[h];
        }
        col.gapFraction = totalWeight > 0 ? gapWeight / totalWeight : 0.0;

        // Entropy comes from the OBSERVED distribution, with no pseudocount, so a
        // perfectly conserved column is exactly 0 bits and a column holding each of
        // the 20 residues equally is exactly log2(20).
        std::vector<double> observedFreq(kAminoAcidCount, 0.0);
        if (observed > 0) {
            for (int i = 0; i < kAminoAcidCount; ++i)
                observedFreq[static_cast<std::size_t>(i)] =
                    counts[static_cast<std::size_t>(i)] / observed;
        }
        col.shannonEntropy = shannonEntropyBits(observedFreq);

        // The pseudocounted posterior is what SIFT and the log-odds column use.
        const double denom = observed + opts.pseudocount * kAminoAcidCount;
        col.frequencies.assign(kAminoAcidCount, 0.0);
        col.pssm.assign(kAminoAcidCount, 0.0);
        for (int i = 0; i < kAminoAcidCount; ++i) {
            const double p = (counts[static_cast<std::size_t>(i)] + opts.pseudocount) / denom;
            col.frequencies[static_cast<std::size_t>(i)] = p;
            col.pssm[static_cast<std::size_t>(i)] = std::log2(p / bg.p[static_cast<std::size_t>(i)]);
        }
        out.columns.push_back(std::move(col));
    }

    out.usable = true;
    out.warnings.push_back(
        "Conservation is estimated from " + std::to_string(n) + " user-supplied homolog(s) "
        "(effective count " + std::to_string(out.homologs.effectiveSequenceCount) +
        " after clustering at " +
        std::to_string(static_cast<int>(std::lround(opts.clusterIdentity * 100))) +
        "% identity). It describes THIS set, not the protein family in general.");
    return out;
}

}  // namespace biocad::bio
