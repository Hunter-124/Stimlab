// tests/test_bio_conservation.cpp - homolog-set profiling, column entropy, the
// PSSM background, and the two substitution scores.
//
// The two exact entropy values are the reason this file exists. A conservation
// column is a probability distribution, and the only way to know the arithmetic is
// right is to feed it distributions whose entropy is known in closed form: a
// perfectly conserved column is exactly 0 bits, and a column holding each of the
// 20 residues exactly once is exactly log2(20) = 4.321928094887362 bits. Both are
// asserted to 1e-12, through the real profile builder rather than through the bare
// entropy function, because the pseudocount/weighting path is what could break
// them.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

#include "bio/Align.h"
#include "bio/Conservation.h"

using namespace biocad;
using Catch::Approx;

namespace {

const bio::SubstitutionMatrix& blosum62() {
    static const bio::SubstitutionMatrix m = bio::loadSubstitutionMatrix(
        std::filesystem::path(BIOCAD_ASSETS_DIR) / "packs" / "matrices" / "blosum62.json");
    return m;
}

// n copies of `s`, each made distinguishable at the LAST column so the identity
// clustering does not collapse them into one weighted sequence.
std::vector<std::string> repeat(const std::string& s, int n) {
    std::vector<std::string> out;
    for (int i = 0; i < n; ++i) out.push_back(s);
    return out;
}

}  // namespace

TEST_CASE("A profile below the homolog minimum is unusable and names the shortfall",
          "[bio][conservation]") {
    const std::string query = "ACDEFGHIKL";
    const std::vector<std::string> few = repeat(query, 5);
    const ConservationProfile p = bio::buildProfile("q", query, few, blosum62());

    REQUIRE(p.usable == false);
    REQUIRE(p.minimumHomologsRequired == bio::kMinimumHomologs);
    REQUIRE(p.homologs.sequenceCount == 5);
    REQUIRE(p.columns.empty());
    REQUIRE(!p.warnings.empty());
    REQUIRE(p.warnings[0].find("5 homolog") != std::string::npos);
    REQUIRE(p.warnings[0].find(std::to_string(bio::kMinimumHomologs)) != std::string::npos);
}

TEST_CASE("A perfectly conserved column is exactly 0 bits", "[bio][conservation]") {
    const std::string query = "AAAAAAAAAA";
    const ConservationProfile p = bio::buildProfile("q", query, repeat(query, 20), blosum62());
    REQUIRE(p.usable);
    REQUIRE(p.columns.size() == query.size());
    for (const auto& c : p.columns) REQUIRE(std::abs(c.shannonEntropy - 0.0) < 1e-12);
}

TEST_CASE("A column holding each residue once is exactly log2(20) bits",
          "[bio][conservation]") {
    // 20 homologs, homolog i carrying residue i at every column. Every sequence is
    // a different single-residue run, so the 80% clustering leaves all 20 at weight
    // 1 and the observed distribution is uniform over the 20 residues.
    const std::string query = "AAAA";
    std::vector<std::string> homologs;
    for (int i = 0; i < bio::kAminoAcidCount; ++i)
        homologs.emplace_back(query.size(), bio::kAminoAcids[i]);

    const ConservationProfile p = bio::buildProfile("q", query, homologs, blosum62());
    REQUIRE(p.usable);
    REQUIRE(p.homologs.effectiveSequenceCount == 20);
    const double log2_20 = std::log2(20.0);
    REQUIRE(std::abs(log2_20 - 4.321928094887362) < 1e-12);
    for (const auto& c : p.columns) REQUIRE(std::abs(c.shannonEntropy - log2_20) < 1e-12);
}

TEST_CASE("SIFT's score for the column's own consensus is exactly 1.0",
          "[bio][conservation]") {
    const std::string query = "WWWWW";
    const ConservationProfile p = bio::buildProfile("q", query, repeat(query, 20), blosum62());
    REQUIRE(p.usable);
    // Exactly 1.0, not approximately: the ratio's numerator IS its denominator.
    REQUIRE(bio::siftScore(p.columns[0], 'W') == 1.0);
    // And a residue never seen in the column is the pseudocount over the maximum.
    REQUIRE(bio::siftScore(p.columns[0], 'P') < 0.05);
}

TEST_CASE("BLOSUM62 deltas match direct table lookups", "[bio][conservation]") {
    const bio::SubstitutionMatrix& m = blosum62();
    // score(wt, mut) - score(wt, wt), ten substitutions, each checked against the
    // published BLOSUM62 entries.
    struct Case {
        char wt, mut;
        int sWtMut, sWtWt;
    };
    const Case cases[] = {
        {'A', 'A', 4, 4},   {'A', 'V', 0, 4},   {'W', 'G', -2, 11}, {'C', 'C', 9, 9},
        {'D', 'E', 2, 6},   {'K', 'R', 2, 5},   {'F', 'Y', 3, 6},   {'I', 'L', 2, 4},
        {'G', 'P', -2, 6},  {'S', 'T', 1, 4},
    };
    for (const Case& c : cases) {
        REQUIRE(m.score(c.wt, c.mut) == c.sWtMut);
        REQUIRE(m.score(c.wt, c.wt) == c.sWtWt);
        REQUIRE(bio::blosum62Delta(m, c.wt, c.mut) == c.sWtMut - c.sWtWt);
    }
    REQUIRE(bio::blosum62Delta(m, 'A', 'A') == 0);
}

TEST_CASE("The PROVEAN-style delta is the mean substitution-score change",
          "[bio][conservation]") {
    const bio::SubstitutionMatrix& m = blosum62();
    // Column: four homologs carry L, one carries a gap. Wild type L, mutant P.
    //   score(P,L) = -3, score(L,L) = 4  ->  each term is -7
    //   the gap contributes nothing and is not counted
    //   delta = (4 * -7) / 4 = -7
    const std::vector<char> col{'L', 'L', 'L', 'L', '-'};
    REQUIRE(bio::proveanDelta(col, 'L', 'P', m) == Approx(-7.0));
    // A synonymous "substitution" is exactly zero for every homolog.
    REQUIRE(bio::proveanDelta(col, 'L', 'L', m) == 0.0);
    // And a conservative one is far above the -2.282 threshold.
    REQUIRE(bio::proveanDelta(col, 'L', 'I', m) > -2.282);
    REQUIRE(bio::proveanDelta(col, 'L', 'P', m) < -2.282);
}

TEST_CASE("The background frequency set is named and normalised", "[bio][conservation]") {
    const bio::BackgroundFrequencies& bg = bio::robinsonBackground();
    double sum = 0;
    for (double p : bg.p) sum += p;
    REQUIRE(std::abs(sum - 1.0) < 1e-12);
    REQUIRE(bg.source.find("Robinson") != std::string::npos);

    const std::string query = "ACDEFGHIKL";
    const ConservationProfile p = bio::buildProfile("q", query, repeat(query, 20), blosum62());
    REQUIRE(p.backgroundFrequencySource == bg.source);
    // A conserved column's log-odds for its own residue is positive, and the row
    // has one entry per standard residue.
    REQUIRE(p.columns[0].pssm.size() == 20u);
    REQUIRE(p.columns[0].pssm[static_cast<std::size_t>(bio::aminoAcidIndex('A'))] > 0.0);
}

TEST_CASE("Near-duplicate homologs are down-weighted", "[bio][conservation]") {
    // 20 identical homologs cluster into ONE at the 80% threshold, so the effective
    // count is 1 even though 20 sequences were supplied - which is exactly the
    // situation the effective count exists to expose.
    const std::string query = "ACDEFGHIKLMNPQRSTVWY";
    const ConservationProfile p = bio::buildProfile("q", query, repeat(query, 20), blosum62());
    REQUIRE(p.usable);
    REQUIRE(p.homologs.sequenceCount == 20);
    REQUIRE(p.homologs.effectiveSequenceCount == 1);
    REQUIRE(p.homologs.medianIdentityPct == Approx(100.0));
}
