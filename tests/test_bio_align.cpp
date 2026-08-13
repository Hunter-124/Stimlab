// tests/test_bio_align.cpp - Gotoh alignment, BLOSUM62 as data, Karlin-Altschul
// statistics, and Kabsch superposition.
//
// Every expected number here is either hand-computed from the published BLOSUM62
// matrix (the arithmetic is written out in the comment, so a reader can audit it
// without running anything) or an exact algebraic identity. Two properties get
// their own cases because getting them wrong is how an alignment or a
// superposition is oversold:
//   - identity and similarity are different numbers, and similarity >= identity;
//   - Kabsch without the determinant correction returns a mirror image with a
//     deceptively low RMSD, so a reflected cloud must NOT superpose well.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <filesystem>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

#include "bio/Align.h"
#include "bio/Superpose.h"
#include "core/Error.h"

using namespace biocad;
using biocad::bio::GapCost;
using biocad::bio::Point3;
using Catch::Approx;

namespace {

std::filesystem::path matrixPath() {
    return std::filesystem::path(BIOCAD_ASSETS_DIR) / "packs" / "matrices" / "blosum62.json";
}

const bio::SubstitutionMatrix& blosum62() {
    static const bio::SubstitutionMatrix m = bio::loadSubstitutionMatrix(matrixPath());
    return m;
}

// The 12-mer used throughout. Self scores from BLOSUM62:
//   A 4, C 9, D 6, E 5, F 6, G 6, H 8, I 4, K 5, L 4, M 5, N 6
// which sum to 4+9+6+5+6+6+8+4+5+4+5+6 = 68.
constexpr const char* kTwelve = "ACDEFGHIKLMN";
constexpr int kTwelveSelfScore = 68;

}  // namespace

TEST_CASE("BLOSUM62 loads from the matrix pack", "[bio][align]") {
    const auto& m = blosum62();
    // The three canonical spot checks against the NCBI distribution.
    REQUIRE(m.score('W', 'W') == 11);
    REQUIRE(m.score('A', 'A') == 4);
    REQUIRE(m.score('A', 'W') == -3);
    REQUIRE(m.score('K', 'R') == 2);   // conservative substitution: positive, not identical
    REQUIRE(m.score('N', 'Q') == 0);   // neither identical nor positive

    REQUIRE(m.defaultGapOpen() == 11);      // BLASTP defaults for BLOSUM62
    REQUIRE(m.defaultGapExtend() == 1);
    REQUIRE(m.alphabet().size() == 24);     // 23 symbols plus '*'
    REQUIRE(m.knows('*'));

    // Lower-case input is the same residue, not an unknown one.
    REQUIRE(m.score('w', 'W') == 11);

    SECTION("statistics travel with the matrix, and a missing row yields no E-value") {
        const auto* ka = m.statisticsFor(11, 1);
        REQUIRE(ka != nullptr);
        REQUIRE(ka->lambda == Approx(0.267));
        REQUIRE(ka->K == Approx(0.041));
        REQUIRE(ka->alpha == Approx(1.9));
        REQUIRE(ka->beta == Approx(-30.0));
        REQUIRE(m.statisticsFor(7, 2) == nullptr);
    }

    SECTION("a wrong schema version is a named error, not a silent default") {
        REQUIRE_THROWS_AS(bio::parseSubstitutionMatrix(nlohmann::json::parse(
                              R"({"schemaVersion": 99, "scores": {"A": {"A": 4}}})")),
                          Error);
    }
}

TEST_CASE("Gotoh global alignment reproduces hand-computed scores", "[bio][align]") {
    const auto& m = blosum62();

    SECTION("identical 12-mers: ungapped, score = sum of self scores = 68") {
        const auto g = bio::alignGlobal(kTwelve, kTwelve, m);
        REQUIRE(g.score == kTwelveSelfScore);
        REQUIRE(g.rows.a == kTwelve);
        REQUIRE(g.rows.b == kTwelve);
        REQUIRE(g.rows.midline == "||||||||||||");
        REQUIRE(g.stats.alignedColumns == 12);
        REQUIRE(g.stats.gapColumns == 0);
        REQUIRE(g.stats.percentIdentity == 100.0);
        REQUIRE(g.stats.percentSimilarity == 100.0);
    }

    SECTION("one conservative substitution K->R: 68 - K/K(5) + K/R(2) = 65") {
        const auto g = bio::alignGlobal(kTwelve, "ACDEFGHIRLMN", m);
        REQUIRE(g.score == 65);
        REQUIRE(g.rows.midline == "||||||||+|||");
        // Identity and similarity are DIFFERENT numbers here: 11 of 12 columns are
        // identical, but all 12 score positively, so similarity is 100%.
        REQUIRE(g.stats.identical == 11);
        REQUIRE(g.stats.positive == 12);
        REQUIRE(g.stats.percentIdentity == Approx(100.0 * 11.0 / 12.0));
        REQUIRE(g.stats.percentSimilarity == 100.0);
        REQUIRE(g.stats.percentSimilarity > g.stats.percentIdentity);
    }

    SECTION("one non-conservative substitution N->Q: 68 - N/N(6) + N/Q(0) = 62") {
        const auto g = bio::alignGlobal(kTwelve, "ACDEFGHIKLMQ", m);
        REQUIRE(g.score == 62);
        REQUIRE(g.rows.midline == "||||||||||| ");
        // N/Q scores 0, which is not positive, so the two numbers coincide.
        REQUIRE(g.stats.percentIdentity == Approx(100.0 * 11.0 / 12.0));
        REQUIRE(g.stats.percentSimilarity == Approx(g.stats.percentIdentity));
    }

    SECTION("one deletion: 11 matches (68 - G/G(6) = 62) minus a length-1 gap (11+1) = 50") {
        const auto g = bio::alignGlobal(kTwelve, "ACDEFHIKLMN", m);
        REQUIRE(g.score == 50);
        REQUIRE(g.rows.a == kTwelve);
        REQUIRE(g.rows.b == "ACDEF-HIKLMN");
        REQUIRE(g.stats.gapOpens == 1);
        REQUIRE(g.stats.gapColumns == 1);
        REQUIRE(g.stats.alignedColumns == 11);
        REQUIRE(g.stats.identical == 11);
        REQUIRE(g.stats.percentIdentity == 100.0);   // over ALIGNED columns only
    }

    SECTION("affine gaps: extending one gap beats opening two") {
        // Deleting FG (self 6 + 6 = 12) leaves 68 - 12 = 56 of matches. One
        // length-2 gap costs 11 + 2*1 = 13, so 56 - 13 = 43. Two separate
        // length-1 gaps would cost 12 + 12 = 24, which is worse by 11.
        const auto g = bio::alignGlobal(kTwelve, "ACDEHIKLMN", m);
        REQUIRE(g.score == 43);
        REQUIRE(g.stats.gapOpens == 1);
        REQUIRE(g.stats.gapColumns == 2);
    }
}

TEST_CASE("Smith-Waterman finds a buried local hit with the right span", "[bio][align]") {
    const auto& m = blosum62();
    // The 12-mer sits at offset 5 in a and offset 3 in b, surrounded by unrelated
    // flanks (Q against W scores -2), so the optimal local hit is exactly the
    // 12-mer and its score is the 68 computed above.
    const std::string a = "QQQQQACDEFGHIKLMNQQQQQ";
    const std::string b = "WWWACDEFGHIKLMNWWWWWW";
    const auto l = bio::alignLocal(a, b, m);

    REQUIRE(l.score == kTwelveSelfScore);
    REQUIRE(l.aBegin == 5);
    REQUIRE(l.aEnd == 17);
    REQUIRE(l.bBegin == 3);
    REQUIRE(l.bEnd == 15);
    REQUIRE(l.rows.a == kTwelve);
    REQUIRE(l.rows.b == kTwelve);
    REQUIRE(l.stats.gapColumns == 0);
    REQUIRE(l.stats.percentIdentity == 100.0);

    // The global alignment of the same pair must align or gap the flanks
    // end-to-end and therefore scores much lower. This is why the two entry
    // points are separate calls with separate result types.
    REQUIRE(bio::alignGlobal(a, b, m).score < l.score);

    SECTION("no positively scoring segment means no hit, not a forced one") {
        const auto none = bio::alignLocal("PPPP", "WWWW", m);
        REQUIRE(none.score == 0);
        REQUIRE(none.rows.a.empty());
    }
}

TEST_CASE("Karlin-Altschul E-values are local-only and behave monotonically", "[bio][align]") {
    const auto& m = blosum62();
    const auto* ka = m.statisticsFor(11, 1);
    REQUIRE(ka != nullptr);

    const std::string a = "QQQQQACDEFGHIKLMNQQQQQ";
    const std::string b = "WWWACDEFGHIKLMNWWWWWW";
    const auto hit = bio::alignLocal(a, b, m);

    const auto sig = bio::evalueOf(hit, *ka, a.size(), 2.0e7, 40000);
    REQUIRE(sig.bitScore == Approx((ka->lambda * hit.score - std::log(ka->K)) / std::log(2.0)));
    REQUIRE(sig.evalue > 0.0);
    REQUIRE(sig.evalue == Approx(ka->K * sig.effectiveQueryLength * sig.effectiveDbLength *
                                 std::exp(-ka->lambda * hit.score)));

    // A larger search space is a larger E-value; a stronger hit is a smaller one.
    REQUIRE(bio::evalueOf(hit, *ka, a.size(), 2.0e8, 400000).evalue > sig.evalue);
    auto strong = hit;
    strong.score = 200;
    REQUIRE(bio::evalueOf(strong, *ka, a.size(), 2.0e7, 40000).evalue < sig.evalue);

    SECTION("the edge-effect length adjustment grows with the query and is bounded") {
        // Values from BLAST_ComputeLengthAdjustment against a 2e7-residue,
        // 40000-sequence database: a 22-residue query is too short for the
        // correction to bite, a 300-residue query loses ~100 residues.
        REQUIRE(bio::computeLengthAdjustment(*ka, 22, 2.0e7, 40000) == 0);
        const long adj300 = bio::computeLengthAdjustment(*ka, 300, 2.0e7, 40000);
        const long adj1000 = bio::computeLengthAdjustment(*ka, 1000, 2.0e7, 40000);
        REQUIRE(adj300 > 50);
        REQUIRE(adj300 < 300);
        REQUIRE(adj1000 >= adj300);
        // The effective lengths actually used are reported, never hidden.
        const auto big = bio::evalueOf(hit, *ka, 300, 2.0e7, 40000);
        REQUIRE(big.effectiveQueryLength == Approx(300.0 - static_cast<double>(adj300)));
    }

    SECTION("a global alignment cannot carry an E-value") {
        // Compile-time, not a runtime warning: GlobalAlignment has no evalue
        // field and evalueOf() only overloads on LocalAlignment. Both of these
        // must fail to compile, which is the whole point:
        //   bio::alignGlobal(a, b, m).evalue;
        //   bio::evalueOf(bio::alignGlobal(a, b, m), *ka, a.size(), 2.0e7);
        static_assert(!std::is_invocable_v<decltype(static_cast<bio::Significance (*)(
                                               const bio::LocalAlignment&,
                                               const bio::KarlinAltschul&, std::size_t, double,
                                               long)>(&bio::evalueOf)),
                                           const bio::GlobalAlignment&, const bio::KarlinAltschul&,
                                           std::size_t, double, long>,
                      "an E-value on a global alignment is meaningless and must not compile");
        SUCCEED();
    }
}

TEST_CASE("Kabsch recovers a known rigid motion", "[bio][superpose]") {
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> u(-20.0, 20.0);
    std::vector<Point3> cloud;
    for (int i = 0; i < 40; ++i) cloud.push_back({u(rng), u(rng), u(rng)});

    // 37 degrees about the normalised axis (1,2,3), then a translation.
    const double ang = 37.0 * 3.14159265358979323846 / 180.0;
    const double norm = std::sqrt(14.0);
    const double ax[3] = {1.0 / norm, 2.0 / norm, 3.0 / norm};
    const double c = std::cos(ang), s = std::sin(ang), t = 1.0 - c;
    const double R[9] = {t * ax[0] * ax[0] + c,         t * ax[0] * ax[1] - s * ax[2],
                         t * ax[0] * ax[2] + s * ax[1], t * ax[0] * ax[1] + s * ax[2],
                         t * ax[1] * ax[1] + c,         t * ax[1] * ax[2] - s * ax[0],
                         t * ax[0] * ax[2] - s * ax[1], t * ax[1] * ax[2] + s * ax[0],
                         t * ax[2] * ax[2] + c};
    const Point3 T = {5.0, -3.0, 11.0};

    std::vector<Point3> moved;
    for (const auto& q : cloud) {
        moved.push_back({R[0] * q[0] + R[1] * q[1] + R[2] * q[2] + T[0],
                         R[3] * q[0] + R[4] * q[1] + R[5] * q[2] + T[1],
                         R[6] * q[0] + R[7] * q[1] + R[8] * q[2] + T[2]});
    }

    const auto fit = bio::kabsch(cloud, moved);
    REQUIRE(fit.pairs == cloud.size());
    REQUIRE(fit.rmsd < 1e-10);
    REQUIRE_FALSE(fit.reflectionCorrected);
    for (std::size_t i = 0; i < 9; ++i) REQUIRE(std::fabs(fit.rotation[i] - R[i]) < 1e-10);
    for (std::size_t i = 0; i < 3; ++i) REQUIRE(std::fabs(fit.translation[i] - T[i]) < 1e-10);

    // The reported fit is the mapping, so applying it reproduces the target.
    const auto q0 = bio::applySuperposition(fit, cloud[0]);
    for (std::size_t i = 0; i < 3; ++i) REQUIRE(std::fabs(q0[i] - moved[0][i]) < 1e-10);

    SECTION("a cloud against itself") {
        // rmsdInPlace is exactly zero: it is a sum of exact zero differences.
        // The FITTED rmsd is an SVD result, so it is zero only to rounding - a
        // test that demanded bit-exact zero there would be testing the FPU.
        REQUIRE(bio::rmsdInPlace(cloud, cloud) == 0.0);
        REQUIRE(bio::kabsch(cloud, cloud).rmsd < 1e-12);
    }

    SECTION("a reflected cloud must NOT superpose well") {
        // Negating z is not reachable by any rotation. A determinant-blind Kabsch
        // reports ~0 RMSD here, which is a wrong answer dressed as a perfect fit.
        std::vector<Point3> mirrored;
        for (const auto& q : cloud) mirrored.push_back({q[0], q[1], -q[2]});
        const auto refl = bio::kabsch(cloud, mirrored);
        REQUIRE(refl.reflectionCorrected);
        REQUIRE(refl.rmsd > 1.0);
        // What comes back is a proper rotation, det = +1.
        const auto& r = refl.rotation;
        const double det = r[0] * (r[4] * r[8] - r[5] * r[7]) -
                           r[1] * (r[3] * r[8] - r[5] * r[6]) +
                           r[2] * (r[3] * r[7] - r[4] * r[6]);
        REQUIRE(det == Approx(1.0).epsilon(1e-12));
    }

    SECTION("degenerate input is refused, not fudged") {
        std::vector<Point3> two{{0, 0, 0}, {1, 0, 0}};
        REQUIRE_THROWS_AS(bio::kabsch(two, two), Error);
        std::vector<Point3> three{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
        REQUIRE_THROWS_AS(bio::kabsch(three, two), Error);
    }
}
