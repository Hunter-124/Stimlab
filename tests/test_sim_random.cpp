// Phase 13.1 - the seeded sampler.
//
// The reproducibility contract is the whole point of pinning PCG64-DXSM and AS241
// instead of using <random>'s distributions, so these cases check the STREAM, not
// just the moments: a known-answer vector, byte-identical repeats, and the exact
// quantiles of the inverse normal CDF.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include "sim/Random.h"

using namespace biocad;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("PCG64-DXSM reproduces the reference stream", "[sim][random]") {
    // Reference: NumPy 2.4.6's PCG64DXSM bit generator, random_raw() from the raw
    // state (123 << 64) | 456 with increment (1 << 64) | 1. Seeding schemes differ
    // between libraries; the state transition and the DXSM output function do not,
    // and those are what the contract is about, so the test starts from a placed
    // state rather than from a seed.
    const std::uint64_t expected[6] = {0x3aa7c032c5f73bbcULL, 0x21083f2d00b154beULL,
                                       0x71a2f46a47ace69bULL, 0xcfbe22186a7ff572ULL,
                                       0x194c5b7018a9fb83ULL, 0x4fe9dfe2a59c5d35ULL};
    auto g = sim::Pcg64Dxsm::fromRawState(123, 456, 1, 1);
    for (std::uint64_t want : expected) REQUIRE(g.nextUInt64() == want);
}

TEST_CASE("The same seed gives a byte-identical stream", "[sim][random]") {
    sim::Pcg64Dxsm a(2026), b(2026), c(2027);
    bool differsSomewhere = false;
    for (int i = 0; i < 500; ++i) {
        const double x = a.nextDouble(), y = b.nextDouble(), z = c.nextDouble();
        REQUIRE(x == y);   // bit equality, not approximate
        if (x != z) differsSomewhere = true;
    }
    REQUIRE(differsSomewhere);
}

TEST_CASE("Uniforms are in range and equidistributed", "[sim][random]") {
    sim::Pcg64Dxsm g(7);
    std::vector<int> bins(10, 0);
    double mean = 0;
    const int n = 200000;
    for (int i = 0; i < n; ++i) {
        const double x = g.nextDouble();
        REQUIRE(x >= 0.0);
        REQUIRE(x < 1.0);
        mean += x;
        bins[static_cast<std::size_t>(x * 10)]++;
    }
    REQUIRE_THAT(mean / n, WithinAbs(0.5, 0.005));
    for (int b : bins) REQUIRE(std::abs(b - n / 10) < n / 100);
}

TEST_CASE("AS241 matches published normal quantiles", "[sim][random]") {
    // Reference values from scipy.stats.norm.ppf (SciPy's own implementation of the
    // inverse normal CDF), quoted to full double precision.
    struct Case { double p, q; };
    const Case cases[] = {{0.975, 1.959963984540054},   {0.025, -1.9599639845400545},
                          {0.95, 1.6448536269514722},   {0.05, -1.6448536269514729},
                          {0.25, -0.6744897501960817},  {0.75, 0.6744897501960817},
                          {1e-10, -6.361340902404056},  {0.999999, 4.753424308817087}};
    for (const auto& c : cases)
        REQUIRE_THAT(sim::inverseNormalCdf(c.p), WithinRel(c.q, 1e-14));
    REQUIRE(sim::inverseNormalCdf(0.5) == 0.0);
    REQUIRE(std::isinf(sim::inverseNormalCdf(0.0)));
    REQUIRE(std::isinf(sim::inverseNormalCdf(1.0)));
}

TEST_CASE("A non-positive-definite Omega is rejected, not repaired", "[sim][random]") {
    std::vector<double> lower;
    // Correlation 2 is not a correlation; nudging the diagonal to make it factor
    // would silently change the variability the caller entered.
    REQUIRE_FALSE(sim::choleskyLower({1.0, 2.0, 2.0, 1.0}, 2, lower));
    REQUIRE(sim::choleskyLower({0.09, 0.03, 0.03, 0.04}, 2, lower));
    REQUIRE(lower.size() == 4);
    REQUIRE_THAT(lower[0], WithinAbs(0.3, 1e-12));
    REQUIRE(lower[1] == 0.0);   // strictly lower triangular
}

TEST_CASE("Multivariate normal draws reproduce Omega", "[sim][random]") {
    std::vector<double> lower;
    REQUIRE(sim::choleskyLower({0.09, 0.03, 0.03, 0.04}, 2, lower));
    sim::Pcg64Dxsm g(5);
    std::vector<double> eta;
    double s11 = 0, s22 = 0, s12 = 0;
    const int n = 100000;
    for (int i = 0; i < n; ++i) {
        sim::multivariateNormal(lower, 2, g, eta);
        s11 += eta[0] * eta[0];
        s22 += eta[1] * eta[1];
        s12 += eta[0] * eta[1];
    }
    REQUIRE_THAT(s11 / n, WithinAbs(0.09, 0.004));
    REQUIRE_THAT(s22 / n, WithinAbs(0.04, 0.003));
    REQUIRE_THAT(s12 / n, WithinAbs(0.03, 0.003));
}

TEST_CASE("Latin hypercube fills every stratum exactly once", "[sim][random]") {
    sim::Pcg64Dxsm g(3);
    const std::size_t samples = 100, dims = 3;
    const auto x = sim::latinHypercube(samples, dims, g);
    REQUIRE(x.size() == samples * dims);
    for (std::size_t d = 0; d < dims; ++d) {
        std::vector<int> hit(samples, 0);
        for (std::size_t s = 0; s < samples; ++s) {
            const double v = x[s * dims + d];
            REQUIRE(v > 0.0);
            REQUIRE(v < 1.0);
            hit[static_cast<std::size_t>(v * static_cast<double>(samples))]++;
        }
        for (int h : hit) REQUIRE(h == 1);
    }
}

TEST_CASE("Iman-Conover induces the target rank correlation", "[sim][random]") {
    sim::Pcg64Dxsm g(99);
    const std::size_t n = 2000, k = 3;
    auto x = sim::latinHypercube(n, k, g);
    const std::vector<double> before = x;
    const std::vector<double> target{1.0, 0.7, 0.0,
                                     0.7, 1.0, 0.0,
                                     0.0, 0.0, 1.0};
    REQUIRE(sim::imanConover(x, n, k, target));

    auto column = [&](const std::vector<double>& m, std::size_t c) {
        std::vector<double> v(n);
        for (std::size_t s = 0; s < n; ++s) v[s] = m[s * k + c];
        return v;
    };
    REQUIRE_THAT(sim::spearman(column(x, 0), column(x, 1)), WithinAbs(0.7, 0.02));
    REQUIRE_THAT(sim::spearman(column(x, 0), column(x, 2)), WithinAbs(0.0, 0.02));

    // The induction only re-pairs values, so every marginal survives unchanged.
    for (std::size_t c = 0; c < k; ++c) {
        auto a = column(before, c), b = column(x, c);
        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());
        REQUIRE(a == b);
    }
}
