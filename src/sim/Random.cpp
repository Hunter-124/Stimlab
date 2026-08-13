#include "sim/Random.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

#include <Eigen/Cholesky>
#include <Eigen/Core>

namespace biocad::sim {
namespace {

// The PCG "cheap multiplier": a 64-bit constant multiplied into the 128-bit state.
// The DXSM variant uses this rather than the full 128-bit multiplier because its
// output permutation, not the state transition, is what supplies the statistical
// quality - and a 64x128 multiply is one instruction pair instead of four.
constexpr std::uint64_t kCheapMultiplier = 0xda942042e4dd58b5ULL;

inline __uint128_t pack128(std::uint64_t hi, std::uint64_t lo) {
    return (static_cast<__uint128_t>(hi) << 64) | lo;
}

}  // namespace

Pcg64Dxsm::Pcg64Dxsm(std::uint64_t seedLow, std::uint64_t seedHigh, std::uint64_t sequence) {
    // The increment must be odd for the LCG to have full period; PCG guarantees that
    // by shifting the stream selector left one bit and setting bit 0.
    increment_ = (pack128(sequence, sequence ^ 0x9e3779b97f4a7c15ULL) << 1) | 1;
    state_ = 0;
    state_ = state_ * kCheapMultiplier + increment_;
    state_ += pack128(seedHigh, seedLow);
    state_ = state_ * kCheapMultiplier + increment_;
}

Pcg64Dxsm Pcg64Dxsm::fromRawState(std::uint64_t stateHigh, std::uint64_t stateLow,
                                  std::uint64_t incHigh, std::uint64_t incLow) {
    Pcg64Dxsm g(0);
    g.state_ = pack128(stateHigh, stateLow);
    g.increment_ = pack128(incHigh, incLow);
    return g;
}

std::uint64_t Pcg64Dxsm::nextUInt64() {
    // DXSM output on the CURRENT state, then advance. The order matters for stream
    // identity with the reference implementation, not for quality.
    std::uint64_t hi = static_cast<std::uint64_t>(state_ >> 64);
    std::uint64_t lo = static_cast<std::uint64_t>(state_) | 1ULL;
    hi ^= hi >> 32;
    hi *= kCheapMultiplier;
    hi ^= hi >> 48;
    hi *= lo;
    state_ = state_ * kCheapMultiplier + increment_;
    return hi;
}

double Pcg64Dxsm::nextDouble() {
    // 53 bits is exactly the mantissa width, so every representable double in the
    // interval is hit with equal probability and no rounding can produce 1.0.
    return static_cast<double>(nextUInt64() >> 11) * (1.0 / 9007199254740992.0);
}

double Pcg64Dxsm::nextNormal() {
    // The uniform is shifted off the open-interval endpoints: AS241 returns an
    // infinity at exactly 0, and a single infinite variate would poison a whole
    // population band.
    double u = nextDouble();
    if (u <= 0.0) u = std::numeric_limits<double>::min();
    return inverseNormalCdf(u);
}

std::uint64_t Pcg64Dxsm::nextBounded(std::uint64_t bound) {
    if (bound == 0) return 0;
    // Lemire's multiply-shift with rejection: unbiased, and it rejects only when the
    // 64-bit draw falls in the short residual window, so the stream advance is
    // deterministic given the draws.
    const std::uint64_t threshold = (-bound) % bound;
    while (true) {
        const std::uint64_t r = nextUInt64();
        if (r >= threshold) return r % bound;
    }
}

double inverseNormalCdf(double p) {
    if (p <= 0.0) return -std::numeric_limits<double>::infinity();
    if (p >= 1.0) return std::numeric_limits<double>::infinity();

    // Wichura (1988), Algorithm AS 241, PPND16. Three rational approximations in
    // |q| < 0.425, in the shoulders, and in the far tail; the published coefficients
    // are reproduced exactly, because a truncated coefficient turns a 1e-16 method
    // into a 1e-9 one.
    const double q = p - 0.5;
    double r = 0.0;
    if (std::fabs(q) <= 0.425) {
        r = 0.180625 - q * q;
        return q *
               (((((((2.5090809287301226727e+3 * r + 3.3430575583588128105e+4) * r +
                     6.7265770927008700853e+4) *
                        r +
                    4.5921953931549871457e+4) *
                       r +
                   1.3731693765509461125e+4) *
                      r +
                  1.9715909503065514427e+3) *
                     r +
                 1.3314166789178437745e+2) *
                    r +
                3.3871328727963666080e0) /
               (((((((5.2264952788528545610e+3 * r + 2.8729085735721942674e+4) * r +
                     3.9307895800092710610e+4) *
                        r +
                    2.1213794301586595867e+4) *
                       r +
                   5.3941960214247511077e+3) *
                      r +
                  6.8718700749205790830e+2) *
                     r +
                 4.2313330701600911252e+1) *
                    r +
                1.0);
    }

    r = (q < 0) ? p : 1.0 - p;
    r = std::sqrt(-std::log(r));
    double value = 0.0;
    if (r <= 5.0) {
        r -= 1.6;
        value = (((((((7.74545014278341407640e-4 * r + 2.27238449892691845833e-2) * r +
                      2.41780725177450611770e-1) *
                         r +
                     1.27045825245236838258e0) *
                        r +
                    3.64784832476320460504e0) *
                       r +
                   5.76949722146069140550e0) *
                      r +
                  4.63033784615654529590e0) *
                     r +
                 1.42343711074968357734e0) /
                (((((((1.05075007164441684324e-9 * r + 5.47593808499534494600e-4) * r +
                      1.51986665636164571966e-2) *
                         r +
                     1.48103976427480074590e-1) *
                        r +
                    6.89767334985100004550e-1) *
                       r +
                   1.67638483018380384940e0) *
                      r +
                  2.05319162663775882187e0) *
                     r +
                 1.0);
    } else {
        r -= 5.0;
        value = (((((((2.01033439929228813265e-7 * r + 2.71155556874348757815e-5) * r +
                      1.24266094738807843860e-3) *
                         r +
                     2.65321895265761230930e-2) *
                        r +
                    2.96560571828504891230e-1) *
                       r +
                   1.78482653991729133580e0) *
                      r +
                  5.46378491116411436990e0) *
                     r +
                 6.65790464350110377720e0) /
                (((((((2.04426310338993978564e-15 * r + 1.42151175831644588870e-7) * r +
                      1.84631831751005468180e-5) *
                         r +
                     7.86869131145613259100e-4) *
                        r +
                    1.48753612908506148525e-2) *
                       r +
                   1.36929880922735805310e-1) *
                      r +
                  5.99832206555887937690e-1) *
                     r +
                 1.0);
    }
    return (q < 0) ? -value : value;
}

bool choleskyLower(const std::vector<double>& covRowMajor, std::size_t n,
                   std::vector<double>& lowerRowMajor) {
    lowerRowMajor.clear();
    if (n == 0 || covRowMajor.size() != n * n) return false;
    Eigen::MatrixXd m(static_cast<Eigen::Index>(n), static_cast<Eigen::Index>(n));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            m(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
                covRowMajor[i * n + j];
    Eigen::LLT<Eigen::MatrixXd> llt(m);
    if (llt.info() != Eigen::Success) return false;
    const Eigen::MatrixXd L = llt.matrixL();
    lowerRowMajor.resize(n * n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            lowerRowMajor[i * n + j] =
                L(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
    return true;
}

void multivariateNormal(const std::vector<double>& lowerRowMajor, std::size_t n,
                        Pcg64Dxsm& rng, std::vector<double>& out) {
    out.assign(n, 0.0);
    if (lowerRowMajor.size() != n * n) return;
    std::vector<double> z(n);
    for (std::size_t i = 0; i < n; ++i) z[i] = rng.nextNormal();
    for (std::size_t i = 0; i < n; ++i) {
        double acc = 0.0;
        for (std::size_t j = 0; j <= i; ++j) acc += lowerRowMajor[i * n + j] * z[j];
        out[i] = acc;
    }
}

std::vector<double> latinHypercube(std::size_t samples, std::size_t dimensions,
                                   Pcg64Dxsm& rng) {
    std::vector<double> out(samples * dimensions, 0.0);
    if (samples == 0 || dimensions == 0) return out;
    std::vector<std::size_t> perm(samples);
    for (std::size_t d = 0; d < dimensions; ++d) {
        std::iota(perm.begin(), perm.end(), std::size_t{0});
        // Fisher-Yates through the pinned bounded generator, not std::shuffle, whose
        // permutation is implementation-defined for a given engine.
        for (std::size_t i = samples; i > 1; --i) {
            const std::size_t j = static_cast<std::size_t>(rng.nextBounded(i));
            std::swap(perm[i - 1], perm[j]);
        }
        for (std::size_t s = 0; s < samples; ++s) {
            const double u = rng.nextDouble();
            out[s * dimensions + d] =
                (static_cast<double>(perm[s]) + u) / static_cast<double>(samples);
        }
    }
    return out;
}

namespace {

// Fractional ranks with ties averaged, 1-based. Averaging is what makes Spearman on
// tied data equal Pearson on the ranks, which is the identity the induction relies
// on.
std::vector<double> ranks(const std::vector<double>& v) {
    const std::size_t n = v.size();
    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), std::size_t{0});
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return v[a] < v[b]; });
    std::vector<double> r(n, 0.0);
    std::size_t i = 0;
    while (i < n) {
        std::size_t j = i;
        while (j + 1 < n && v[idx[j + 1]] == v[idx[i]]) ++j;
        const double avg = 0.5 * (static_cast<double>(i) + static_cast<double>(j)) + 1.0;
        for (std::size_t k = i; k <= j; ++k) r[idx[k]] = avg;
        i = j + 1;
    }
    return r;
}

double pearson(const std::vector<double>& a, const std::vector<double>& b) {
    const std::size_t n = a.size();
    if (n < 2 || b.size() != n) return 0.0;
    const double ma = std::accumulate(a.begin(), a.end(), 0.0) / static_cast<double>(n);
    const double mb = std::accumulate(b.begin(), b.end(), 0.0) / static_cast<double>(n);
    double sab = 0.0, saa = 0.0, sbb = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double da = a[i] - ma, db = b[i] - mb;
        sab += da * db;
        saa += da * da;
        sbb += db * db;
    }
    if (saa <= 0.0 || sbb <= 0.0) return 0.0;
    return sab / std::sqrt(saa * sbb);
}

}  // namespace

double spearman(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size() || a.size() < 2) return 0.0;
    return pearson(ranks(a), ranks(b));
}

bool imanConover(std::vector<double>& samplesRowMajor, std::size_t samples,
                 std::size_t dimensions, const std::vector<double>& targetCorrRowMajor) {
    if (samples < 3 || dimensions < 2) return false;
    if (samplesRowMajor.size() != samples * dimensions) return false;
    if (targetCorrRowMajor.size() != dimensions * dimensions) return false;

    const Eigen::Index n = static_cast<Eigen::Index>(samples);
    const Eigen::Index k = static_cast<Eigen::Index>(dimensions);

    // Van der Waerden scores: the deterministic normal order statistics
    // Phi^{-1}(i/(N+1)). Iman-Conover works on these rather than on the data because
    // their correlation matrix can be driven to the target exactly by a linear map,
    // and only the resulting RANK pattern is transplanted back onto the data.
    Eigen::VectorXd score(n);
    for (Eigen::Index i = 0; i < n; ++i)
        score(i) = inverseNormalCdf(static_cast<double>(i + 1) / static_cast<double>(n + 1));

    // Each column gets the same scores in an independently permuted order, so the
    // score matrix starts near-uncorrelated and the linear map below has something
    // to work with; a deterministic shift of a monotone sequence would start at
    // correlation ~1 and the induction would have no room to move.
    //
    // The permutation stream is INTERNAL and fixed-seed: the induction must be
    // reproducible, and it must not consume numbers from the caller's stream, or
    // switching correlation on would shift every subsequent variate.
    Pcg64Dxsm perm(0x9e3779b97f4a7c15ULL, 0xbf58476d1ce4e5b9ULL);
    Eigen::MatrixXd m(n, k);
    std::vector<std::size_t> order0(samples);
    for (Eigen::Index c = 0; c < k; ++c) {
        std::iota(order0.begin(), order0.end(), std::size_t{0});
        for (std::size_t i = samples; i > 1; --i)
            std::swap(order0[i - 1], order0[static_cast<std::size_t>(perm.nextBounded(i))]);
        for (Eigen::Index i = 0; i < n; ++i)
            m(i, c) = score(static_cast<Eigen::Index>(order0[static_cast<std::size_t>(i)]));
    }

    // The target is a SPEARMAN rank correlation, but the linear map below controls
    // the PEARSON correlation of the normal scores. For normal scores the two are
    // related exactly by r_s = (6/pi) asin(rho/2), so the target is inverted through
    // rho = 2 sin(pi * r_s / 6) first. Skipping this step is the classic
    // Iman-Conover error: asking for 0.70 and measuring 0.68.
    Eigen::MatrixXd target(k, k);
    for (Eigen::Index i = 0; i < k; ++i) {
        for (Eigen::Index j = 0; j < k; ++j) {
            const double rs = targetCorrRowMajor[static_cast<std::size_t>(i) * dimensions +
                                                 static_cast<std::size_t>(j)];
            target(i, j) = (i == j) ? 1.0 : 2.0 * std::sin(M_PI * rs / 6.0);
        }
    }

    // E = corr(M), C = target. With E = Q Q' and C = P P', the map M -> M Q'^-1 P'
    // has correlation exactly C in the limit and to within sampling error here.
    Eigen::MatrixXd centered = m.rowwise() - m.colwise().mean();
    Eigen::MatrixXd cov = (centered.transpose() * centered) / static_cast<double>(n - 1);
    Eigen::VectorXd sd = cov.diagonal().array().sqrt();
    Eigen::MatrixXd e(k, k);
    for (Eigen::Index i = 0; i < k; ++i)
        for (Eigen::Index j = 0; j < k; ++j) e(i, j) = cov(i, j) / (sd(i) * sd(j));

    Eigen::LLT<Eigen::MatrixXd> lltE(e);
    Eigen::LLT<Eigen::MatrixXd> lltC(target);
    if (lltE.info() != Eigen::Success || lltC.info() != Eigen::Success) return false;
    const Eigen::MatrixXd q = lltE.matrixL();
    const Eigen::MatrixXd p = lltC.matrixL();
    // M Q'^-1 by triangular back-substitution rather than an explicit inverse: Q is
    // a Cholesky factor, so the solve is exact to working precision and cheaper.
    Eigen::MatrixXd mq = m;
    q.transpose().triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheRight>(mq);
    const Eigen::MatrixXd t = mq * p.transpose();   // n x k target-rank pattern

    // Transplant the ranks: sort each data column and place its s-th smallest value
    // at the row whose T-rank is s. The marginal is therefore preserved exactly -
    // no value is altered, only re-paired.
    for (std::size_t c = 0; c < dimensions; ++c) {
        std::vector<double> col(samples), tcol(samples);
        for (std::size_t s = 0; s < samples; ++s) {
            col[s] = samplesRowMajor[s * dimensions + c];
            tcol[s] = t(static_cast<Eigen::Index>(s), static_cast<Eigen::Index>(c));
        }
        std::vector<double> sorted = col;
        std::sort(sorted.begin(), sorted.end());
        std::vector<std::size_t> order(samples);
        std::iota(order.begin(), order.end(), std::size_t{0});
        std::sort(order.begin(), order.end(),
                  [&](std::size_t a, std::size_t b) { return tcol[a] < tcol[b]; });
        for (std::size_t s = 0; s < samples; ++s)
            samplesRowMajor[order[s] * dimensions + c] = sorted[s];
    }
    return true;
}

}  // namespace biocad::sim
