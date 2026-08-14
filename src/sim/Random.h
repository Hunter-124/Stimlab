// sim/Random.h - the one seeded random stream, its inverse normal CDF, and the
// correlated-sampling machinery built on them.
//
// REPRODUCIBILITY CONTRACT, load-bearing: for a given seed and a given sequence of
// calls, every generator here produces a BYTE-IDENTICAL stream on every platform,
// compiler and standard library. That is why std::mt19937 is not used through
// std::normal_distribution and std::uniform_real_distribution: the distributions
// are not algorithmically specified by the standard, so libstdc++ and MSVC return
// different numbers from the same engine, and a percentile band or a CI-coverage
// figure computed through them silently changes when the toolchain changes. A
// number that cannot be reproduced cannot be checked, so the algorithm is pinned
// here instead: PCG64-DXSM for uniforms and Wichura AS241 for normals.
//
// PCG64-DXSM is O'Neill's 2019 "double xorshift multiply" output permutation over a
// 128-bit LCG state, the generator NumPy adopted as its recommended bit generator.
// AS241 is Wichura's algorithm for the inverse normal CDF, accurate to about 1e-16
// in the central region and to full double precision in the tails BioCAD samples.
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>
#include "core/Uint128.h"

namespace biocad::sim {

// The 128-bit-state PCG64-DXSM bit generator.
//
// State advance is the standard LCG state = state*M + increment over unsigned
// 128-bit arithmetic with the PCG multiplier; the increment is forced odd, which is
// what guarantees the full 2^128 period. The DXSM output function folds the high
// half of the state through a multiply and two xorshifts, which is what removes the
// high-bit correlation a raw truncated LCG has.
class Pcg64Dxsm {
public:
    // Seeded from a 128-bit (seed, sequence) pair exactly as PCG specifies: the
    // state is zeroed, the stream increment is set from `sequence`, one step is
    // taken, the seed is added, and one more step is taken.
    explicit Pcg64Dxsm(std::uint64_t seedLow, std::uint64_t seedHigh = 0,
                       std::uint64_t sequence = 0xda3e39cb94b95bdbULL);

    // Construct directly from the raw 128-bit LCG state and increment. This exists
    // so the known-answer test can start where the reference implementation
    // (NumPy's PCG64DXSM) was placed: seeding schemes differ between libraries,
    // the state transition and output function do not, and it is those that the
    // reproducibility contract is about.
    static Pcg64Dxsm fromRawState(std::uint64_t stateHigh, std::uint64_t stateLow,
                                  std::uint64_t incHigh, std::uint64_t incLow);

    // The next 64 raw bits.
    std::uint64_t nextUInt64();

    // A double in [0, 1) built from the top 53 bits, which is the only division-free
    // construction that is exactly uniform over the representable doubles in the
    // interval.
    double nextDouble();

    // Alias for nextDouble(); the two names exist because a caller sampling a
    // uniform variate and a caller consuming raw randomness read differently.
    double nextUniform() { return nextDouble(); }

    // A standard normal variate by inverse transform through AS241. Inverse
    // transform, not Box-Muller or ziggurat, because it consumes exactly one
    // uniform per variate: that makes a stream position depend only on how many
    // variates were drawn, so a simulation stays reproducible when a caller changes
    // which variates it uses.
    double nextNormal();

    // An unbiased integer in [0, bound) by Lemire-style rejection. Returns 0 for a
    // zero bound rather than dividing by zero.
    std::uint64_t nextBounded(std::uint64_t bound);

private:
    core::Uint128 state_{};
    core::Uint128 increment_{0, 1};
};

// Wichura's AS241 PPND16: the inverse of the standard normal CDF.
// Returns -inf at p <= 0 and +inf at p >= 1 rather than throwing, so a caller that
// clamps its own probabilities is not forced into exception handling.
double inverseNormalCdf(double p);

// Cholesky (LLT) factor of a symmetric positive-definite covariance, row-major in
// and row-major lower-triangular out. Returns false when the matrix is not positive
// definite - an Omega that is not a covariance is a user error to report, not one to
// repair by nudging the diagonal, because the repaired matrix is no longer the
// variability that was entered.
bool choleskyLower(const std::vector<double>& covRowMajor, std::size_t n,
                   std::vector<double>& lowerRowMajor);

// One draw of eta ~ MVN(0, Omega) given the lower Cholesky factor: eta = L * z with
// z standard normal. `out` is resized to n.
void multivariateNormal(const std::vector<double>& lowerRowMajor, std::size_t n,
                        Pcg64Dxsm& rng, std::vector<double>& out);

// Latin hypercube sample of `samples` points in `dimensions`, returned row-major
// (row = sample, column = dimension) on the open unit hypercube. Each dimension is
// stratified into `samples` equal bins with one point per bin, then permuted
// independently, which is what makes LHS cover a marginal in N points as well as
// Monte Carlo covers it in many more.
std::vector<double> latinHypercube(std::size_t samples, std::size_t dimensions,
                                   Pcg64Dxsm& rng);

// Iman-Conover rank-correlation induction, in place, on a row-major
// samples x dimensions matrix. `targetCorrRowMajor` is the desired Spearman rank
// correlation, dimensions x dimensions.
//
// The method rearranges the EXISTING values within each column, so every marginal
// distribution is preserved exactly; only the pairing between columns changes. That
// is the property that lets a Latin hypercube keep its stratification while
// acquiring the correlation structure of a parameter covariance.
bool imanConover(std::vector<double>& samplesRowMajor, std::size_t samples,
                 std::size_t dimensions, const std::vector<double>& targetCorrRowMajor);

// Spearman rank correlation between two equal-length samples, exposed because the
// only honest way to report an induced correlation is to measure it afterwards.
double spearman(const std::vector<double>& a, const std::vector<double>& b);

}  // namespace biocad::sim
