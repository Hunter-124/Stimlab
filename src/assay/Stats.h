// assay/Stats.h - the small descriptive statistics both Dataset and Qc need.
//
// Header-only and inline because these are four-line definitions that would
// otherwise be duplicated in two translation units; sharing them means the QC
// median and the B-score median cannot drift apart.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace biocad::assay {

// 1/Phi^-1(3/4): scales the MAD so that it estimates sigma for normal data, which
// is what makes a robust Z-prime comparable to the classical one.
inline constexpr double kMadToSigma = 1.4826022185056018;

inline double mean(const std::vector<double>& x) {
    if (x.empty()) return std::nan("");
    double s = 0.0;
    for (double v : x) s += v;
    return s / static_cast<double>(x.size());
}

// Sample standard deviation (n-1). One observation has no spread, not a spread of
// zero, so it returns NaN.
inline double stdDev(const std::vector<double>& x) {
    if (x.size() < 2) return std::nan("");
    const double m = mean(x);
    double s = 0.0;
    for (double v : x) s += (v - m) * (v - m);
    return std::sqrt(s / static_cast<double>(x.size() - 1));
}

inline double stdError(const std::vector<double>& x) {
    if (x.size() < 2) return std::nan("");
    return stdDev(x) / std::sqrt(static_cast<double>(x.size()));
}

inline double median(std::vector<double> x) {
    if (x.empty()) return std::nan("");
    const std::size_t n = x.size();
    std::sort(x.begin(), x.end());
    return (n % 2 == 1) ? x[n / 2] : 0.5 * (x[n / 2 - 1] + x[n / 2]);
}

inline double medianAbsoluteDeviation(const std::vector<double>& x) {
    if (x.empty()) return std::nan("");
    const double m = median(x);
    std::vector<double> dev;
    dev.reserve(x.size());
    for (double v : x) dev.push_back(std::abs(v - m));
    return median(std::move(dev));
}

// Linear-interpolation quantile, the R type-7 / numpy default definition. Stated
// explicitly because the Tukey fences move by a few percent between definitions.
inline double quantileType7(std::vector<double> x, double p) {
    if (x.empty()) return std::nan("");
    std::sort(x.begin(), x.end());
    const double h = (static_cast<double>(x.size()) - 1.0) * p;
    const std::size_t lo = static_cast<std::size_t>(std::floor(h));
    const std::size_t hi = static_cast<std::size_t>(std::ceil(h));
    return x[lo] + (h - static_cast<double>(lo)) * (x[hi] - x[lo]);
}

}  // namespace biocad::assay
