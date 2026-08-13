// Tests for the shared fitter and integrator. Every case here is one that would
// actually catch a bug: exact parameter recovery from noiseless data, refusal to
// invent an error bar, and the integrator checked against a closed form.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#include "numeric/Ode.h"
#include "numeric/Optimize.h"

using namespace biocad::numeric;
using Catch::Matchers::WithinAbs;

TEST_CASE("levenbergMarquardt recovers a noiseless linear model", "[numeric]") {
    const std::vector<double> xs{0, 1, 2, 3, 4, 5, 6};
    std::vector<double> ys;
    for (double x : xs) {
        ys.push_back(3.0 * x + 7.0);
    }

    auto evaluate = [&](const std::vector<double>& p, std::vector<double>& r,
                        std::vector<double>& j) {
        for (std::size_t i = 0; i < xs.size(); ++i) {
            r[i] = p[0] * xs[i] + p[1] - ys[i];
            j[i * 2 + 0] = xs[i];
            j[i * 2 + 1] = 1.0;
        }
    };

    const LmResult fit = levenbergMarquardt({0.0, 0.0}, xs.size(), evaluate);
    REQUIRE(fit.converged);
    REQUIRE_THAT(fit.params[0], WithinAbs(3.0, 1e-9));
    REQUIRE_THAT(fit.params[1], WithinAbs(7.0, 1e-9));

    // A perfect fit has a vanishing error bar - but the vector must still exist.
    REQUIRE(fit.standardErrors.size() == 2);
    REQUIRE(fit.standardErrors[0] < 1e-6);
    REQUIRE(fit.standardErrors[1] < 1e-6);

    std::vector<double> fitted;
    for (double x : xs) {
        fitted.push_back(fit.params[0] * x + fit.params[1]);
    }
    REQUIRE_THAT(rSquared(ys, fitted), WithinAbs(1.0, 1e-12));
}

TEST_CASE("levenbergMarquardt recovers a noiseless exponential decay", "[numeric]") {
    const double amplitude = 5.0;
    const double rate = 0.35;
    const std::vector<double> ts{0, 0.5, 1, 2, 3, 5, 8, 12};
    std::vector<double> ys;
    for (double t : ts) {
        ys.push_back(amplitude * std::exp(-rate * t));
    }

    auto evaluate = [&](const std::vector<double>& p, std::vector<double>& r,
                        std::vector<double>& j) {
        for (std::size_t i = 0; i < ts.size(); ++i) {
            const double e = std::exp(-p[1] * ts[i]);
            r[i] = p[0] * e - ys[i];
            j[i * 2 + 0] = e;
            j[i * 2 + 1] = -p[0] * ts[i] * e;
        }
    };

    const LmResult fit = levenbergMarquardt({1.0, 1.0}, ts.size(), evaluate);
    REQUIRE(fit.converged);
    REQUIRE_THAT(fit.params[0], WithinAbs(amplitude, 1e-9));
    REQUIRE_THAT(fit.params[1], WithinAbs(rate, 1e-9));
}

TEST_CASE("levenbergMarquardt recovers a noiseless four-parameter logistic", "[numeric]") {
    // Same parameterisation the PD slice fits:
    //   E(x) = Bottom + (Top - Bottom) / (1 + 10^(nH * (log10 EC50 - x))), x = log10[A]
    const double kTop = 100.0;
    const double kBottom = 0.0;
    const double kLogEc50 = -7.0;  // EC50 = 1e-7 M
    const double kSlope = 1.2;

    std::vector<double> xs;
    for (int i = 0; i < 11; ++i) {
        xs.push_back(-10.0 + 0.5 * i);
    }
    std::vector<double> ys;
    for (double x : xs) {
        ys.push_back(kBottom
                     + (kTop - kBottom) / (1.0 + std::pow(10.0, kSlope * (kLogEc50 - x))));
    }

    auto evaluate = [&](const std::vector<double>& p, std::vector<double>& r,
                        std::vector<double>& j) {
        const double ln10 = std::log(10.0);
        for (std::size_t i = 0; i < xs.size(); ++i) {
            const double u = std::pow(10.0, p[3] * (p[2] - xs[i]));
            const double delta = p[0] - p[1];
            const double den = 1.0 + u;
            r[i] = p[1] + delta / den - ys[i];
            j[i * 4 + 0] = 1.0 / den;
            j[i * 4 + 1] = u / den;
            j[i * 4 + 2] = -delta * p[3] * ln10 * u / (den * den);
            j[i * 4 + 3] = -delta * (p[2] - xs[i]) * ln10 * u / (den * den);
        }
    };

    const LmResult fit = levenbergMarquardt({90.0, 5.0, -6.0, 1.0}, xs.size(), evaluate);
    REQUIRE(fit.converged);
    REQUIRE_THAT(fit.params[0], WithinAbs(kTop, 1e-9));
    REQUIRE_THAT(fit.params[1], WithinAbs(kBottom, 1e-9));
    REQUIRE_THAT(fit.params[2], WithinAbs(kLogEc50, 1e-9));
    REQUIRE_THAT(fit.params[3], WithinAbs(kSlope, 1e-9));
    REQUIRE(fit.standardErrors.size() == 4);
    for (double e : fit.standardErrors) {
        REQUIRE(e < 1e-5);
    }
}

TEST_CASE("levenbergMarquardt refuses an underdetermined fit", "[numeric]") {
    auto evaluate = [](const std::vector<double>& p, std::vector<double>& r,
                       std::vector<double>& j) {
        for (std::size_t i = 0; i < r.size(); ++i) {
            r[i] = p[0] - 1.0;
            for (std::size_t k = 0; k < 3; ++k) {
                j[i * 3 + k] = (k == 0) ? 1.0 : 0.0;
            }
        }
    };

    const LmResult fit = levenbergMarquardt({0.0, 0.0, 0.0}, 2, evaluate);
    REQUIRE_FALSE(fit.converged);
    REQUIRE(fit.standardErrors.empty());
    REQUIRE_FALSE(fit.note.empty());
}

TEST_CASE("levenbergMarquardt emits no error bars for a rank-deficient system", "[numeric]") {
    // Two parameters that only ever appear as their sum: the fit still finds a
    // minimum, but the covariance matrix is singular and an error bar would be a
    // fabrication.
    const std::vector<double> xs{0, 1, 2, 3, 4, 5};
    auto evaluate = [&](const std::vector<double>& p, std::vector<double>& r,
                        std::vector<double>& j) {
        for (std::size_t i = 0; i < xs.size(); ++i) {
            r[i] = (p[0] + p[1]) * xs[i] - 2.0 * xs[i];
            j[i * 2 + 0] = xs[i];
            j[i * 2 + 1] = xs[i];
        }
    };

    const LmResult fit = levenbergMarquardt({0.0, 0.0}, xs.size(), evaluate);
    REQUIRE_THAT(fit.params[0] + fit.params[1], WithinAbs(2.0, 1e-9));
    REQUIRE(fit.standardErrors.empty());
    REQUIRE(fit.note.find("rank-deficient") != std::string::npos);
}

TEST_CASE("levenbergMarquardt rejects empty and non-finite input", "[numeric]") {
    auto noop = [](const std::vector<double>&, std::vector<double>&, std::vector<double>&) {};
    const LmResult empty = levenbergMarquardt({1.0}, 0, noop);
    REQUIRE_FALSE(empty.converged);
    REQUIRE(empty.standardErrors.empty());
    REQUIRE_FALSE(empty.note.empty());

    auto nanModel = [](const std::vector<double>&, std::vector<double>& r,
                       std::vector<double>& j) {
        for (std::size_t i = 0; i < r.size(); ++i) {
            r[i] = std::nan("");
            j[i] = 1.0;
        }
    };
    const LmResult bad = levenbergMarquardt({1.0}, 4, nanModel);
    REQUIRE_FALSE(bad.converged);
    REQUIRE(bad.standardErrors.empty());
    REQUIRE_FALSE(bad.note.empty());
}

TEST_CASE("rSquared handles a variance-free observation set", "[numeric]") {
    REQUIRE(rSquared({2.0, 2.0, 2.0}, {2.0, 2.0, 2.0}) == 0.0);
    REQUIRE(rSquared({}, {}) == 0.0);
}

TEST_CASE("rk4Integrate matches the analytic solutions", "[numeric]") {
    const double k = 0.7;
    const double y0 = 3.0;
    std::vector<double> y{y0};
    double worst = 0.0;
    rk4Integrate(0.0, 5.0, 0.001, y,
                 [&](double, const std::vector<double>& s, std::vector<double>& d) {
                     d[0] = -k * s[0];
                 },
                 [&](double t, const std::vector<double>& s) {
                     worst = std::max(worst, std::fabs(s[0] - y0 * std::exp(-k * t)));
                 });
    REQUIRE(worst < 1e-9);

    std::vector<double> z{1.0};
    double growthWorst = 0.0;
    rk4Integrate(0.0, 2.0, 0.001, z,
                 [](double, const std::vector<double>& s, std::vector<double>& d) {
                     d[0] = s[0];
                 },
                 [&](double t, const std::vector<double>& s) {
                     growthWorst = std::max(growthWorst, std::fabs(s[0] - std::exp(t)));
                 });
    REQUIRE(growthWorst < 1e-9);
}

TEST_CASE("rk4Integrate lands its final observation exactly on t1", "[numeric]") {
    std::vector<double> y{1.0};
    std::vector<double> times;
    rk4Integrate(0.0, 1.0, 0.3, y,
                 [](double, const std::vector<double>&, std::vector<double>& d) { d[0] = 0.0; },
                 [&](double t, const std::vector<double>&) { times.push_back(t); });

    REQUIRE(times.size() == 5);  // t0 plus ceil(1.0 / 0.3) = 4 steps
    REQUIRE(times.front() == 0.0);
    REQUIRE(times.back() == 1.0);  // exactly, not within a tolerance
}

TEST_CASE("trapezoid is exact on a line and empty-safe", "[numeric]") {
    const std::vector<double> x{0, 1, 2, 4, 7};  // deliberately non-uniform
    std::vector<double> y;
    for (double v : x) {
        y.push_back(2.0 * v + 1.0);
    }
    REQUIRE_THAT(trapezoid(x, y), WithinAbs(56.0, 1e-12));  // 7^2 + 7

    REQUIRE(trapezoid({1.0}, {5.0}) == 0.0);
    REQUIRE(trapezoid({}, {}) == 0.0);
}

// ---------------------------------------------------------------------------
// Phase 10 prerequisites: the rank guard, robust regression, profile intervals
// and the model-selection criterion, all in the one fitter rather than a second.
// ---------------------------------------------------------------------------

TEST_CASE("the SVD guard refuses error bars for an unidentifiable parameter pair", "[numeric]") {
    // y = (a + b) * x: the data constrain the sum and nothing else, so a Wald
    // error bar on either parameter alone would be pure fiction.
    const std::vector<double> xs{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    auto evaluate = [&](const std::vector<double>& p, std::vector<double>& r,
                        std::vector<double>& j) {
        for (std::size_t i = 0; i < xs.size(); ++i) {
            r[i] = (p[0] + p[1]) * xs[i] - 3.0 * xs[i];
            j[i * 2 + 0] = xs[i];
            j[i * 2 + 1] = xs[i];
        }
    };

    const LmResult fit = levenbergMarquardt({0.5, 0.5}, xs.size(), evaluate);
    REQUIRE_THAT(fit.params[0] + fit.params[1], WithinAbs(3.0, 1e-9));
    REQUIRE(fit.rank == 1);
    REQUIRE(fit.standardErrors.empty());
    REQUIRE(fit.covariance.empty());
    REQUIRE(fit.note.find("rank 1 of 2") != std::string::npos);
}

TEST_CASE("a full-rank fit reports its rank, condition number and covariance", "[numeric]") {
    std::vector<double> xs, ys;
    for (int i = 0; i < 20; ++i) {
        xs.push_back(0.1 * i);
        ys.push_back(2.5 * std::exp(-0.7 * xs.back()));
    }
    auto evaluate = [&](const std::vector<double>& p, std::vector<double>& r,
                        std::vector<double>& j) {
        for (std::size_t i = 0; i < xs.size(); ++i) {
            const double e = std::exp(p[1] * xs[i]);
            r[i] = p[0] * e - ys[i];
            j[i * 2 + 0] = e;
            j[i * 2 + 1] = p[0] * xs[i] * e;
        }
    };

    const LmResult fit = levenbergMarquardt({1.0, -0.1}, xs.size(), evaluate);
    REQUIRE_THAT(fit.params[0], WithinAbs(2.5, 1e-9));
    REQUIRE_THAT(fit.params[1], WithinAbs(-0.7, 1e-9));
    REQUIRE(fit.rank == 2);
    REQUIRE(std::isfinite(fit.conditionNumber));
    REQUIRE(fit.covariance.size() == 4);
}

TEST_CASE("Tukey-biweight IRLS rejects one ruined point without deleting it", "[numeric]") {
    std::vector<double> xs, ys;
    for (int i = 0; i < 21; ++i) {
        xs.push_back(i);
        ys.push_back(1.0 + 2.0 * i);
    }
    ys[10] = 500.0;  // one ruined well, sitting at the mean of x

    auto evaluate = [&](const std::vector<double>& p, std::vector<double>& r,
                        std::vector<double>& j) {
        for (std::size_t i = 0; i < xs.size(); ++i) {
            r[i] = p[0] + p[1] * xs[i] - ys[i];
            j[i * 2 + 0] = 1.0;
            j[i * 2 + 1] = xs[i];
        }
    };

    const LmResult ols = levenbergMarquardt({0.0, 1.0}, xs.size(), evaluate);
    const IrlsResult robust = tukeyBiweight({0.0, 1.0}, xs.size(), evaluate);

    // Least squares absorbs 479/21 = 22.8 of the outlier into the intercept.
    REQUIRE(std::abs(ols.params[0] - 1.0) > 20.0);
    REQUIRE_THAT(robust.fit.params[0], WithinAbs(1.0, 1e-6));
    REQUIRE_THAT(robust.fit.params[1], WithinAbs(2.0, 1e-6));
    REQUIRE(robust.converged);

    // The outlier is at weight zero, and it is still in the data set: a hollow
    // marker on the plot, not a deleted row.
    REQUIRE(robust.robustWeights.size() == xs.size());
    REQUIRE(robust.robustWeights[10] == 0.0);
    for (std::size_t i = 0; i < robust.robustWeights.size(); ++i) {
        if (i != 10) {
            REQUIRE(robust.robustWeights[i] > 0.9);
        }
    }
}

TEST_CASE("a profile interval brackets the estimate and refuses an untabulated level",
          "[numeric]") {
    // Deterministic pseudo-noise: the interval must be a real interval, so the
    // residual variance cannot be zero.
    const double noise[20] = {0.03, -0.02, 0.01, 0.04, -0.03, 0.02, -0.01, 0.00,
                              0.03, -0.04, 0.02, 0.01, -0.02, 0.03, -0.01, 0.02,
                              -0.03, 0.01, 0.00, -0.02};
    std::vector<double> xs, ys;
    for (int i = 0; i < 20; ++i) {
        xs.push_back(0.1 * i);
        ys.push_back(2.5 * std::exp(-0.7 * xs.back()) + noise[i]);
    }
    auto evaluate = [&](const std::vector<double>& p, std::vector<double>& r,
                        std::vector<double>& j) {
        for (std::size_t i = 0; i < xs.size(); ++i) {
            const double e = std::exp(p[1] * xs[i]);
            r[i] = p[0] * e - ys[i];
            j[i * 2 + 0] = e;
            j[i * 2 + 1] = p[0] * xs[i] * e;
        }
    };

    const LmResult fit = levenbergMarquardt({1.0, -0.1}, xs.size(), evaluate);
    REQUIRE(fit.standardErrors.size() == 2);

    const ProfileInterval pi = profileLikelihood(fit.params, xs.size(), evaluate, 1, 0.95);
    REQUIRE(pi.lowerFound);
    REQUIRE(pi.upperFound);
    REQUIRE(pi.lower < fit.params[1]);
    REQUIRE(pi.upper > fit.params[1]);

    // This parameter is close enough to linear that the profile and Wald widths
    // should agree to well within 15%; a large disagreement here would mean the
    // profile is not profiling.
    const double wald = 2.0 * 1.959964 * fit.standardErrors[1];
    REQUIRE(std::abs((pi.upper - pi.lower) / wald - 1.0) < 0.15);

    const ProfileInterval bad = profileLikelihood(fit.params, xs.size(), evaluate, 1, 0.925);
    REQUIRE_FALSE(bad.lowerFound);
    REQUIRE(bad.note.find("unsupported confidence level") != std::string::npos);
}

TEST_CASE("aicc matches its formula and is infinite where the correction is undefined",
          "[numeric]") {
    const double expected = 10.0 * std::log(0.1) + 4.0 + 2.0 * 2.0 * 3.0 / 7.0;
    REQUIRE_THAT(aicc(1.0, 10, 2), WithinAbs(expected, 1e-12));
    REQUIRE(std::isinf(aicc(1.0, 3, 2)));   // n <= k + 1
    REQUIRE(std::isinf(aicc(0.0, 10, 2)));  // log(0) is not a model comparison
}
