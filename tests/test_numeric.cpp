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
