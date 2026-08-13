#include "numeric/Ode.h"

#include <cmath>

namespace biocad::numeric {

void rk4Step(double t0, double h, std::vector<double>& y, const OdeDerivative& f) {
    const std::size_t n = y.size();
    if (n == 0 || h == 0.0) {
        return;
    }

    std::vector<double> k1(n, 0.0), k2(n, 0.0), k3(n, 0.0), k4(n, 0.0), tmp(n, 0.0);

    f(t0, y, k1);
    for (std::size_t i = 0; i < n; ++i) {
        tmp[i] = y[i] + 0.5 * h * k1[i];
    }
    f(t0 + 0.5 * h, tmp, k2);
    for (std::size_t i = 0; i < n; ++i) {
        tmp[i] = y[i] + 0.5 * h * k2[i];
    }
    f(t0 + 0.5 * h, tmp, k3);
    for (std::size_t i = 0; i < n; ++i) {
        tmp[i] = y[i] + h * k3[i];
    }
    f(t0 + h, tmp, k4);

    for (std::size_t i = 0; i < n; ++i) {
        y[i] += (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}

void rk4Integrate(double t0, double t1, double h, std::vector<double>& y, const OdeDerivative& f,
                  const std::function<void(double, const std::vector<double>&)>& observe) {
    if (observe) {
        observe(t0, y);
    }
    const double span = t1 - t0;
    if (!(span > 0.0) || !(h > 0.0) || y.empty()) {
        return;
    }

    // The step count is derived once and the time of each observation is computed
    // from it, so no float error accumulates and the final sample lands exactly on
    // t1 - simulations from the same spec must be byte-reproducible.
    const long long steps = static_cast<long long>(std::ceil(span / h - 1e-12));
    for (long long i = 0; i < steps; ++i) {
        const double tStart = (i == 0) ? t0 : t0 + static_cast<double>(i) * h;
        const double tEnd = (i == steps - 1) ? t1 : t0 + static_cast<double>(i + 1) * h;
        rk4Step(tStart, tEnd - tStart, y, f);
        if (observe) {
            observe(tEnd, y);
        }
    }
}

double trapezoid(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() < 2 || x.size() != y.size()) {
        return 0.0;
    }
    double area = 0.0;
    for (std::size_t i = 1; i < x.size(); ++i) {
        area += 0.5 * (x[i] - x[i - 1]) * (y[i] + y[i - 1]);
    }
    return area;
}

}  // namespace biocad::numeric
