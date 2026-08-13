// numeric/Ode.h - fixed-step Runge-Kutta 4 integration.
//
// The PK models are integrated rather than evaluated in closed form because the
// closed forms are numerically singular when ka == ke and cannot express nonlinear
// (Michaelis-Menten) elimination. The closed forms still exist - as the unit-test
// oracle and the summary readout - and the integrator is asserted against them.
//
// Fixed step, not adaptive: the systems here are small, smooth and short-horizon,
// and a deterministic step count makes a simulation byte-reproducible.
#pragma once

#include <cstddef>
#include <functional>
#include <vector>

namespace biocad::numeric {

// dy/dt = f(t, y). `dydt` is pre-sized to y.size() and must be fully written.
using OdeDerivative = std::function<void(double t, const std::vector<double>& y,
                                         std::vector<double>& dydt)>;

// Advance `y` from t0 to t0 + h with one classical RK4 step.
void rk4Step(double t0, double h, std::vector<double>& y, const OdeDerivative& f);

// Integrate from t0 to t1 with step h, invoking `observe(t, y)` at t0 and after each
// step. The final step is shortened so the last observation lands exactly on t1.
void rk4Integrate(double t0, double t1, double h, std::vector<double>& y,
                  const OdeDerivative& f,
                  const std::function<void(double, const std::vector<double>&)>& observe);

// Linear-up trapezoidal area under a sampled curve. Used for AUC over a simulated
// horizon, where the sampling is dense and uniform.
double trapezoid(const std::vector<double>& x, const std::vector<double>& y);

}  // namespace biocad::numeric
