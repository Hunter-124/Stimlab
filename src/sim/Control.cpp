#include "sim/Control.h"

#include <algorithm>
#include <cmath>

#include <Eigen/Dense>

#include "sim/Solvers.h"

namespace biocad::sim {
namespace {

// The parameters that scale a reaction's RATE without changing its equilibrium
// constant: both directions of a reversible law move together, which is what an
// enzyme-amount perturbation does. Scaling only kf would perturb thermodynamics and
// silently violate Wegscheider.
std::vector<std::size_t> rateScalingParameters(RateLaw law) {
    switch (law) {
        case RateLaw::MassAction: return {0};
        case RateLaw::ReversibleMassAction: return {0, 1};
        case RateLaw::MichaelisMenten: return {0};
        case RateLaw::ReversibleMichaelisMenten: return {0, 2};
        case RateLaw::Hill: return {0};
    }
    return {0};
}

NetworkSpec scaled(const NetworkSpec& spec, std::size_t reaction, double factor) {
    NetworkSpec out = spec;
    if (reaction >= out.reactions.size()) return out;
    for (std::size_t i : rateScalingParameters(out.reactions[reaction].law))
        if (i < out.reactions[reaction].parameters.size())
            out.reactions[reaction].parameters[i] *= factor;
    return out;
}

}  // namespace

SteadyState steadyState(const NetworkSpec& spec, double horizon) {
    SteadyState out;
    Network net;
    std::string error;
    if (!Network::compile(spec, net, &error)) {
        out.note = error;
        return out;
    }
    IntegrationOptions o;
    o.horizon = horizon;
    o.outputPoints = 2;
    o.relativeTolerance = 1e-10;
    o.absoluteTolerance = 1e-14;
    o.initialStep = 1e-8;
    const TimeCourse tc = integrate(spec, o);
    if (tc.times.empty()) {
        out.note = tc.solver.note;
        return out;
    }
    const std::size_t ns = net.speciesCount();
    std::vector<double> c(ns);
    for (std::size_t i = 0; i < ns; ++i) c[i] = tc.trajectories[i].back();

    // Newton polish. The Jacobian is singular by exactly the number of conservation
    // laws, so the step is taken as a least-squares solution: it moves within the
    // reachable subspace and leaves the conserved totals alone.
    std::vector<double> f, j;
    for (int iteration = 0; iteration < 60; ++iteration) {
        net.derivatives(c, f);
        double residual = 0;
        for (double v : f) residual = std::max(residual, std::abs(v));
        out.residual = residual;
        if (residual < 1e-13) { out.converged = true; break; }
        net.jacobian(c, j);
        Eigen::MatrixXd jm(static_cast<Eigen::Index>(ns), static_cast<Eigen::Index>(ns));
        Eigen::VectorXd fv(static_cast<Eigen::Index>(ns));
        for (std::size_t r = 0; r < ns; ++r) {
            fv[static_cast<Eigen::Index>(r)] = -f[r];
            for (std::size_t k = 0; k < ns; ++k)
                jm(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(k)) = j[r * ns + k];
        }
        const Eigen::VectorXd step =
            jm.completeOrthogonalDecomposition().solve(fv);
        double norm = 0;
        for (std::size_t i = 0; i < ns; ++i)
            norm = std::max(norm, std::abs(step[static_cast<Eigen::Index>(i)]));
        // A damped step: a full Newton step from a partially converged transient can
        // leave the positive cone entirely.
        const double damping = norm > 0.25 ? 0.25 / norm : 1.0;
        for (std::size_t i = 0; i < ns; ++i)
            if (!net.boundary(i)) c[i] += damping * step[static_cast<Eigen::Index>(i)];
    }
    out.concentrations = c;
    net.rates(c, out.fluxes);
    if (!out.converged)
        out.note = "steady state not reached: worst |dc/dt| is " + std::to_string(out.residual) +
                   ". A control analysis from a non-steady state is meaningless, so the "
                   "residual is reported rather than hidden.";
    return out;
}

std::vector<SweepPoint> sweep(const NetworkSpec& spec, std::size_t reaction,
                              const std::vector<double>& factors, double horizon) {
    std::vector<SweepPoint> out;
    for (double f : factors) out.push_back({f, steadyState(scaled(spec, reaction, f), horizon)});
    return out;
}

ControlAnalysis controlAnalysis(const NetworkSpec& spec, double horizon, double perturbation) {
    ControlAnalysis out;
    Network net;
    std::string error;
    if (!Network::compile(spec, net, &error)) {
        out.warnings.push_back("network does not compile: " + error);
        return out;
    }
    out.reactionIds = net.reactionIds();
    out.speciesIds = net.speciesIds();
    const std::size_t ns = net.speciesCount(), nr = net.reactionCount();

    const SteadyState base = steadyState(spec, horizon);
    if (!base.converged) {
        out.warnings.push_back(base.note.empty() ? "no steady state was reached" : base.note);
        return out;
    }

    // Elasticities eps_{k,j} = (c_j / v_k) * dv_k/dc_j at the reference state.
    std::vector<double> dvdc;
    net.rateJacobian(base.concentrations, dvdc);
    out.elasticities.assign(nr, std::vector<double>(ns, 0.0));
    for (std::size_t k = 0; k < nr; ++k)
        for (std::size_t j = 0; j < ns; ++j)
            out.elasticities[k][j] = std::abs(base.fluxes[k]) > 1e-300
                                         ? base.concentrations[j] * dvdc[k * ns + j] /
                                               base.fluxes[k]
                                         : 0.0;

    out.fluxControlCoefficients.assign(nr, std::vector<double>(nr, 0.0));
    out.concentrationControlCoefficients.assign(ns, std::vector<double>(nr, 0.0));
    for (std::size_t k = 0; k < nr; ++k) {
        const SteadyState up = steadyState(scaled(spec, k, 1.0 + perturbation), horizon);
        const SteadyState down = steadyState(scaled(spec, k, 1.0 - perturbation), horizon);
        if (!up.converged || !down.converged) {
            out.warnings.push_back("the perturbation of reaction '" + out.reactionIds[k] +
                                   "' did not reach a steady state; its control coefficients are "
                                   "left at zero");
            continue;
        }
        const double dlnp = std::log(1.0 + perturbation) - std::log(1.0 - perturbation);
        for (std::size_t r = 0; r < nr; ++r) {
            if (std::abs(base.fluxes[r]) < 1e-12) continue;
            if (up.fluxes[r] <= 0 || down.fluxes[r] <= 0) {
                // A flux that changes sign under the perturbation has no log
                // derivative; the linear one is used and the fact is recorded.
                out.fluxControlCoefficients[r][k] =
                    (up.fluxes[r] - down.fluxes[r]) / (2.0 * perturbation) / base.fluxes[r];
                continue;
            }
            out.fluxControlCoefficients[r][k] =
                (std::log(up.fluxes[r]) - std::log(down.fluxes[r])) / dlnp;
        }
        for (std::size_t i = 0; i < ns; ++i) {
            if (net.boundary(i) || base.concentrations[i] <= 1e-12) continue;
            if (up.concentrations[i] <= 0 || down.concentrations[i] <= 0) continue;
            out.concentrationControlCoefficients[i][k] =
                (std::log(up.concentrations[i]) - std::log(down.concentrations[i])) / dlnp;
        }
    }

    // The theorems. Flux control coefficients sum to 1 for every flux that is
    // nonzero; concentration control coefficients sum to 0.
    double summation = 0;
    for (std::size_t r = 0; r < nr; ++r) {
        if (std::abs(base.fluxes[r]) < 1e-12) continue;
        double s = 0;
        for (std::size_t k = 0; k < nr; ++k) s += out.fluxControlCoefficients[r][k];
        summation = std::max(summation, std::abs(s - 1.0));
    }
    for (std::size_t i = 0; i < ns; ++i) {
        if (net.boundary(i) || base.concentrations[i] <= 1e-12) continue;
        double s = 0;
        for (std::size_t k = 0; k < nr; ++k) s += out.concentrationControlCoefficients[i][k];
        summation = std::max(summation, std::abs(s));
    }
    out.summationResidual = summation;

    double connectivity = 0;
    for (std::size_t r = 0; r < nr; ++r) {
        if (std::abs(base.fluxes[r]) < 1e-12) continue;
        for (std::size_t j = 0; j < ns; ++j) {
            if (net.boundary(j) || base.concentrations[j] <= 1e-12) continue;
            double s = 0;
            for (std::size_t k = 0; k < nr; ++k)
                s += out.fluxControlCoefficients[r][k] * out.elasticities[k][j];
            connectivity = std::max(connectivity, std::abs(s));
        }
    }
    out.connectivityResidual = connectivity;
    if (summation > 1e-4 || connectivity > 1e-4)
        out.warnings.push_back(
            "the summation residual is " + std::to_string(summation) +
            " and the connectivity residual " + std::to_string(connectivity) +
            "; both are theorems, so a residual this large means the coefficients below are "
            "not trustworthy");
    return out;
}

}  // namespace biocad::sim
