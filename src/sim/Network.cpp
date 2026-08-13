#include "sim/Network.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <unordered_map>

#include <Eigen/Dense>

namespace biocad::sim {
namespace {

constexpr double kNullTol = 1e-9;

const char* const kMassAction[]        = {"k", nullptr};
const char* const kReversibleMass[]    = {"kf", "kr", nullptr};
const char* const kMichaelisMenten[]   = {"Vmax", "Km", nullptr};
const char* const kReversibleMm[]      = {"Vf", "Kms", "Vr", "Kmp", nullptr};
const char* const kHill[]              = {"Vmax", "K", "n", nullptr};

// c^a for a non-negative concentration, written so that 0^0 == 1 and 0^a == 0
// without calling std::pow on the degenerate pair.
double powc(double c, double a) {
    if (a == 1.0) return c;
    if (a == 0.0) return 1.0;
    if (c == 0.0) return a > 0.0 ? 0.0 : std::numeric_limits<double>::infinity();
    return std::pow(c, a);
}

std::string number(double v) {
    if (std::abs(v - std::round(v)) < 1e-9) {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%lld", static_cast<long long>(std::llround(v)));
        return buf;
    }
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.4g", v);
    return buf;
}

// Scale a null-space vector into a readable, sign-canonical form: tiny entries are
// zeroed, the vector is divided by its smallest nonzero magnitude (which turns the
// usual stoichiometric kernels into small integers) and the first nonzero entry is
// made positive.
void canonicalize(std::vector<double>& v) {
    double biggest = 0;
    for (double x : v) biggest = std::max(biggest, std::abs(x));
    if (biggest <= 0) return;
    for (double& x : v)
        if (std::abs(x) < kNullTol * biggest) x = 0;
    double smallest = biggest;
    for (double x : v)
        if (x != 0) smallest = std::min(smallest, std::abs(x));
    for (double& x : v) x /= smallest;
    for (double x : v) {
        if (x != 0) {
            if (x < 0)
                for (double& y : v) y = -y;
            break;
        }
    }
    for (double& x : v)
        if (std::abs(x - std::round(x)) < 1e-9) x = std::round(x);
}

}  // namespace

std::size_t rateLawParameterCount(RateLaw law) {
    switch (law) {
        case RateLaw::MassAction: return 1;
        case RateLaw::ReversibleMassAction: return 2;
        case RateLaw::MichaelisMenten: return 2;
        case RateLaw::ReversibleMichaelisMenten: return 4;
        case RateLaw::Hill: return 3;
    }
    return 0;
}

const char* const* rateLawParameterNames(RateLaw law) {
    switch (law) {
        case RateLaw::MassAction: return kMassAction;
        case RateLaw::ReversibleMassAction: return kReversibleMass;
        case RateLaw::MichaelisMenten: return kMichaelisMenten;
        case RateLaw::ReversibleMichaelisMenten: return kReversibleMm;
        case RateLaw::Hill: return kHill;
    }
    return kMassAction;
}

bool Network::compile(const NetworkSpec& spec, Network& out, std::string* error) {
    out = Network{};
    std::unordered_map<std::string, std::size_t> index;
    for (const auto& s : spec.species) {
        if (index.count(s.id)) {
            if (error) *error = "duplicate species id '" + s.id + "'";
            return false;
        }
        index.emplace(s.id, out.species_.size());
        out.species_.push_back({s.boundary, s.initialConcentration});
        out.speciesIds_.push_back(s.id);
    }
    const std::size_t ns = out.species_.size();
    if (ns == 0) {
        if (error) *error = "network has no species";
        return false;
    }

    out.s_.assign(ns * spec.reactions.size(), 0.0);
    for (std::size_t r = 0; r < spec.reactions.size(); ++r) {
        const ReactionSpec& rs = spec.reactions[r];
        Reaction rx;
        rx.law = rs.law;
        rx.parameters = rs.parameters;
        const std::size_t want = rateLawParameterCount(rs.law);
        if (rx.parameters.size() != want) {
            if (error)
                *error = "reaction '" + rs.id + "' needs " + std::to_string(want) +
                         " parameter(s) for its rate law but was given " +
                         std::to_string(rx.parameters.size());
            return false;
        }
        auto resolve = [&](const std::pair<std::string, double>& p, std::vector<Term>& into,
                           double sign) -> bool {
            auto it = index.find(p.first);
            if (it == index.end()) {
                if (error)
                    *error = "reaction '" + rs.id + "' references unknown species '" + p.first + "'";
                return false;
            }
            into.push_back({it->second, p.second});
            out.s_[it->second * spec.reactions.size() + r] += sign * p.second;
            return true;
        };
        for (const auto& p : rs.reactants)
            if (!resolve(p, rx.reactants, -1.0)) return false;
        for (const auto& p : rs.products)
            if (!resolve(p, rx.products, +1.0)) return false;

        const bool needsSubstrate = rs.law == RateLaw::MichaelisMenten ||
                                    rs.law == RateLaw::ReversibleMichaelisMenten ||
                                    rs.law == RateLaw::Hill;
        if (needsSubstrate && rx.reactants.empty()) {
            if (error) *error = "reaction '" + rs.id + "' has a saturable rate law but no substrate";
            return false;
        }
        if (rs.law == RateLaw::ReversibleMichaelisMenten && rx.products.empty()) {
            if (error) *error = "reaction '" + rs.id + "' is a reversible Michaelis-Menten step "
                                "but has no product";
            return false;
        }
        out.reactions_.push_back(std::move(rx));
        out.reactionIds_.push_back(rs.id);
    }
    return true;
}

std::vector<double> Network::initialState() const {
    std::vector<double> c(species_.size());
    for (std::size_t i = 0; i < species_.size(); ++i) c[i] = species_[i].initial;
    return c;
}

void Network::rates(const std::vector<double>& c, std::vector<double>& v) const {
    v.assign(reactions_.size(), 0.0);
    for (std::size_t r = 0; r < reactions_.size(); ++r) {
        const Reaction& rx = reactions_[r];
        const std::vector<double>& p = rx.parameters;
        switch (rx.law) {
            case RateLaw::MassAction: {
                double prod = p[0];
                for (const Term& t : rx.reactants) prod *= powc(c[t.species], t.stoichiometry);
                v[r] = prod;
                break;
            }
            case RateLaw::ReversibleMassAction: {
                double f = p[0], b = p[1];
                for (const Term& t : rx.reactants) f *= powc(c[t.species], t.stoichiometry);
                for (const Term& t : rx.products) b *= powc(c[t.species], t.stoichiometry);
                v[r] = f - b;
                break;
            }
            case RateLaw::MichaelisMenten: {
                const double s = c[rx.reactants.front().species];
                v[r] = p[0] * s / (p[1] + s);
                break;
            }
            case RateLaw::ReversibleMichaelisMenten: {
                const double s = c[rx.reactants.front().species];
                const double q = c[rx.products.front().species];
                const double num = p[0] * s / p[1] - p[2] * q / p[3];
                const double den = 1.0 + s / p[1] + q / p[3];
                v[r] = num / den;
                break;
            }
            case RateLaw::Hill: {
                const double s = c[rx.reactants.front().species];
                const double kn = powc(p[1], p[2]);
                const double sn = powc(s, p[2]);
                v[r] = p[0] * sn / (kn + sn);
                break;
            }
        }
    }
}

void Network::derivatives(const std::vector<double>& c, std::vector<double>& dcdt) const {
    std::vector<double> v;
    rates(c, v);
    const std::size_t nr = reactions_.size();
    dcdt.assign(species_.size(), 0.0);
    for (std::size_t i = 0; i < species_.size(); ++i) {
        if (species_[i].boundary) continue;  // a clamped pool is constant by definition
        double sum = 0;
        for (std::size_t r = 0; r < nr; ++r) {
            const double sij = s_[i * nr + r];
            if (sij != 0) sum += sij * v[r];
        }
        dcdt[i] = sum;
    }
}

void Network::rateJacobian(const std::vector<double>& c, std::vector<double>& dvdc) const {
    const std::size_t ns = species_.size();
    dvdc.assign(reactions_.size() * ns, 0.0);
    // The EXCLUDING product: for term t on species j with stoichiometry a, the
    // derivative is a * c_j^(a-1) * prod over the OTHER terms. Never v * a / c_j -
    // that form is undefined at c_j = 0, where the true derivative is finite.
    auto accumulateMassAction = [&](std::size_t r, const std::vector<Term>& terms, double k,
                                    double sign) {
        for (std::size_t t = 0; t < terms.size(); ++t) {
            double other = k;
            for (std::size_t u = 0; u < terms.size(); ++u)
                if (u != t) other *= powc(c[terms[u].species], terms[u].stoichiometry);
            const double a = terms[t].stoichiometry;
            const double self = a * powc(c[terms[t].species], a - 1.0);
            dvdc[r * ns + terms[t].species] += sign * other * self;
        }
    };
    for (std::size_t r = 0; r < reactions_.size(); ++r) {
        const Reaction& rx = reactions_[r];
        const std::vector<double>& p = rx.parameters;
        switch (rx.law) {
            case RateLaw::MassAction:
                accumulateMassAction(r, rx.reactants, p[0], +1.0);
                break;
            case RateLaw::ReversibleMassAction:
                accumulateMassAction(r, rx.reactants, p[0], +1.0);
                accumulateMassAction(r, rx.products, p[1], -1.0);
                break;
            case RateLaw::MichaelisMenten: {
                const std::size_t j = rx.reactants.front().species;
                const double s = c[j];
                const double d = p[1] + s;
                dvdc[r * ns + j] += p[0] * p[1] / (d * d);
                break;
            }
            case RateLaw::ReversibleMichaelisMenten: {
                const std::size_t js = rx.reactants.front().species;
                const std::size_t jp = rx.products.front().species;
                const double s = c[js], q = c[jp];
                const double num = p[0] * s / p[1] - p[2] * q / p[3];
                const double den = 1.0 + s / p[1] + q / p[3];
                dvdc[r * ns + js] += (p[0] / p[1] * den - num / p[1]) / (den * den);
                dvdc[r * ns + jp] += (-p[2] / p[3] * den - num / p[3]) / (den * den);
                break;
            }
            case RateLaw::Hill: {
                const std::size_t j = rx.reactants.front().species;
                const double s = c[j], n = p[2];
                const double kn = powc(p[1], n);
                const double sn = powc(s, n);
                const double den = kn + sn;
                // n * s^(n-1) is the excluding form again; at s = 0 with n > 1 the
                // derivative is genuinely 0, and with n == 1 it is Vmax/K.
                const double dsn = n * powc(s, n - 1.0);
                dvdc[r * ns + j] += p[0] * kn * dsn / (den * den);
                break;
            }
        }
    }
}

void Network::jacobian(const std::vector<double>& c, std::vector<double>& j) const {
    const std::size_t ns = species_.size(), nr = reactions_.size();
    std::vector<double> dvdc;
    rateJacobian(c, dvdc);
    j.assign(ns * ns, 0.0);
    for (std::size_t i = 0; i < ns; ++i) {
        if (species_[i].boundary) continue;
        for (std::size_t r = 0; r < nr; ++r) {
            const double sij = s_[i * nr + r];
            if (sij == 0) continue;
            for (std::size_t k = 0; k < ns; ++k) {
                const double d = dvdc[r * ns + k];
                if (d != 0) j[i * ns + k] += sij * d;
            }
        }
    }
}

NetworkSpec analyze(const NetworkSpec& spec) {
    NetworkSpec out = spec;
    out.conservationLaws.clear();
    out.conservationLabels.clear();
    out.thermodynamicCycles.clear();
    out.wegscheiderViolations.clear();

    Network net;
    std::string error;
    if (!Network::compile(spec, net, &error)) {
        out.warnings.push_back("network does not compile: " + error);
        return out;
    }
    const std::size_t ns = net.speciesCount(), nr = net.reactionCount();
    if (nr == 0) {
        out.warnings.push_back("network has no reactions; nothing to analyze");
        return out;
    }

    // Boundary species are held constant by the CLAMP, not by stoichiometry, so
    // they are excluded from the conservation analysis: a moiety that flows into a
    // clamped pool is not conserved, and reporting it as such would be a lie the
    // integrator immediately contradicts.
    std::vector<std::size_t> variable;
    for (std::size_t i = 0; i < ns; ++i)
        if (!net.boundary(i)) variable.push_back(i);

    const std::vector<double>& s = net.stoichiometry();
    Eigen::MatrixXd sv(static_cast<Eigen::Index>(variable.size()), static_cast<Eigen::Index>(nr));
    for (std::size_t a = 0; a < variable.size(); ++a)
        for (std::size_t r = 0; r < nr; ++r)
            sv(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(r)) =
                s[variable[a] * nr + r];

    // Left null space: y^T . S = 0  <=>  y in ker(S^T). Each basis vector is a
    // conserved moiety, and its value over time is the integration's own audit.
    {
        Eigen::MatrixXd st = sv.transpose();
        Eigen::FullPivLU<Eigen::MatrixXd> lu(st);
        lu.setThreshold(1e-10);
        if (lu.dimensionOfKernel() > 0) {
            Eigen::MatrixXd k = lu.kernel();
            for (Eigen::Index col = 0; col < k.cols(); ++col) {
                std::vector<double> law(ns, 0.0);
                for (std::size_t a = 0; a < variable.size(); ++a)
                    law[variable[a]] = k(static_cast<Eigen::Index>(a), col);
                canonicalize(law);
                std::string label;
                for (std::size_t i = 0; i < ns; ++i) {
                    if (law[i] == 0) continue;
                    if (!label.empty()) label += law[i] < 0 ? " - " : " + ";
                    else if (law[i] < 0) label += "-";
                    const double m = std::abs(law[i]);
                    if (m != 1.0) label += number(m) + " ";
                    label += net.speciesIds()[i];
                }
                out.conservationLaws.push_back(std::move(law));
                out.conservationLabels.push_back(label.empty() ? "(degenerate)" : label);
            }
        }
    }

    // Right null space: S . k = 0. A basis vector supported on reversible
    // mass-action steps is a thermodynamic cycle, and Wegscheider requires
    // prod (kf/kr)^k_i = 1 around it.
    {
        Eigen::FullPivLU<Eigen::MatrixXd> lu(sv);
        lu.setThreshold(1e-10);
        if (lu.dimensionOfKernel() > 0) {
            Eigen::MatrixXd k = lu.kernel();
            for (Eigen::Index col = 0; col < k.cols(); ++col) {
                std::vector<double> cycle(nr, 0.0);
                for (std::size_t r = 0; r < nr; ++r)
                    cycle[r] = k(static_cast<Eigen::Index>(r), col);
                canonicalize(cycle);
                bool checkable = true;
                double logProduct = 0;
                std::string uncheckable;
                for (std::size_t r = 0; r < nr; ++r) {
                    if (cycle[r] == 0) continue;
                    const ReactionSpec& rs = spec.reactions[r];
                    if (rs.law != RateLaw::ReversibleMassAction) {
                        checkable = false;
                        uncheckable = rs.id;
                        break;
                    }
                    const double kf = rs.parameters[0], kr2 = rs.parameters[1];
                    if (!(kf > 0) || !(kr2 > 0)) {
                        checkable = false;
                        uncheckable = rs.id;
                        break;
                    }
                    logProduct += cycle[r] * (std::log(kf) - std::log(kr2));
                }
                if (checkable) {
                    if (std::abs(logProduct) > 1e-8) {
                        out.wegscheiderViolations.push_back(
                            "cycle " + std::to_string(out.thermodynamicCycles.size() + 1) +
                            ": product of equilibrium constants is " +
                            number(std::exp(logProduct)) +
                            ", not 1 - the rate constants describe a cycle that does net work "
                            "at equilibrium and no choice of solver can fix that");
                    }
                } else if (!uncheckable.empty()) {
                    out.warnings.push_back(
                        "cycle " + std::to_string(out.thermodynamicCycles.size() + 1) +
                        " involves '" + uncheckable +
                        "', which is not a reversible mass-action step with two positive rate "
                        "constants, so its Wegscheider condition is not checkable");
                }
                out.thermodynamicCycles.push_back(std::move(cycle));
            }
        }
    }
    return out;
}

}  // namespace biocad::sim
