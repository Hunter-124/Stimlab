#include "sim/Flux.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <numeric>
#include <unordered_map>

#include "chem/Formula.h"
#include "sim/Network.h"

namespace biocad::sim {
namespace {

constexpr double kTol = 1e-9;

std::string fmt(double v) {
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%.6g", v);
    return buf;
}

// Standard-form tableau simplex, min c'x s.t. Ax = b, x >= 0.
// Bland's rule (lowest index among eligible) is used for BOTH the entering and the
// leaving choice: FBA problems are massively degenerate, and Dantzig's rule can cycle
// on them forever. Bland's rule is slower per iteration count and provably finite,
// which is the correct trade for a solver whose answers are published as fluxes.
struct Standard {
    std::size_t                       n = 0, m = 0;
    std::vector<double>               a;   // m x n
    std::vector<double>               b;
    std::vector<double>               c;
};

struct SimplexOut {
    bool                feasible = false;
    bool                bounded = true;
    std::vector<double> x;
    int                 iterations = 0;
};

// Pivot the tableau so that column `col` becomes basic in row `row`.
void pivot(std::vector<double>& t, std::size_t rows, std::size_t cols, std::size_t row,
           std::size_t col) {
    const double p = t[row * cols + col];
    for (std::size_t j = 0; j < cols; ++j) t[row * cols + j] /= p;
    for (std::size_t i = 0; i < rows; ++i) {
        if (i == row) continue;
        const double f = t[i * cols + col];
        if (f == 0) continue;
        for (std::size_t j = 0; j < cols; ++j) t[i * cols + j] -= f * t[row * cols + j];
    }
}

// One simplex phase over a tableau with `m` constraint rows plus one objective row.
// `basis` is updated in place. Returns false when the problem is unbounded.
bool simplexPhase(std::vector<double>& t, std::size_t m, std::size_t cols,
                  std::vector<std::size_t>& basis, int& iterations, std::size_t columnLimit) {
    const std::size_t objRow = m;
    for (int guard = 0; guard < 200000; ++guard) {
        std::size_t enter = columnLimit;
        for (std::size_t j = 0; j < columnLimit; ++j) {
            if (t[objRow * cols + j] < -kTol) { enter = j; break; }   // Bland: lowest index
        }
        if (enter == columnLimit) return true;   // optimal
        std::size_t leave = m;
        double best = 0;
        for (std::size_t i = 0; i < m; ++i) {
            const double aij = t[i * cols + enter];
            if (aij <= kTol) continue;
            const double ratio = t[i * cols + cols - 1] / aij;
            if (leave == m || ratio < best - 1e-12 ||
                (std::abs(ratio - best) <= 1e-12 && basis[i] < basis[leave])) {
                best = ratio;
                leave = i;
            }
        }
        if (leave == m) return false;   // unbounded
        pivot(t, m + 1, cols, leave, enter);
        basis[leave] = enter;
        ++iterations;
    }
    return true;
}

SimplexOut simplexStandard(const Standard& s) {
    SimplexOut out;
    const std::size_t n = s.n, m = s.m;
    // Columns: n structural + m artificial + 1 right-hand side.
    const std::size_t cols = n + m + 1;
    std::vector<double> t((m + 1) * cols, 0.0);
    std::vector<double> b = s.b;
    std::vector<double> a = s.a;
    for (std::size_t i = 0; i < m; ++i) {
        if (b[i] < 0) {   // keep the right-hand side non-negative
            b[i] = -b[i];
            for (std::size_t j = 0; j < n; ++j) a[i * n + j] = -a[i * n + j];
        }
        for (std::size_t j = 0; j < n; ++j) t[i * cols + j] = a[i * n + j];
        t[i * cols + n + i] = 1.0;
        t[i * cols + cols - 1] = b[i];
    }
    std::vector<std::size_t> basis(m);
    for (std::size_t i = 0; i < m; ++i) basis[i] = n + i;

    // Phase 1: minimise the sum of artificials. The objective row is the negated sum
    // of the constraint rows over the structural columns.
    for (std::size_t j = 0; j < n; ++j) {
        double sum = 0;
        for (std::size_t i = 0; i < m; ++i) sum += t[i * cols + j];
        t[m * cols + j] = -sum;
    }
    {
        double sum = 0;
        for (std::size_t i = 0; i < m; ++i) sum += t[i * cols + cols - 1];
        t[m * cols + cols - 1] = -sum;
    }
    if (!simplexPhase(t, m, cols, basis, out.iterations, n + m)) return out;
    if (-t[m * cols + cols - 1] > 1e-7) return out;   // phase-1 optimum > 0: infeasible
    out.feasible = true;

    // Drive any artificial still in the basis out; a row that cannot be pivoted is
    // linearly dependent and is simply left alone (its artificial sits at zero).
    for (std::size_t i = 0; i < m; ++i) {
        if (basis[i] < n) continue;
        for (std::size_t j = 0; j < n; ++j) {
            if (std::abs(t[i * cols + j]) > kTol) {
                pivot(t, m + 1, cols, i, j);
                basis[i] = j;
                break;
            }
        }
    }

    // Phase 2: the real objective, with artificial columns forbidden.
    for (std::size_t j = 0; j < cols; ++j) t[m * cols + j] = 0.0;
    for (std::size_t j = 0; j < n; ++j) t[m * cols + j] = s.c[j];
    for (std::size_t i = 0; i < m; ++i) {
        const std::size_t bj = basis[i];
        if (bj >= n) continue;
        const double f = t[m * cols + bj];
        if (f == 0) continue;
        for (std::size_t j = 0; j < cols; ++j) t[m * cols + j] -= f * t[i * cols + j];
    }
    if (!simplexPhase(t, m, cols, basis, out.iterations, n)) {
        out.bounded = false;
        return out;
    }
    out.x.assign(n, 0.0);
    for (std::size_t i = 0; i < m; ++i)
        if (basis[i] < n) out.x[basis[i]] = t[i * cols + cols - 1];
    return out;
}

}  // namespace

LpResult solveLp(const LpProblem& problem) {
    LpResult out;
    const std::size_t n = problem.variables;
    // Shift every variable to its lower bound, x = l + y with y >= 0, and add one
    // slack row per finite upper bound: y_j + s_j = u_j - l_j. A free variable
    // (lower <= -kInfiniteBound) is split into a difference of two non-negative
    // variables, which is the standard exact transformation.
    std::vector<std::size_t> plusIndex(n), minusIndex(n, static_cast<std::size_t>(-1));
    std::vector<double>      shift(n, 0.0);
    std::size_t              columns = 0;
    for (std::size_t j = 0; j < n; ++j) {
        if (problem.lower[j] <= -kInfiniteBound) {
            plusIndex[j] = columns++;
            minusIndex[j] = columns++;
            shift[j] = 0.0;
        } else {
            plusIndex[j] = columns++;
            shift[j] = problem.lower[j];
        }
    }
    std::vector<std::size_t> upperRows;
    for (std::size_t j = 0; j < n; ++j)
        if (problem.upper[j] < kInfiniteBound) upperRows.push_back(j);
    const std::size_t slackStart = columns;
    columns += upperRows.size();

    Standard s;
    s.n = columns;
    s.m = problem.constraints + upperRows.size();
    s.a.assign(s.m * s.n, 0.0);
    s.b.assign(s.m, 0.0);
    s.c.assign(s.n, 0.0);
    for (std::size_t i = 0; i < problem.constraints; ++i) {
        double rhs = problem.b[i];
        for (std::size_t j = 0; j < n; ++j) {
            const double aij = problem.a[i * n + j];
            if (aij == 0) continue;
            s.a[i * s.n + plusIndex[j]] += aij;
            if (minusIndex[j] != static_cast<std::size_t>(-1))
                s.a[i * s.n + minusIndex[j]] -= aij;
            rhs -= aij * shift[j];
        }
        s.b[i] = rhs;
    }
    for (std::size_t k = 0; k < upperRows.size(); ++k) {
        const std::size_t row = problem.constraints + k;
        const std::size_t j = upperRows[k];
        s.a[row * s.n + plusIndex[j]] = 1.0;
        if (minusIndex[j] != static_cast<std::size_t>(-1))
            s.a[row * s.n + minusIndex[j]] = -1.0;
        s.a[row * s.n + slackStart + k] = 1.0;
        s.b[row] = problem.upper[j] - shift[j];
    }
    double constant = 0;
    for (std::size_t j = 0; j < n; ++j) {
        s.c[plusIndex[j]] += problem.c[j];
        if (minusIndex[j] != static_cast<std::size_t>(-1)) s.c[minusIndex[j]] -= problem.c[j];
        constant += problem.c[j] * shift[j];
    }

    const SimplexOut sol = simplexStandard(s);
    out.iterations = sol.iterations;
    out.feasible = sol.feasible;
    out.bounded = sol.bounded;
    if (!sol.feasible) {
        out.status = "infeasible";
        return out;
    }
    if (!sol.bounded) {
        out.status = "unbounded";
        return out;
    }
    out.x.assign(n, 0.0);
    for (std::size_t j = 0; j < n; ++j) {
        double v = shift[j] + sol.x[plusIndex[j]];
        if (minusIndex[j] != static_cast<std::size_t>(-1)) v -= sol.x[minusIndex[j]];
        out.x[j] = v;
    }
    // The objective is recomputed from the recovered x rather than read off the
    // tableau: the shift and split transformations both move the constant term, and
    // one honest multiplication is cheaper than keeping three of them in agreement.
    out.objective = 0;
    for (std::size_t j = 0; j < n; ++j) out.objective += problem.c[j] * out.x[j];
    (void)constant;
    out.status = "optimal";
    return out;
}

FluxSolution balance(const NetworkSpec& network) {
    FluxSolution out;
    out.objectiveReactionId = "";
    out.solverStatus = "balance check only";
    for (const ReactionSpec& r : network.reactions) out.reactionIds.push_back(r.id);

    std::unordered_map<std::string, const SpeciesSpec*> byId;
    for (const SpeciesSpec& s : network.species) byId[s.id] = &s;

    bool all = true;
    for (const ReactionSpec& r : network.reactions) {
        std::map<int, double> atoms;   // atomic number -> signed count
        double charge = 0;
        bool unchecked = false;
        std::string missing;
        auto accumulate = [&](const std::vector<std::pair<std::string, double>>& side,
                              double sign) {
            for (const auto& [id, stoich] : side) {
                auto it = byId.find(id);
                if (it == byId.end()) { unchecked = true; missing = id; return; }
                const SpeciesSpec& s = *it->second;
                if (s.formula.empty()) { unchecked = true; missing = id; continue; }
                const auto parsed = chem::parseFormula(s.formula);
                if (!parsed || !parsed->ok) { unchecked = true; missing = id; continue; }
                for (const chem::FormulaTerm& t : parsed->terms)
                    atoms[t.z] += sign * stoich * t.count;
                // The species' declared charge wins over one written on the formula:
                // a reconstruction states charge as a field, and a formula that also
                // carries a sign would double-count it.
                charge += sign * stoich * (s.charge != 0 ? s.charge : parsed->charge);
            }
        };
        accumulate(r.reactants, -1.0);
        accumulate(r.products, +1.0);

        if (unchecked) {
            all = false;
            out.unbalancedReactions.push_back(r.id + ": UNCHECKED - species '" + missing +
                                              "' has no elemental formula, so this reaction is "
                                              "counted as unbalanced rather than assumed to "
                                              "balance");
            continue;
        }
        std::string offence;
        for (const auto& [z, delta] : atoms) {
            if (std::abs(delta) > 1e-9) {
                if (!offence.empty()) offence += ", ";
                offence += "element Z=" + std::to_string(z) + " off by " + fmt(delta);
            }
        }
        if (std::abs(charge) > 1e-9) {
            if (!offence.empty()) offence += ", ";
            offence += "charge off by " + fmt(charge);
        }
        if (!offence.empty()) {
            all = false;
            out.unbalancedReactions.push_back(r.id + ": " + offence);
        }
    }
    out.massBalanced = all;
    out.feasible = all;
    if (!all)
        out.warnings.push_back(
            "mass/charge balance FAILED for " + std::to_string(out.unbalancedReactions.size()) +
            " reaction(s); fba() will refuse to run until every reaction balances");
    return out;
}

namespace {

// The shared setup for every LP over a network: S v = 0 with per-reaction bounds.
struct FluxModel {
    Network                                    net;
    std::vector<double>                        lower, upper;
    std::vector<FluxBound>                     applied;
    std::size_t                                objective = 0;
    bool                                       ok = false;
    std::string                                error;
};

FluxModel prepare(const NetworkSpec& network, const std::string& objectiveReactionId,
                  const std::vector<FluxBound>& bounds) {
    FluxModel m;
    if (!Network::compile(network, m.net, &m.error)) return m;
    const std::size_t nr = m.net.reactionCount();
    m.lower.assign(nr, 0.0);
    m.upper.assign(nr, 1000.0);
    for (std::size_t r = 0; r < nr; ++r) {
        if (network.reactions[r].reversible) m.lower[r] = -1000.0;
    }
    for (const FluxBound& b : bounds) {
        const auto it = std::find(m.net.reactionIds().begin(), m.net.reactionIds().end(),
                                  b.reactionId);
        if (it == m.net.reactionIds().end()) {
            m.error = "bound names unknown reaction '" + b.reactionId + "'";
            return m;
        }
        const std::size_t idx =
            static_cast<std::size_t>(std::distance(m.net.reactionIds().begin(), it));
        m.lower[idx] = b.lower;
        m.upper[idx] = b.upper;
    }
    // Every bound actually in force is reported, not just the ones the caller passed:
    // a flux distribution without the medium it was allowed is not interpretable.
    for (std::size_t r = 0; r < nr; ++r)
        m.applied.push_back({m.net.reactionIds()[r], m.lower[r], m.upper[r]});
    const auto obj = std::find(m.net.reactionIds().begin(), m.net.reactionIds().end(),
                               objectiveReactionId);
    if (obj == m.net.reactionIds().end()) {
        m.error = "objective names unknown reaction '" + objectiveReactionId + "'";
        return m;
    }
    m.objective = static_cast<std::size_t>(std::distance(m.net.reactionIds().begin(), obj));
    m.ok = true;
    return m;
}

// S v = 0 over the NON-boundary species only: a boundary species is an exchange pool
// and must not be steady-state constrained, or every exchange reaction is forced to
// zero and the model is trivially infeasible.
LpProblem steadyState(const FluxModel& m, std::size_t extraVariables) {
    LpProblem lp;
    const std::size_t nr = m.net.reactionCount();
    std::vector<std::size_t> rows;
    for (std::size_t i = 0; i < m.net.speciesCount(); ++i)
        if (!m.net.boundary(i)) rows.push_back(i);
    lp.variables = nr + extraVariables;
    lp.constraints = rows.size();
    lp.a.assign(lp.constraints * lp.variables, 0.0);
    lp.b.assign(lp.constraints, 0.0);
    lp.c.assign(lp.variables, 0.0);
    lp.lower.assign(lp.variables, 0.0);
    lp.upper.assign(lp.variables, kInfiniteBound);
    const std::vector<double>& s = m.net.stoichiometry();
    for (std::size_t k = 0; k < rows.size(); ++k)
        for (std::size_t r = 0; r < nr; ++r)
            lp.a[k * lp.variables + r] = s[rows[k] * nr + r];
    for (std::size_t r = 0; r < nr; ++r) {
        lp.lower[r] = m.lower[r];
        lp.upper[r] = m.upper[r];
    }
    return lp;
}

FluxSolution refuseUnbalanced(const NetworkSpec& network, const std::string& objectiveId) {
    FluxSolution b = balance(network);
    b.objectiveReactionId = objectiveId;
    b.feasible = false;
    b.solverStatus = "refused: mass/charge balance did not pass";
    b.warnings.push_back(
        "fba() refused to run. A reaction that does not conserve elements and charge can carry "
        "flux from nothing, so an objective optimised over it measures the arithmetic, not the "
        "metabolism.");
    return b;
}

}  // namespace

FluxSolution fba(const NetworkSpec& network, const std::string& objectiveReactionId,
                 const std::vector<FluxBound>& bounds) {
    const FluxSolution check = balance(network);
    if (!check.massBalanced) return refuseUnbalanced(network, objectiveReactionId);

    FluxSolution out;
    out.objectiveReactionId = objectiveReactionId;
    out.massBalanced = true;
    const FluxModel m = prepare(network, objectiveReactionId, bounds);
    if (!m.ok) {
        out.solverStatus = m.error;
        out.warnings.push_back(m.error);
        return out;
    }
    out.reactionIds = m.net.reactionIds();
    out.exchangeBounds = m.applied;
    LpProblem lp = steadyState(m, 0);
    lp.c[m.objective] = -1.0;   // maximise the objective flux
    const LpResult r = solveLp(lp);
    out.solverStatus = r.status + " (dense two-phase primal simplex, " +
                       std::to_string(r.iterations) + " iterations)";
    out.feasible = r.feasible && r.bounded;
    if (!out.feasible) {
        out.warnings.push_back("the LP is " + r.status +
                               "; no flux distribution is reported, because reporting one from an "
                               + r.status + " problem would be reporting a number with no meaning");
        return out;
    }
    out.fluxes = r.x;
    out.fluxes.resize(m.net.reactionCount());
    out.objectiveValue = out.fluxes[m.objective];
    if (m.net.reactionCount() > 300)
        out.warnings.push_back(
            "this network has " + std::to_string(m.net.reactionCount()) +
            " reactions; the shipped LP is a DENSE tableau simplex and is not a genome-scale "
            "solver - treat the runtime, not the answer, as the limitation");
    return out;
}

std::vector<FluxRange> fva(const NetworkSpec& network, const std::string& objectiveReactionId,
                           const std::vector<FluxBound>& bounds, double objectiveFraction) {
    std::vector<FluxRange> out;
    const FluxSolution base = fba(network, objectiveReactionId, bounds);
    if (!base.feasible) return out;
    const FluxModel m = prepare(network, objectiveReactionId, bounds);
    if (!m.ok) return out;
    const std::size_t nr = m.net.reactionCount();
    const double floorValue = objectiveFraction * base.objectiveValue;

    for (std::size_t target = 0; target < nr; ++target) {
        FluxRange range;
        range.reactionId = m.net.reactionIds()[target];
        for (int direction = 0; direction < 2; ++direction) {
            // One extra variable: the slack on "objective >= fraction * optimum".
            LpProblem lp = steadyState(m, 1);
            const std::size_t slack = nr;
            lp.constraints += 1;
            lp.a.resize(lp.constraints * lp.variables, 0.0);
            lp.b.resize(lp.constraints, 0.0);
            const std::size_t row = lp.constraints - 1;
            lp.a[row * lp.variables + m.objective] = 1.0;
            lp.a[row * lp.variables + slack] = -1.0;
            lp.b[row] = floorValue;
            lp.lower[slack] = 0.0;
            lp.upper[slack] = kInfiniteBound;
            lp.c.assign(lp.variables, 0.0);
            lp.c[target] = direction == 0 ? 1.0 : -1.0;
            const LpResult r = solveLp(lp);
            if (!r.feasible || !r.bounded) continue;
            if (direction == 0) range.minimum = r.x[target];
            else range.maximum = r.x[target];
        }
        out.push_back(range);
    }
    return out;
}

FluxSolution parsimonious(const NetworkSpec& network, const std::string& objectiveReactionId,
                          const std::vector<FluxBound>& bounds) {
    FluxSolution out = fba(network, objectiveReactionId, bounds);
    if (!out.feasible) return out;
    const FluxModel m = prepare(network, objectiveReactionId, bounds);
    if (!m.ok) return out;
    const std::size_t nr = m.net.reactionCount();
    const double optimum = out.objectiveValue;

    // Variables: nr fluxes, then nr absolute-value surrogates t_j, then 2*nr slacks
    // for t_j - v_j >= 0 and t_j + v_j >= 0, then one slack pinning the objective.
    LpProblem lp = steadyState(m, 3 * nr);
    const std::size_t tStart = nr, sStart = 2 * nr;
    const std::size_t base = lp.constraints;
    lp.constraints += 2 * nr + 1;
    lp.a.assign(lp.constraints * lp.variables, 0.0);
    // Rebuild: steadyState already filled the first `base` rows, so redo them.
    {
        LpProblem seed = steadyState(m, 3 * nr);
        for (std::size_t i = 0; i < base; ++i)
            for (std::size_t j = 0; j < lp.variables; ++j)
                lp.a[i * lp.variables + j] = seed.a[i * seed.variables + j];
        lp.b.assign(lp.constraints, 0.0);
        for (std::size_t i = 0; i < base; ++i) lp.b[i] = seed.b[i];
        lp.lower = seed.lower;
        lp.upper = seed.upper;
    }
    for (std::size_t j = 0; j < nr; ++j) {
        const std::size_t rowA = base + 2 * j, rowB = base + 2 * j + 1;
        lp.a[rowA * lp.variables + tStart + j] = 1.0;
        lp.a[rowA * lp.variables + j] = -1.0;
        lp.a[rowA * lp.variables + sStart + 2 * j] = -1.0;
        lp.a[rowB * lp.variables + tStart + j] = 1.0;
        lp.a[rowB * lp.variables + j] = 1.0;
        lp.a[rowB * lp.variables + sStart + 2 * j + 1] = -1.0;
        lp.lower[tStart + j] = 0.0;
        lp.upper[tStart + j] = kInfiniteBound;
        lp.lower[sStart + 2 * j] = 0.0;
        lp.upper[sStart + 2 * j] = kInfiniteBound;
        lp.lower[sStart + 2 * j + 1] = 0.0;
        lp.upper[sStart + 2 * j + 1] = kInfiniteBound;
    }
    const std::size_t pinRow = lp.constraints - 1;
    lp.a[pinRow * lp.variables + m.objective] = 1.0;
    lp.b[pinRow] = optimum;
    lp.c.assign(lp.variables, 0.0);
    for (std::size_t j = 0; j < nr; ++j) lp.c[tStart + j] = 1.0;

    const LpResult r = solveLp(lp);
    if (!r.feasible || !r.bounded) {
        out.warnings.push_back("the parsimonious step is " + r.status +
                               "; the plain FBA flux is returned unchanged");
        return out;
    }
    out.fluxes.assign(r.x.begin(), r.x.begin() + static_cast<std::ptrdiff_t>(nr));
    out.objectiveValue = out.fluxes[m.objective];
    double total = 0;
    for (double v : out.fluxes) total += std::abs(v);
    out.solverStatus = "optimal, parsimonious (total |flux| = " + fmt(total) + ")";
    return out;
}

std::vector<FluxRange> deletions(const NetworkSpec& network,
                                 const std::string& objectiveReactionId,
                                 const std::vector<FluxBound>& bounds, int order) {
    std::vector<FluxRange> out;
    const FluxSolution base = fba(network, objectiveReactionId, bounds);
    if (!base.feasible) return out;
    const std::vector<std::string>& ids = base.reactionIds;

    auto knockOut = [&](const std::vector<std::string>& off) {
        std::vector<FluxBound> b = bounds;
        for (const std::string& id : off) {
            auto it = std::find_if(b.begin(), b.end(),
                                   [&](const FluxBound& x) { return x.reactionId == id; });
            if (it != b.end()) { it->lower = 0; it->upper = 0; }
            else b.push_back({id, 0.0, 0.0});
        }
        const FluxSolution s = fba(network, objectiveReactionId, b);
        return s.feasible ? s.objectiveValue : 0.0;
    };

    if (order <= 1) {
        for (const std::string& id : ids) {
            if (id == objectiveReactionId) continue;
            out.push_back({id, knockOut({id}), base.objectiveValue});
        }
        return out;
    }
    for (std::size_t i = 0; i < ids.size(); ++i) {
        if (ids[i] == objectiveReactionId) continue;
        for (std::size_t j = i + 1; j < ids.size(); ++j) {
            if (ids[j] == objectiveReactionId) continue;
            out.push_back({ids[i] + " + " + ids[j], knockOut({ids[i], ids[j]}),
                           base.objectiveValue});
        }
    }
    return out;
}

}  // namespace biocad::sim
