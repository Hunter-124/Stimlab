#include "chem/Speciation.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>

#include <Eigen/Dense>

namespace biocad::chem {
namespace {

// Davies (1962): log10 gamma = -A z^2 (sqrt(I)/(1+sqrt(I)) - 0.3 I), A = 0.5085
// at 25 C in water. The 0.3 I term is what buys it the extra half-molar of range
// over Debye-Huckel, and it is also why it must not be used past I = 0.5 M: the
// term is empirical and the expression turns back up on itself.
constexpr double kDaviesA = 0.5085;
constexpr double kDaviesMaxIonicStrength = 0.5;

constexpr double kLn10 = 2.302585092994045684;
constexpr int    kMaxNewtonIterations = 200;
constexpr double kRelativeResidualTarget = 1e-12;

const char* const kCitation =
    "log Kw = -14.00 at 25 C (CODATA-consistent water ion product); "
    "activity coefficients by Davies (1962), A = 0.5085 at 25 C, valid to I = 0.5 mol/L";

std::string formatDouble(double v, int decimals) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*f", decimals, v);
    return buf;
}

// The proton component is identified by name, because the tableau is otherwise
// symmetric in its components and water's chemistry has to attach to a specific
// one. `fixedComponent` is the fallback: "at pH x" can only mean the proton.
int protonComponentIndex(const SpeciationProblem& p) {
    for (std::size_t j = 0; j < p.components.size(); ++j) {
        const std::string& n = p.components[j];
        if (n == "H" || n == "H+" || n == "H3O+") return static_cast<int>(j);
    }
    if (p.fixedComponent >= 0 && p.fixedComponent < static_cast<int>(p.components.size())) {
        return p.fixedComponent;
    }
    return -1;
}

// The caller's tableau plus water, as flat rows the solver can walk quickly.
struct Tableau {
    std::size_t                      nComponents = 0;
    std::vector<std::string>         species;
    std::vector<std::vector<double>> a;      // species x components
    std::vector<double>              logK;
    std::vector<double>              charge;
    std::vector<double>              gamma;  // activity coefficient per species
    int                              proton = -1;
    std::vector<double>              totals;
};

bool isLoneProtonRow(const std::vector<double>& row, int proton, double sign) {
    if (proton < 0) return false;
    for (std::size_t j = 0; j < row.size(); ++j) {
        const double want = (static_cast<int>(j) == proton) ? sign : 0.0;
        if (std::fabs(row[j] - want) > 1e-12) return false;
    }
    return true;
}

Tableau buildTableau(const SpeciationProblem& p, std::vector<std::string>& warnings) {
    Tableau t;
    t.nComponents = p.components.size();
    t.proton = protonComponentIndex(p);
    t.totals = p.totals;
    t.totals.resize(t.nComponents, 0.0);

    bool haveProtonSpecies = false;
    bool haveHydroxide = false;
    for (std::size_t i = 0; i < p.species.size(); ++i) {
        std::vector<double> row = i < p.stoichiometry.size() ? p.stoichiometry[i]
                                                            : std::vector<double>{};
        row.resize(t.nComponents, 0.0);
        t.species.push_back(p.species[i]);
        t.logK.push_back(i < p.logK.size() ? p.logK[i] : 0.0);
        t.charge.push_back(i < p.charges.size() ? p.charges[i] : 0.0);
        if (isLoneProtonRow(row, t.proton, 1.0)) haveProtonSpecies = true;
        if (isLoneProtonRow(row, t.proton, -1.0)) haveHydroxide = true;
        t.a.push_back(std::move(row));
    }

    if (t.proton < 0) {
        warnings.push_back(
            "no proton component named H/H+/H3O+: water self-ionization was not added and "
            "the reported pH is meaningless");
    } else {
        if (!haveProtonSpecies) {
            std::vector<double> row(t.nComponents, 0.0);
            row[static_cast<std::size_t>(t.proton)] = 1.0;
            t.species.push_back("H+");
            t.a.push_back(std::move(row));
            t.logK.push_back(0.0);
            t.charge.push_back(1.0);
        }
        if (!haveHydroxide) {
            std::vector<double> row(t.nComponents, 0.0);
            row[static_cast<std::size_t>(t.proton)] = -1.0;
            t.species.push_back("OH-");
            t.a.push_back(std::move(row));
            t.logK.push_back(kLogKwAt25C);
            t.charge.push_back(-1.0);
        }
    }

    t.gamma.assign(t.species.size(), 1.0);
    if (p.daviesActivities && p.ionicStrength > 0.0) {
        const double s = std::sqrt(p.ionicStrength);
        const double bracket = s / (1.0 + s) - 0.3 * p.ionicStrength;
        for (std::size_t i = 0; i < t.species.size(); ++i) {
            const double z = t.charge[i];
            t.gamma[i] = std::pow(10.0, -kDaviesA * z * z * bracket);
        }
    }
    return t;
}

// One evaluation of the tableau: concentrations from ln(component activities).
void concentrationsFrom(const Tableau& t, const std::vector<double>& lnx,
                        const std::vector<bool>& inactive, std::vector<double>& c) {
    c.assign(t.species.size(), 0.0);
    for (std::size_t i = 0; i < t.species.size(); ++i) {
        double lna = kLn10 * t.logK[i];
        bool   dead = false;
        for (std::size_t j = 0; j < t.nComponents; ++j) {
            const double aij = t.a[i][j];
            if (aij == 0.0) continue;
            if (inactive[j]) {  // that component's free concentration is exactly zero
                dead = true;
                break;
            }
            lna += aij * lnx[j];
        }
        if (dead) continue;
        lna = std::clamp(lna, -700.0, 700.0);
        c[i] = std::exp(lna) / t.gamma[i];
    }
}

struct Residual {
    std::vector<double> r;    // absolute, per component
    double              rel = 0.0;
};

Residual residualFrom(const Tableau& t, const std::vector<double>& c,
                      const std::vector<int>& unknownOf) {
    Residual out;
    out.r.assign(t.nComponents, 0.0);
    for (std::size_t j = 0; j < t.nComponents; ++j) {
        if (unknownOf[j] < 0) continue;
        double sum = 0.0;
        double scale = 0.0;
        for (std::size_t i = 0; i < t.species.size(); ++i) {
            const double aij = t.a[i][j];
            if (aij == 0.0) continue;
            sum += aij * c[i];
            scale += std::fabs(aij) * c[i];
        }
        out.r[j] = sum - t.totals[j];
        const double den = std::max({std::fabs(t.totals[j]), scale, 1e-300});
        out.rel = std::max(out.rel, std::fabs(out.r[j]) / den);
    }
    return out;
}

SpeciationResult finish(const SpeciationProblem& p, const Tableau& t,
                        const std::vector<double>& lnx, const std::vector<bool>& inactive,
                        const std::vector<double>& c, const Residual& res, int iterations,
                        bool converged, std::vector<std::string> warnings) {
    SpeciationResult out;
    out.concentrations = c;
    out.activityCoefficients = t.gamma;
    out.iterations = iterations;
    out.converged = converged;
    out.massBalanceResidual = res.rel;
    out.warnings = std::move(warnings);

    out.componentFree.assign(t.nComponents, 0.0);
    for (std::size_t j = 0; j < t.nComponents; ++j) {
        if (inactive[j]) continue;
        out.componentFree[j] = std::exp(std::clamp(lnx[j], -700.0, 700.0));
    }

    // Fractions are of the species' own principal component - the non-proton
    // component it draws most of its stoichiometry from - so a distribution
    // diagram's bars sum to 1 per component rather than across the whole system.
    out.fractions.assign(t.species.size(), 0.0);
    for (std::size_t i = 0; i < t.species.size(); ++i) {
        int    principal = -1;
        double best = 0.0;
        for (std::size_t j = 0; j < t.nComponents; ++j) {
            if (static_cast<int>(j) == t.proton) continue;
            const double w = std::fabs(t.a[i][j]);
            if (w > best) {
                best = w;
                principal = static_cast<int>(j);
            }
        }
        if (principal < 0) continue;
        const double total = t.totals[static_cast<std::size_t>(principal)];
        if (total > 0.0) {
            out.fractions[i] = t.a[i][static_cast<std::size_t>(principal)] * c[i] / total;
        }
    }

    double netCharge = 0.0;
    double chargeScale = 0.0;
    double ionic = 0.0;
    for (std::size_t i = 0; i < t.species.size(); ++i) {
        netCharge += t.charge[i] * c[i];
        chargeScale += std::fabs(t.charge[i]) * c[i];
        ionic += 0.5 * t.charge[i] * t.charge[i] * c[i];
    }
    out.netCharge = netCharge;
    out.chargeBalanceResidual = std::fabs(netCharge) / std::max(chargeScale, 1e-300);

    if (t.proton >= 0 && !inactive[static_cast<std::size_t>(t.proton)]) {
        out.pH = -lnx[static_cast<std::size_t>(t.proton)] / kLn10;
    }
    // The ionic strength ACTUALLY used: the caller's value when Davies is on
    // (activity coefficients were built from it), the solution's own value when
    // activities are ideal, where it is reported for information only.
    out.ionicStrength = p.daviesActivities ? p.ionicStrength : ionic;
    return out;
}

}  // namespace

const char* speciationCitation() { return kCitation; }

SpeciationResult solveSpeciation(const SpeciationProblem& p) {
    std::vector<std::string> warnings;

    if (p.daviesActivities && p.ionicStrength > kDaviesMaxIonicStrength) {
        SpeciationResult out;
        out.converged = false;
        out.ionicStrength = p.ionicStrength;
        out.warnings.push_back(
            "Davies activity model refused at ionic strength " +
            formatDouble(p.ionicStrength, 3) +
            " mol/L: its validity domain is I <= 0.5 mol/L at 25 C, and it is not "
            "extrapolated. Supply a model valid at this ionic strength or solve ideally.");
        return out;
    }

    const Tableau t = buildTableau(p, warnings);

    // A component with no total and only positive stoichiometry cannot exist in
    // solution at all; keeping it as an unknown would send ln x to -infinity.
    std::vector<bool> inactive(t.nComponents, false);
    for (std::size_t j = 0; j < t.nComponents; ++j) {
        if (t.totals[j] > 0.0) continue;
        bool signedColumn = false;
        for (std::size_t i = 0; i < t.species.size(); ++i) {
            if (t.a[i][j] < 0.0) signedColumn = true;
        }
        if (!signedColumn) inactive[j] = true;
    }

    std::vector<int> unknownOf(t.nComponents, -1);
    std::vector<std::size_t> unknowns;
    for (std::size_t j = 0; j < t.nComponents; ++j) {
        if (inactive[j]) continue;
        if (static_cast<int>(j) == p.fixedComponent) continue;
        unknownOf[j] = static_cast<int>(unknowns.size());
        unknowns.push_back(j);
    }

    // Several proton starts, because a strongly acidic or strongly basic system
    // is many decades away from neutrality and the per-step clamp is one decade.
    const double protonStarts[] = {-7.0, -3.0, -11.0, -1.0, -13.0, -5.0, -9.0};

    SpeciationResult best;
    double           bestRel = std::numeric_limits<double>::infinity();

    for (double startLogH : protonStarts) {
        std::vector<double> lnx(t.nComponents, 0.0);
        for (std::size_t j = 0; j < t.nComponents; ++j) {
            if (static_cast<int>(j) == p.fixedComponent) {
                lnx[j] = kLn10 * p.fixedLog10Activity;
            } else if (static_cast<int>(j) == t.proton) {
                lnx[j] = kLn10 * startLogH;
            } else {
                lnx[j] = std::log(std::max(t.totals[j], 1e-12) * 0.5);
            }
        }

        std::vector<double> c;
        concentrationsFrom(t, lnx, inactive, c);
        Residual res = residualFrom(t, c, unknownOf);
        int      iter = 0;
        bool     converged = res.rel < kRelativeResidualTarget;

        const std::size_t nu = unknowns.size();
        Eigen::MatrixXd   J(static_cast<Eigen::Index>(nu), static_cast<Eigen::Index>(nu));
        Eigen::VectorXd   rhs(static_cast<Eigen::Index>(nu));

        while (!converged && iter < kMaxNewtonIterations && nu > 0) {
            ++iter;
            J.setZero();
            for (std::size_t uj = 0; uj < nu; ++uj) {
                const std::size_t j = unknowns[uj];
                rhs(static_cast<Eigen::Index>(uj)) = -res.r[j];
                for (std::size_t uk = 0; uk <= uj; ++uk) {
                    const std::size_t k = unknowns[uk];
                    double sum = 0.0;
                    for (std::size_t i = 0; i < t.species.size(); ++i) {
                        const double aij = t.a[i][j];
                        const double aik = t.a[i][k];
                        if (aij == 0.0 || aik == 0.0) continue;
                        sum += aij * aik * c[i];
                    }
                    J(static_cast<Eigen::Index>(uj), static_cast<Eigen::Index>(uk)) = sum;
                    J(static_cast<Eigen::Index>(uk), static_cast<Eigen::Index>(uj)) = sum;
                }
            }

            Eigen::VectorXd step;
            Eigen::LDLT<Eigen::MatrixXd> ldlt(J);
            bool ok = ldlt.info() == Eigen::Success && ldlt.isPositive();
            if (ok) {
                step = ldlt.solve(rhs);
                ok = step.allFinite();
            }
            if (!ok) {  // a rank-deficient tableau (linearly dependent components)
                step = Eigen::FullPivLU<Eigen::MatrixXd>(J).solve(rhs);
                if (!step.allFinite()) break;
            }

            std::vector<double> d(t.nComponents, 0.0);
            for (std::size_t uj = 0; uj < nu; ++uj) {
                d[unknowns[uj]] =
                    std::clamp(step(static_cast<Eigen::Index>(uj)), -kLn10, kLn10);
            }

            double              lambda = 1.0;
            std::vector<double> trialX(t.nComponents);
            std::vector<double> trialC;
            Residual            trialRes;
            for (int back = 0; back < 7; ++back) {
                for (std::size_t j = 0; j < t.nComponents; ++j) {
                    trialX[j] = lnx[j] + lambda * d[j];
                }
                concentrationsFrom(t, trialX, inactive, trialC);
                trialRes = residualFrom(t, trialC, unknownOf);
                if (trialRes.rel < res.rel) break;
                lambda *= 0.5;
            }
            lnx = trialX;
            c = trialC;
            res = trialRes;
            converged = res.rel < kRelativeResidualTarget;
        }
        if (nu == 0) converged = true;  // everything was fixed; nothing to solve

        if (res.rel < bestRel) {
            bestRel = res.rel;
            best = finish(p, t, lnx, inactive, c, res, iter, converged, warnings);
        }
        if (converged) break;
    }

    if (!best.converged) {
        best.warnings.push_back(
            "speciation did not converge: best relative mass-balance residual " +
            formatDouble(bestRel, 3) + " after " + std::to_string(kMaxNewtonIterations) +
            " damped Newton iterations from " + std::to_string(std::size(protonStarts)) +
            " starting points; the distribution below is NOT a solution");
    }
    return best;
}

SpeciationResult solveSpeciationPh(const SpeciationProblem& p) {
    const int proton = protonComponentIndex(p);
    if (proton < 0) {
        SpeciationResult out;
        out.converged = false;
        out.warnings.push_back(
            "solveSpeciationPh needs a proton component named H/H+/H3O+ to solve for");
        return out;
    }

    SpeciationProblem q = p;
    q.fixedComponent = proton;

    const auto chargeAt = [&q](double pH, SpeciationResult& r) {
        q.fixedLog10Activity = -pH;
        r = solveSpeciation(q);
        return r.converged ? r.netCharge : std::numeric_limits<double>::quiet_NaN();
    };

    // Coarse scan then bisection: the charge-balance function is monotone in pH
    // for a real solution, but it flattens over several decades inside a buffer
    // plateau, where an unguarded Newton step leaves the bracket entirely.
    SpeciationResult probe;
    double loPh = 0.0, hiPh = 14.0;
    double loF = chargeAt(loPh, probe), hiF = 0.0;
    bool   bracketed = false;
    if (std::isfinite(loF)) {
        for (double pH = 0.25; pH <= 14.0 + 1e-9; pH += 0.25) {
            const double f = chargeAt(pH, probe);
            if (!std::isfinite(f)) continue;
            if ((loF > 0.0) != (f > 0.0)) {
                hiPh = pH;
                hiF = f;
                bracketed = true;
                break;
            }
            loPh = pH;
            loF = f;
        }
    }

    if (!bracketed) {
        SpeciationResult out = probe;
        out.converged = false;
        out.warnings.push_back(
            "charge balance has no root over pH 0-14: every counter-ion must appear as a "
            "component with its own total (e.g. Cl in 0.1 M NH4Cl), otherwise the stated "
            "composition is not electrically neutral");
        return out;
    }

    SpeciationResult inner;
    for (int i = 0; i < 200 && (hiPh - loPh) > 1e-12; ++i) {
        const double mid = 0.5 * (loPh + hiPh);
        const double f = chargeAt(mid, inner);
        if (!std::isfinite(f)) break;
        if ((f > 0.0) == (loF > 0.0)) {
            loPh = mid;
            loF = f;
        } else {
            hiPh = mid;
            hiF = f;
        }
    }
    (void)hiF;

    SpeciationResult out;
    chargeAt(0.5 * (loPh + hiPh), out);
    return out;
}

namespace {

// One microstate of the independent-group ladder: bit g set means group g holds
// its proton (HA for an acid, BH+ for a base).
std::string microstateLabel(const std::vector<IonizableGroup>& groups, unsigned mask,
                            int charge) {
    std::string sites;
    for (std::size_t g = 0; g < groups.size(); ++g) {
        if ((mask >> g) & 1u) {
            if (!sites.empty()) sites += "+";
            sites += groups[g].label;
        }
    }
    if (sites.empty()) sites = "none";
    std::string label = "H@[" + sites + "] z=";
    if (charge > 0) label += "+";
    label += std::to_string(charge);
    return label;
}

int microstateCharge(const std::vector<IonizableGroup>& groups, unsigned mask) {
    int z = 0;
    for (std::size_t g = 0; g < groups.size(); ++g) {
        const bool protonated = ((mask >> g) & 1u) != 0u;
        if (groups[g].acidic) {
            if (!protonated) --z;  // HA -> A-
        } else {
            if (protonated) ++z;   // B + H+ -> BH+
        }
    }
    return z;
}

}  // namespace

SpeciationCurve titrateGroups(const std::vector<IonizableGroup>& groups,
                              const Quantity& logP, double pHmin, double pHmax,
                              double step) {
    SpeciationCurve curve;
    curve.logP = logP;
    curve.assumptions.push_back(
        "25 C, thermodynamic (zero ionic strength) constants; every pKa is an input, "
        "never predicted");
    curve.assumptions.push_back(
        "groups titrate independently: the 2^N microstate distribution is the product of "
        "N two-state distributions, so site-site interaction (of order a pKa unit for "
        "adjacent groups) is neglected");
    curve.assumptions.push_back(
        "logD = logP + log10(f_neutral), summing every microstate of net charge zero; "
        "ion-pair partitioning of the charged microstates is not modelled");

    Provenance pKaTier = Provenance::Measured;
    for (const IonizableGroup& g : groups) {
        pKaTier = weakest(pKaTier, g.pKa.provenance);
        if (g.pKa.provenance == Provenance::NotComputed) {
            curve.isoelectricPoint = notComputed("pKa of " + g.label);
            curve.logDAtPh74 = notComputed("pKa of " + g.label);
            curve.warnings.push_back("no curve: pKa of " + g.label + " was not supplied");
            return curve;
        }
    }

    if (step <= 0.0 || pHmax < pHmin) {
        curve.isoelectricPoint = notComputed("a valid pH range and step");
        curve.logDAtPh74 = notComputed("a valid pH range and step");
        curve.warnings.push_back("pH range or step is not usable");
        return curve;
    }

    const std::size_t n = groups.size();
    const unsigned    nStates = 1u << n;
    std::vector<unsigned> order;                 // most protonated first
    order.reserve(nStates);
    for (int bits = static_cast<int>(n); bits >= 0; --bits) {
        for (unsigned mask = 0; mask < nStates; ++mask) {
            int popcount = 0;
            for (std::size_t g = 0; g < n; ++g) popcount += static_cast<int>((mask >> g) & 1u);
            if (popcount == bits) order.push_back(mask);
        }
    }
    std::vector<int> stateCharge(order.size(), 0);
    for (std::size_t s = 0; s < order.size(); ++s) {
        stateCharge[s] = microstateCharge(groups, order[s]);
        curve.labels.push_back(microstateLabel(groups, order[s], stateCharge[s]));
    }

    // Fraction of each microstate, and the two aggregates the panel plots.
    const auto evaluate = [&](double pH, std::vector<double>& fractions, double& netCharge,
                              double& fNeutral) {
        fractions.assign(order.size(), 0.0);
        netCharge = 0.0;
        fNeutral = 0.0;
        for (std::size_t s = 0; s < order.size(); ++s) {
            double f = 1.0;
            for (std::size_t g = 0; g < n; ++g) {
                // Protonated fraction of one site: 1/(1 + 10^(pH - pKa)).
                const double fH = 1.0 / (1.0 + std::pow(10.0, pH - groups[g].pKa.value));
                f *= ((order[s] >> g) & 1u) ? fH : (1.0 - fH);
            }
            fractions[s] = f;
            netCharge += f * stateCharge[s];
            if (stateCharge[s] == 0) fNeutral += f;
        }
    };

    std::vector<double> fractions;
    double              netCharge = 0.0;
    double              fNeutral = 0.0;
    for (double pH = pHmin; pH <= pHmax + 1e-9; pH += step) {
        evaluate(pH, fractions, netCharge, fNeutral);
        SpeciationPoint pt;
        pt.pH = pH;
        pt.microspeciesFractions = fractions;
        pt.netCharge = netCharge;
        pt.logD = logP.provenance == Provenance::NotComputed
                      ? 0.0
                      : logP.value + std::log10(std::max(fNeutral, 1e-300));
        curve.points.push_back(std::move(pt));
    }

    const auto chargeAtPh = [&](double pH) {
        std::vector<double> f;
        double              z = 0.0, fn = 0.0;
        evaluate(pH, f, z, fn);
        return z;
    };

    const double zLo = chargeAtPh(pHmin);
    const double zHi = chargeAtPh(pHmax);
    if ((zLo > 0.0) == (zHi > 0.0)) {
        curve.isoelectricPoint = notComputed("a pH of zero net charge within the scanned range");
    } else {
        double lo = pHmin, hi = pHmax;
        const bool loPositive = zLo > 0.0;
        while ((hi - lo) > 1e-9) {
            const double mid = 0.5 * (lo + hi);
            if ((chargeAtPh(mid) > 0.0) == loPositive) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        curve.isoelectricPoint =
            makeQuantity(0.5 * (lo + hi), "pH", 0.0, pKaTier,
                         std::string("zero net charge of the independent-group ladder; ") +
                             kCitation);
    }

    if (logP.provenance == Provenance::NotComputed) {
        curve.logDAtPh74 = notComputed("logP");
    } else {
        evaluate(7.4, fractions, netCharge, fNeutral);
        curve.logDAtPh74 = makeQuantity(
            logP.value + std::log10(std::max(fNeutral, 1e-300)), "log10 D", logP.error,
            weakest(pKaTier, logP.provenance),
            "logP (" + logP.source + ") plus log10 of the neutral microspecies fraction");
    }
    return curve;
}

}  // namespace biocad::chem
