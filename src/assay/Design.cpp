#include "assay/Design.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "assay/Dataset.h"
#include "assay/Fits.h"
#include "assay/Qc.h"
#include "data/Domain.h"

namespace biocad::assay {
namespace {

// PCG64-DXSM constants. kMul is the "cheap multiplier" the DXSM variant uses on
// both the state advance and the output permutation; kIncrement is an arbitrary
// odd stream constant, fixed so the stream is a property of the seed alone.
constexpr std::uint64_t kMul = 15750249268501108917ULL;
constexpr __uint128_t   kIncrement =
    (static_cast<__uint128_t>(6364136223846793005ULL) << 64) | 1442695040888963407ULL;

constexpr double kTwoPi = 6.283185307179586476925286766559;

// Weight one design point's information contribution. Numeric derivatives keep the
// D-optimality machinery model-agnostic: the analytic 4PL Jacobian would have to be
// re-derived for every truth model this simulator can generate from, and the
// determinant ratio is insensitive to a 1e-6 relative step.
std::vector<double> sensitivity(AssayModel model, const std::vector<double>& p, double x) {
    std::vector<double> g(p.size(), 0.0);
    for (std::size_t k = 0; k < p.size(); ++k) {
        std::vector<double> lo = p, hi = p;
        const double h = std::max(1e-8, std::fabs(p[k]) * 1e-6);
        lo[k] -= h;
        hi[k] += h;
        double a = 0.0, b = 0.0;
        if (!truthValue(model, lo, x, a) || !truthValue(model, hi, x, b)) return {};
        g[k] = (b - a) / (2.0 * h);
    }
    return g;
}

// log|F'F| for a design. The logarithm rather than the determinant because a 5-
// parameter information matrix over nanomolar concentrations overflows a double's
// exponent long before the design is unreasonable.
double logDetInformation(AssayModel model, const std::vector<double>& p,
                         const std::vector<double>& design, int replicates) {
    const std::size_t np = p.size();
    if (np == 0 || design.empty()) return -std::numeric_limits<double>::infinity();
    std::vector<double> m(np * np, 0.0);
    for (double x : design) {
        const std::vector<double> g = sensitivity(model, p, x);
        if (g.size() != np) return -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < np; ++i)
            for (std::size_t j = 0; j < np; ++j)
                m[i * np + j] += static_cast<double>(replicates) * g[i] * g[j];
    }
    // Cholesky; a non-positive pivot means the design cannot identify the model,
    // which is exactly the case a D-criterion must reject rather than smooth over.
    double logDet = 0.0;
    for (std::size_t i = 0; i < np; ++i) {
        double d = m[i * np + i];
        for (std::size_t k = 0; k < i; ++k) d -= m[i * np + k] * m[i * np + k];
        if (!(d > 0.0)) return -std::numeric_limits<double>::infinity();
        const double l = std::sqrt(d);
        m[i * np + i] = l;
        for (std::size_t j = i + 1; j < np; ++j) {
            double s = m[j * np + i];
            for (std::size_t k = 0; k < i; ++k) s -= m[j * np + k] * m[i * np + k];
            m[j * np + i] = s / l;
        }
        logDet += 2.0 * std::log(l);
    }
    return logDet;
}

double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();
    return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

std::string num(double v) {
    char buf[40];
    std::snprintf(buf, sizeof buf, "%.10g", v);
    return buf;
}

std::string wellName(int row, int column) {
    std::string s;
    if (row >= 26) s.push_back(static_cast<char>('A' + row / 26 - 1));
    s.push_back(static_cast<char>('A' + row % 26));
    return s + std::to_string(column + 1);
}

// The log10 interval the fit itself reported for the half-maximal concentration.
// Preference order is deliberate: the profile-likelihood interval when the fitter
// computed one, then the Wald interval from the standard error, and NOTHING when
// neither exists - a run with no interval is a run that cannot contribute to a
// coverage count, and dropping it from the denominator would inflate coverage.
bool reportedLog10Interval(const FitResult& f, double& lo, double& hi) {
    for (const auto& p : f.parameters) {
        if (p.name != "log10EC50" && p.name != "log10C") continue;
        if (p.profileComputed && p.profileUpper > p.profileLower) {
            lo = p.profileLower;
            hi = p.profileUpper;
            return true;
        }
        if (p.value.error > 0.0) {
            lo = p.value.value - 1.959963984540054 * p.value.error;
            hi = p.value.value + 1.959963984540054 * p.value.error;
            return true;
        }
        return false;
    }
    return false;
}

}  // namespace

// ---------------------------------------------------------------------------
// DesignRng
// ---------------------------------------------------------------------------

DesignRng::DesignRng(std::uint64_t seed) {
    inc_ = (static_cast<__uint128_t>(seed) << 1) | 1;
    state_ = kIncrement + inc_;
    next();
    state_ += (static_cast<__uint128_t>(seed) << 64) | seed;
    next();
}

std::uint64_t DesignRng::next() {
    state_ = state_ * kMul + inc_;
    // DXSM output: fold the high half, multiply by the cheap multiplier, xor the
    // shifted result back in. This is the 2019 permutation, not the older XSL-RR.
    std::uint64_t hi = static_cast<std::uint64_t>(state_ >> 64);
    const std::uint64_t lo = static_cast<std::uint64_t>(state_) | 1ULL;
    hi ^= hi >> 32;
    hi *= kMul;
    hi ^= hi >> 48;
    hi *= lo;
    return hi;
}

double DesignRng::uniform() {
    return static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0);
}

double DesignRng::normal() {
    if (hasSpare_) {
        hasSpare_ = false;
        return spare_;
    }
    // Box-Muller. u1 is nudged off zero rather than rejected so the number of
    // uniforms consumed per normal is fixed at two, which is what makes a stream
    // position reproducible without tracking rejections.
    double u1 = uniform();
    if (u1 < 1e-300) u1 = 1e-300;
    const double u2 = uniform();
    const double r = std::sqrt(-2.0 * std::log(u1));
    spare_ = r * std::sin(kTwoPi * u2);
    hasSpare_ = true;
    return r * std::cos(kTwoPi * u2);
}

double DesignRng::lognormal(double sigma) {
    if (sigma <= 0.0) return 1.0;
    return std::exp(sigma * normal());
}

// ---------------------------------------------------------------------------
// Truth models
// ---------------------------------------------------------------------------

bool truthValue(AssayModel model, const std::vector<double>& p, double x, double& out) {
    switch (model) {
        case AssayModel::FourParameterLogistic:
        case AssayModel::FiveParameterLogistic: {
            const std::size_t need =
                model == AssayModel::FourParameterLogistic ? 4u : 5u;
            if (p.size() < need) return false;
            const double a = p[0], b = p[1], c = p[2], d = p[3];
            const double g = need == 5u ? p[4] : 1.0;
            if (!(c > 0.0) || !(g > 0.0)) return false;
            if (x <= 0.0) {  // the x -> 0 limit, which is the untreated signal
                out = a;
                return true;
            }
            const double u = std::pow(x / c, b);
            out = d + (a - d) / std::pow(1.0 + u, g);
            return true;
        }
        case AssayModel::MichaelisMenten: {
            if (p.size() < 2) return false;
            if (!(p[1] > 0.0)) return false;
            out = p[0] * x / (p[1] + x);
            return true;
        }
        case AssayModel::Hill: {
            if (p.size() < 3) return false;
            if (!(p[1] > 0.0)) return false;
            const double xn = std::pow(x, p[2]);
            out = p[0] * xn / (std::pow(p[1], p[2]) + xn);
            return true;
        }
        case AssayModel::SubstrateInhibition: {
            if (p.size() < 3) return false;
            if (!(p[1] > 0.0) || !(p[2] > 0.0)) return false;
            out = p[0] * x / (p[1] + x * (1.0 + x / p[2]));
            return true;
        }
        default:
            // The biophysics models are fitted from a trace over time or over
            // injections, not from a concentration ladder on a plate. Refusing is
            // correct: simulating "a Wiseman isotherm plate" would be fiction.
            return false;
    }
}

double truthEc50(AssayModel model, const std::vector<double>& p) {
    switch (model) {
        case AssayModel::FourParameterLogistic:
            return p.size() >= 4 ? p[2] : 0.0;
        case AssayModel::FiveParameterLogistic:
            if (p.size() < 5 || !(p[1] != 0.0) || !(p[4] > 0.0)) return 0.0;
            // The 5PL midpoint is NOT C. Reporting C as an EC50 is the single most
            // common 5PL error and the reason this helper exists.
            return p[2] * std::pow(std::pow(2.0, 1.0 / p[4]) - 1.0, 1.0 / p[1]);
        case AssayModel::Hill:
            return p.size() >= 3 ? p[1] : 0.0;
        case AssayModel::MichaelisMenten:
        case AssayModel::SubstrateInhibition:
            return p.size() >= 2 ? p[1] : 0.0;
        default:
            return 0.0;
    }
}

// ---------------------------------------------------------------------------
// Plate simulation
// ---------------------------------------------------------------------------

std::string simulatePlateCsv(const AssayDesignSpec& spec, std::uint64_t runSeed,
                             const std::string& plateId) {
    DesignRng rng(runSeed);

    const int rows = std::max(1, spec.rows);
    const int cols = std::max(1, spec.columns);
    const double maxConc =
        spec.concentrations.empty()
            ? 0.0
            : *std::max_element(spec.concentrations.begin(), spec.concentrations.end());
    const double pipetteSigma =
        spec.pipettingCv > 0.0 ? std::log(1.0 + spec.pipettingCv) : 0.0;

    double untreated = 0.0, plateau = 0.0;
    const bool haveControls =
        truthValue(spec.truthModel, spec.truthParameters, 0.0, untreated) &&
        truthValue(spec.truthModel, spec.truthParameters, 1e12, plateau);

    std::string csv =
        "plate_id,well,role,sample_id,series_id,concentration,conc_unit,replicate,readout,"
        "readout_unit\n";

    // One well's readout: the truth value, scaled by the plate gradient and the
    // DMSO tolerance loss, then perturbed by proportional and additive noise. The
    // order matters: a gradient is a multiplicative property of the well, and
    // detector noise is added after everything the liquid did.
    auto emit = [&](int row, int col, const char* role, const std::string& series,
                    double reportedConc, double actualConc, int replicate, double truth) {
        const double gx = cols > 1 ? static_cast<double>(col) / (cols - 1) : 0.5;
        const double gy = rows > 1 ? static_cast<double>(row) / (rows - 1) : 0.5;
        const double gradient =
            1.0 + (spec.plateGradientPct / 100.0) * (0.5 * (gx + gy) - 0.5);
        double dmso = 1.0;
        if (spec.dmsoTolerancePct > 0.0 && maxConc > 0.0 && actualConc > 0.0)
            dmso = 1.0 - (spec.dmsoTolerancePct / 100.0) * (actualConc / maxConc);
        double y = truth * gradient * dmso;
        if (spec.proportionalNoiseCv > 0.0) y *= 1.0 + spec.proportionalNoiseCv * rng.normal();
        if (spec.additiveNoiseSd > 0.0) y += spec.additiveNoiseSd * rng.normal();

        csv += plateId;
        csv += ',';
        csv += wellName(row, col);
        csv += ',';
        csv += role;
        csv += ',';
        csv += series.empty() ? std::string("control") : series;
        csv += ',';
        csv += series;
        csv += ',';
        csv += num(reportedConc);
        csv += ",M,";
        csv += std::to_string(replicate);
        csv += ',';
        csv += num(y);
        csv += ",RFU\n";
    };

    // Controls first (columns 1 and 2), so a plate that cannot hold both controls
    // and the ladder fails QC visibly instead of quietly losing the controls.
    if (haveControls && cols >= 3) {
        for (int r = 0; r < rows; ++r) {
            emit(r, 0, "negative", "", 0.0, 0.0, r, plateau);
            emit(r, 1, "positive", "", 0.0, 0.0, r, untreated);
        }
    }

    const int firstSampleCol = (haveControls && cols >= 3) ? 2 : 0;
    for (int rep = 0; rep < std::max(1, spec.replicates); ++rep) {
        // The compounding pipetting error: each transfer down the ladder multiplies
        // the accumulated factor, so the bottom of the ladder carries sqrt(n) times
        // the CV of a single transfer. That compounding IS the phenomenon being
        // simulated - drawing one independent error per well would understate the
        // uncertainty at exactly the concentrations that set the EC50.
        double accumulated = 1.0;
        std::size_t slot = 0;
        // Descending concentration is the physical order of a serial dilution.
        std::vector<double> ladder = spec.concentrations;
        std::sort(ladder.begin(), ladder.end(), std::greater<double>());
        for (double nominal : ladder) {
            accumulated *= rng.lognormal(pipetteSigma);
            const double actual = nominal * accumulated;
            double truth = 0.0;
            if (!truthValue(spec.truthModel, spec.truthParameters, actual, truth)) return {};
            const int linear = static_cast<int>(slot) +
                               rep * static_cast<int>(ladder.size());
            const int row = linear / (cols - firstSampleCol);
            const int col = firstSampleCol + linear % (cols - firstSampleCol);
            if (row >= rows) break;  // the ladder does not fit; QC will see the gap
            // The REPORTED concentration is the nominal one, because that is what
            // the experimenter writes on the plate map. The pipetting error is
            // invisible to the analysis, which is the whole reason it hurts.
            emit(row, col, "sample", "series" + std::to_string(rep + 1), nominal, actual,
                 rep, truth);
            ++slot;
        }
    }
    return csv;
}

// ---------------------------------------------------------------------------
// D-optimal ladder
// ---------------------------------------------------------------------------

DOptimalLadder dOptimalLadder(const AssayDesignSpec& spec) {
    DOptimalLadder out;
    std::vector<double> entered = spec.concentrations;
    std::sort(entered.begin(), entered.end());
    if (entered.size() < 2 || spec.truthParameters.empty()) {
        out.concentrations = entered;
        out.notes.push_back("D-optimal placement needs at least two ladder points and a "
                            "stated truth model");
        return out;
    }

    // The achievable candidate set: the entered ladder's own dilution factor, plus
    // the one intermediate point a half-step second dilution can reach. Anything
    // finer is not a serial dilution any more.
    const double top = entered.back();
    const double bottom = entered.front();
    const int k = static_cast<int>(entered.size());
    const double fold = std::pow(top / bottom, 1.0 / (k - 1));
    if (!(fold > 1.0)) {
        out.concentrations = entered;
        out.notes.push_back("the entered ladder is not monotonically diluted");
        return out;
    }
    const double halfFold = std::sqrt(fold);
    for (int j = 0; j <= 2 * (k - 1); ++j)
        out.candidates.push_back(top / std::pow(halfFold, j));
    std::sort(out.candidates.begin(), out.candidates.end());

    const double logDetEntered =
        logDetInformation(spec.truthModel, spec.truthParameters, entered, spec.replicates);

    // Fedorov coordinate exchange: one design slot at a time, swap in the candidate
    // that most improves log|F'F|, repeat until a full pass changes nothing.
    std::vector<double> design = entered;
    double best = logDetEntered;
    for (int pass = 0; pass < 64; ++pass) {
        bool improved = false;
        for (std::size_t s = 0; s < design.size(); ++s) {
            const double keep = design[s];
            double bestX = keep;
            for (double cand : out.candidates) {
                design[s] = cand;
                const double d = logDetInformation(spec.truthModel, spec.truthParameters,
                                                   design, spec.replicates);
                if (d > best + 1e-12) {
                    best = d;
                    bestX = cand;
                    improved = true;
                }
            }
            design[s] = bestX;
        }
        if (!improved) break;
    }
    std::sort(design.begin(), design.end());
    out.concentrations = design;

    const std::size_t np = spec.truthParameters.size();
    if (std::isfinite(best) && std::isfinite(logDetEntered) && np > 0)
        out.dEfficiencyGain =
            std::exp((best - logDetEntered) / static_cast<double>(np));
    else
        out.notes.push_back("one of the two designs cannot identify the model, so no "
                            "D-efficiency ratio is reported");
    out.notes.push_back("candidates are ladder points only: top / sqrt(fold)^j");
    return out;
}

// ---------------------------------------------------------------------------
// The measurement
// ---------------------------------------------------------------------------

AssayDesignReport simulateDesign(const AssayDesignSpec& spec) {
    AssayDesignReport rep;
    rep.spec = spec;
    rep.assumptions.push_back(
        "the truth model and every noise term are the USER's stated belief; this report "
        "measures that belief, not the assay");
    rep.assumptions.push_back(
        "pipetting error compounds down the ladder (lognormal per transfer) and is invisible "
        "to the analysis, exactly as on the bench");
    rep.assumptions.push_back(
        "every simulated plate is serialised to long CSV and pushed through the same "
        "import -> QC -> fit path real data takes");
    rep.assumptions.push_back(
        "the RNG is PCG64-DXSM with a Box-Muller normal transform, both implemented in "
        "BioCAD, so the same seed reproduces this report byte for byte");

    double truthMid = truthEc50(spec.truthModel, spec.truthParameters);
    double probe = 0.0;
    if (!truthValue(spec.truthModel, spec.truthParameters, 1.0, probe)) {
        rep.medianZPrime = notComputed("a truth model this simulator can generate from");
        rep.medianEc50 = notComputed("a truth model this simulator can generate from");
        rep.ec50CiWidthLog10 = notComputed("a truth model this simulator can generate from");
        rep.empiricalCoveragePct = notComputed("a truth model this simulator can generate from");
        rep.convergenceRatePct = notComputed("a truth model this simulator can generate from");
        rep.warnings.push_back(
            "the biophysics models (SPR, DSF, ITC) are fitted from a trace, not from a "
            "concentration ladder on a plate, so there is no plate to simulate");
        return rep;
    }
    if (spec.concentrations.size() < 4) {
        rep.warnings.push_back("fewer than four concentrations: a four-parameter curve cannot "
                               "be identified, so recovery here is not a design result");
    }

    const int runs = std::max(1, spec.replicateRuns);
    std::vector<double> zPrimes, widths;
    int converged = 0, covered = 0, withInterval = 0;
    FitOptions opts;
    opts.profileLikelihood = false;  // Wald from the covariance; the profile version
                                     // is a per-fit user choice, not 1000x default

    for (int run = 0; run < runs; ++run) {
        // The run seed is derived arithmetically from the spec seed, not drawn from a
        // shared stream, so run 700 is reproducible without replaying runs 0-699.
        const std::uint64_t runSeed =
            spec.seed + 0x9E3779B97F4A7C15ULL * static_cast<std::uint64_t>(run + 1);
        const std::string csv = simulatePlateCsv(spec, runSeed, "SIM" + std::to_string(run));
        if (csv.empty()) continue;

        std::string err;
        const AssayDataset ds = importText(csv, &err);
        if (ds.plates.empty()) {
            if (rep.warnings.empty() && !err.empty())
                rep.warnings.push_back("simulated plate failed import: " + err);
            continue;
        }
        const Plate& plate = ds.plates.front();

        const QcReport qc = plateQc(plate);
        if (qc.zPrime.provenance != Provenance::NotComputed) zPrimes.push_back(qc.zPrime.value);

        std::vector<Well> samples;
        for (const auto& w : plate.wells)
            if (w.role == WellRole::Sample && !w.excluded) samples.push_back(w);

        const FitResult fit = fitSeries(samples, spec.truthModel, opts);
        if (!fit.converged) continue;
        ++converged;
        if (fit.derivedEc50.provenance != Provenance::NotComputed && fit.derivedEc50.value > 0.0)
            rep.recoveredEc50.push_back(fit.derivedEc50.value);

        double lo = 0.0, hi = 0.0;
        if (reportedLog10Interval(fit, lo, hi)) {
            ++withInterval;
            widths.push_back(hi - lo);
            const double logTruth = truthMid > 0.0 ? std::log10(truthMid) : 0.0;
            if (truthMid > 0.0 && logTruth >= lo && logTruth <= hi) ++covered;
        }
    }

    const char* src = "forward simulation of the stated truth model and noise structure";
    if (!zPrimes.empty())
        rep.medianZPrime = makeQuantity(median(zPrimes), "", 0.0, Provenance::Heuristic,
                                        "Z-prime over simulated plates; bands >= 0.5 "
                                        "excellent, 0-0.5 marginal, <= 0 unusable");
    else
        rep.medianZPrime = notComputed("both a positive and a negative control on the plate");

    if (!rep.recoveredEc50.empty())
        rep.medianEc50 = makeQuantity(median(rep.recoveredEc50), "mol/L", 0.0,
                                      Provenance::Model, src);
    else
        rep.medianEc50 = notComputed("at least one converged fit");

    if (!widths.empty())
        rep.ec50CiWidthLog10 = makeQuantity(median(widths), "log10 mol/L", 0.0,
                                            Provenance::Model, src);
    else
        rep.ec50CiWidthLog10 = notComputed("a fit that reported an interval");

    rep.convergenceRatePct = makeQuantity(100.0 * converged / runs, "%", 0.0, Provenance::Model,
                                          src);
    if (withInterval > 0) {
        rep.empiricalCoveragePct = makeQuantity(
            100.0 * covered / withInterval, "%", 0.0, Provenance::Model,
            "fraction of runs whose OWN reported interval contained the truth; the "
            "nominal level is 95%");
        const double cov = rep.empiricalCoveragePct.value;
        if (cov < 90.0 || cov > 99.0)
            rep.warnings.push_back(
                "empirical coverage is " + num(cov) +
                "% against a nominal 95%: the interval this design reports does not mean "
                "what it says, and widening the ladder or the replicate count is the fix, "
                "not reporting it anyway");
    } else {
        rep.empiricalCoveragePct = notComputed("a fit that reported an EC50 interval");
    }

    const DOptimalLadder opt = dOptimalLadder(spec);
    rep.optimalConcentrations = opt.concentrations;
    rep.dOptimalityGain = opt.dEfficiencyGain;
    for (const auto& n : opt.notes) rep.assumptions.push_back(n);
    return rep;
}

// ---------------------------------------------------------------------------
// Dilution and mass arithmetic
// ---------------------------------------------------------------------------

double massForStock(double molarMassGPerMol, double molarity, double volumeL) {
    return molarMassGPerMol * molarity * volumeL;
}

DilutionPlan serialDilution(double stockMolarity, double topMolarity, double fold, int steps,
                            double wellVolumeUl, double minTransferUl) {
    DilutionPlan plan;
    plan.stockMolarity = stockMolarity;
    plan.wellVolumeUl = wellVolumeUl;
    if (!(stockMolarity > 0.0) || !(topMolarity > 0.0) || !(fold > 1.0) || steps < 1 ||
        !(wellVolumeUl > 0.0)) {
        plan.warnings.push_back("a serial ladder needs a positive stock, a top concentration, "
                                "a fold > 1, at least one step and a positive well volume");
        return plan;
    }
    if (topMolarity > stockMolarity) {
        plan.warnings.push_back("the requested top concentration exceeds the stock, so the "
                                "first well cannot be made by dilution");
        return plan;
    }
    double conc = topMolarity;
    for (int i = 0; i < steps; ++i) {
        DilutionStep s;
        s.index = i;
        s.concentration = conc;
        if (i == 0) {
            // From the stock: C_stock * V_transfer = C_top * V_well.
            s.fold = stockMolarity / topMolarity;
            s.transferVolumeUl = wellVolumeUl / s.fold;
        } else {
            s.fold = fold;
            s.transferVolumeUl = wellVolumeUl / fold;
        }
        s.diluentVolumeUl = wellVolumeUl - s.transferVolumeUl;
        if (s.transferVolumeUl < minTransferUl)
            plan.warnings.push_back(
                "step " + std::to_string(i + 1) + " transfers " + num(s.transferVolumeUl) +
                " uL, below the " + num(minTransferUl) +
                " uL floor: the pipetting CV you assumed no longer applies there");
        plan.steps.push_back(s);
        conc /= fold;
    }
    // n independent lognormal transfers compound to sqrt(n) * cv in the log domain;
    // reported per ladder rather than per well because it is the ladder's property.
    plan.compoundedCvAtLastStep = std::sqrt(static_cast<double>(steps));
    return plan;
}

}  // namespace biocad::assay
