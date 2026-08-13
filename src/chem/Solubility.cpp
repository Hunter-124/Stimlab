#include "chem/Solubility.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>

#include "numeric/Ode.h"

namespace biocad::chem {
namespace {

// Formatting helper: a number in an assumptions string must be the real number the
// code used, so the reader can check it against the formula.
std::string num(double v, int decimals = 3) {
    std::string s(64, '\0');
    const int n = std::snprintf(s.data(), s.size(), "%.*f", decimals, v);
    s.resize(n > 0 ? static_cast<std::size_t>(n) : 0);
    return s;
}

constexpr double kPi = std::numbers::pi;

// BCS reference numbers. These are assumptions, not measurements, and every one of
// them is echoed into SolubilityReport::assumptions as a real number.
constexpr double kRefVolumeL = 0.250;       // fasted-state gastric volume
constexpr double kResidenceTimeS = 3600.0;  // small-intestinal transit, 1 h
constexpr double kIntestinalRadiusCm = 1.0;

}  // namespace

// --------------------------------------------------------------------- buffers

double vanSlykeBeta(const std::vector<BufferComponent>& components, double pH, double kw) {
    const double h = std::pow(10.0, -pH);
    double sum = kw / h + h;
    for (const auto& c : components) {
        const double ka = std::pow(10.0, -c.pKa);
        const double d = ka + h;
        sum += c.totalMolar * ka * h / (d * d);
    }
    return 2.302585092994046 * sum;
}

BufferReport bufferCapacity(const BufferSpec& spec) {
    BufferReport r;
    r.components = spec.components;

    const double step = spec.pHStep > 0.0 ? spec.pHStep : 0.02;
    double bestBeta = -1.0, bestPh = 0.0;
    for (double pH = spec.pHMin; pH <= spec.pHMax + 1e-12; pH += step) {
        const double beta = vanSlykeBeta(spec.components, pH, spec.kw);
        r.curve.push_back({pH, beta});
        if (beta > bestBeta) {
            bestBeta = beta;
            bestPh = pH;
        }
    }

    const std::string src = "Van Slyke buffer value beta = 2.303 (Kw/[H+] + [H+] + "
                            "sum C Ka [H+]/(Ka+[H+])^2), ideal activities";
    r.betaAtPh74 = makeQuantity(vanSlykeBeta(spec.components, 7.4, spec.kw), "mol/L/pH", 0.0,
                                Provenance::Model, src);
    if (r.curve.empty()) {
        r.maxCapacity = notComputed("a non-empty pH range");
        r.maxCapacityPh = notComputed("a non-empty pH range");
    } else {
        r.maxCapacity = makeQuantity(bestBeta, "mol/L/pH", 0.0, Provenance::Model, src);
        // The maximum is located on the sampled grid, so its pH is grid-resolved:
        // the half-step is a real uncertainty and is reported as the error bar.
        r.maxCapacityPh = makeQuantity(bestPh, "pH", 0.5 * step, Provenance::Model,
                                       src + "; maximum located on the sampled grid");
    }
    r.assumptions.push_back("Water ion product Kw = " + num(spec.kw, 16) +
                            " (25 C); the Kw/[H+] and [H+] terms are always included, so "
                            "unbuffered water has a nonzero buffer value.");
    r.assumptions.push_back("pH grid " + num(spec.pHMin, 2) + " to " + num(spec.pHMax, 2) +
                            " in steps of " + num(step, 4) + " pH units.");
    r.assumptions.push_back("Ideal activities (activity coefficients 1); each component is a "
                            "conjugate pair at its stated total concentration.");
    return r;
}

// ------------------------------------------------------------------ solubility

const char* gseCitation() {
    return "General Solubility Equation, log S0 = 0.5 - 0.01 (MP_C - 25) - logP "
           "(Jain & Yalkowsky revision of the Yalkowsky-Valvani equation); reported "
           "average absolute error about 0.5 log10 units on drug-like sets, which is "
           "the error bar propagated onto S0";
}

Quantity gseIntrinsicSolubility(double logP, double meltingPointC) {
    const double logS0 = 0.5 - 0.01 * (meltingPointC - 25.0) - logP;
    // DOMAIN GUARD. The GSE is linear in -logP with no upper bound, so a strongly
    // hydrophilic solute (beta-alanine, logP -3.05, MP 200 C) drives it to 63 mol/L
    // - about 5.6 kg/L, which is not a solubility. The equation was regressed on
    // sparingly soluble drug-like solids, so log S0 > 0 (above ~1 mol/L) is outside
    // the range it was fitted in, and extrapolating there is exactly the kind of
    // confident nonsense the provenance rule exists to prevent. Reported by the
    // Ionization panel's harness against beta-alanine.
    if (logS0 > 0.0) {
        return notComputed("an intrinsic solubility inside the General Solubility "
                           "Equation's fitted domain (it extrapolates to " +
                           num(std::pow(10.0, logS0), 3) +
                           " mol/L here, above the ~1 mol/L ceiling of the sparingly "
                           "soluble solids it was regressed on)");
    }
    const double s0 = std::pow(10.0, logS0);
    // 0.5 log10 units, propagated to mol/L: dS = S ln(10) d(logS).
    const double err = s0 * 2.302585092994046 * 0.5;
    return makeQuantity(s0, "mol/L", err, Provenance::Predicted, gseCitation());
}

namespace {

// Total dissolved concentration of the pH-dependent branch, mol/L. Monoprotic
// Henderson-Hasselbalch: an acid gains A-, a base gains BH+.
double phBranch(IonizationKind kind, double s0, double pKa, double pH) {
    const double h = std::pow(10.0, -pH);
    const double ka = std::pow(10.0, -pKa);
    switch (kind) {
        case IonizationKind::MonoproticAcid:  return s0 * (1.0 + ka / h);
        case IonizationKind::MonoproticBase:  return s0 * (1.0 + h / ka);
        case IonizationKind::Neutral:         break;
    }
    return s0;
}

}  // namespace

SolubilityReport phSolubility(const SolubilityInput& in) {
    SolubilityReport r;
    r.moleculeId = in.moleculeId;

    // ---- S0. A measured value wins; otherwise the GSE, which needs BOTH a melting
    // point and a logP. Neither is guessable, so their absence is named.
    bool haveS0 = false;
    double s0 = 0.0;
    if (in.hasMeasuredS0) {
        s0 = in.measuredS0Molar;
        haveS0 = true;
        r.intrinsic = makeQuantity(s0, "mol/L", 0.0, Provenance::Measured,
                                   in.measuredS0Source.empty()
                                       ? std::string("measured intrinsic solubility (user input)")
                                       : in.measuredS0Source);
    } else if (!in.hasMeltingPoint) {
        r.intrinsic = notComputed("melting point");
    } else if (!in.hasLogP) {
        r.intrinsic = notComputed("logP");
    } else {
        r.intrinsic = gseIntrinsicSolubility(in.logP, in.meltingPointC);
        if (r.intrinsic.provenance == Provenance::NotComputed) {
            // Out of the GSE's domain: no S0, therefore no curve and no BCS numbers.
            r.warnings.push_back("The General Solubility Equation was not used: " +
                                 r.intrinsic.source + ". Enter a measured intrinsic "
                                 "solubility to get the pH-solubility profile.");
        } else {
            s0 = r.intrinsic.value;
            haveS0 = true;
        }
    }

    const bool ionizable = in.kind != IonizationKind::Neutral;
    if (ionizable && !in.hasPKa) {
        r.warnings.push_back("An ionizable compound was declared without a pKa, so only the "
                             "neutral-species solubility is shown; pKa is an input, never "
                             "predicted here.");
    }
    const bool haveBranch = haveS0 && (!ionizable || in.hasPKa);

    // ---- Salt plateau. Ksp AND a counterion concentration are both required: the
    // plateau height is Ksp/[counterion], and without it there is no locatable kink.
    bool haveSalt = false;
    double saltPlateau = 0.0;
    if (!in.hasKsp) {
        r.pHmax = notComputed("salt solubility product");
    } else if (!in.hasCounterion || in.counterionMolar <= 0.0) {
        r.pHmax = notComputed("counterion concentration");
    } else if (!haveBranch) {
        r.pHmax = notComputed(haveS0 ? "pKa" : "intrinsic solubility S0");
    } else if (!ionizable) {
        r.pHmax = notComputed("an ionizable group (a neutral solid has no salt branch)");
    } else {
        saltPlateau = s0 + in.ksp / in.counterionMolar;
        haveSalt = true;
        // The two branches intersect where the ionized term equals the salt term:
        //   acid: s0 Ka/[H+] = Ksp/[Cn]  ->  [H+] = s0 Ka [Cn] / Ksp
        //   base: s0 [H+]/Ka = Ksp/[Cn]  ->  [H+] = Ka Ksp / (s0 [Cn])
        const double ka = std::pow(10.0, -in.pKa);
        const double ionized = in.ksp / in.counterionMolar;
        const double h = in.kind == IonizationKind::MonoproticAcid
                             ? s0 * ka * in.counterionMolar / in.ksp
                             : ka * ionized / s0;
        r.pHmax = makeQuantity(-std::log10(h), "pH", 0.0, Provenance::Model,
                               "intersection of the pH-dependent branch with the salt "
                               "plateau Ksp/[counterion]");
        r.assumptions.push_back("Salt plateau " + num(saltPlateau, 8) +
                                " mol/L = S0 + Ksp/[counterion] with Ksp = " + num(in.ksp, 12) +
                                " (mol/L)^2 and [counterion] = " + num(in.counterionMolar, 6) +
                                " mol/L (common-ion suppression is this input, not a default).");
    }

    // ---- Curve.
    const double step = in.pHStep > 0.0 ? in.pHStep : 0.02;
    if (haveBranch) {
        for (double pH = in.pHMin; pH <= in.pHMax + 1e-12; pH += step) {
            const double branch = phBranch(in.kind, s0, in.pKa, pH);
            const bool limited = haveSalt && branch > saltPlateau;
            const double s = limited ? saltPlateau : branch;
            r.curve.push_back({pH, std::log10(s), limited});
        }
        const double b74 = phBranch(in.kind, s0, in.pKa, 7.4);
        const bool lim74 = haveSalt && b74 > saltPlateau;
        r.solubilityAtPh74 =
            makeQuantity(lim74 ? saltPlateau : b74, "mol/L", r.intrinsic.error,
                         weakest(r.intrinsic.provenance, Provenance::Model),
                         (lim74 ? std::string("salt-limited plateau; ")
                                : std::string("Henderson-Hasselbalch pH-solubility; ")) +
                             r.intrinsic.source);
    } else {
        r.solubilityAtPh74 = notComputed(haveS0 ? "pKa" : "intrinsic solubility S0");
    }

    // ---- BCS numbers. Each is notComputed when any of its inputs is missing.
    const bool haveS74 = r.solubilityAtPh74.provenance != Provenance::NotComputed;
    const double s74 = r.solubilityAtPh74.value;

    if (!haveS74) {
        r.doseNumber = notComputed("solubility at the reference pH");
    } else if (!in.hasDose) {
        r.doseNumber = notComputed("dose");
    } else if (!in.hasMolWeight) {
        r.doseNumber = notComputed("molecular weight");
    } else {
        const double doseMolar = (in.doseMg / 1000.0 / in.molWeight) / kRefVolumeL;
        r.doseNumber = makeQuantity(doseMolar / s74, "", 0.0, Provenance::Model,
                                    "BCS dose number Do = (dose/V0)/S with V0 = " +
                                        num(kRefVolumeL, 3) + " L");
    }

    if (!haveS74) {
        r.dissolutionNumber = notComputed("solubility at the reference pH");
    } else if (!in.hasParticleRadius) {
        r.dissolutionNumber = notComputed("particle radius");
    } else if (!in.hasDiffusivity) {
        r.dissolutionNumber = notComputed("aqueous diffusivity");
    } else if (!in.hasDensity) {
        r.dissolutionNumber = notComputed("true density");
    } else if (!in.hasMolWeight) {
        r.dissolutionNumber = notComputed("molecular weight");
    } else {
        const double rCm = in.particleRadiusUm * 1.0e-4;
        const double sGPerCm3 = s74 * in.molWeight / 1000.0;   // mol/L * g/mol -> g/cm3
        const double dn = (3.0 * in.diffusivityCm2PerS /
                           (rCm * rCm * in.densityGPerCm3)) * sGPerCm3 * kResidenceTimeS;
        r.dissolutionNumber = makeQuantity(dn, "", 0.0, Provenance::Model,
                                          "BCS dissolution number Dn = (3 D / (r^2 rho)) S "
                                          "t_res with t_res = " + num(kResidenceTimeS, 1) + " s");
    }

    if (!in.hasPeff) {
        r.absorptionNumber = notComputed("effective permeability Peff");
    } else {
        r.absorptionNumber = makeQuantity(in.peffCmPerS * kResidenceTimeS / kIntestinalRadiusCm,
                                          "", 0.0, Provenance::Model,
                                          "BCS absorption number An = Peff t_res / R with "
                                          "t_res = " + num(kResidenceTimeS, 1) +
                                          " s and R = " + num(kIntestinalRadiusCm, 2) + " cm");
    }

    r.assumptions.push_back("BCS reference volume V0 = " + num(kRefVolumeL, 3) +
                            " L, small-intestinal residence time t_res = " +
                            num(kResidenceTimeS, 1) + " s, intestinal radius R = " +
                            num(kIntestinalRadiusCm, 2) +
                            " cm. These are stated assumptions, not measurements.");
    r.assumptions.push_back("Ideal activities; the pH-dependent branch is monoprotic "
                            "Henderson-Hasselbalch, S(pH) = S0 (1 + Ka/[H+]) for an acid and "
                            "S0 (1 + [H+]/Ka) for a base.");
    r.assumptions.push_back("pH grid " + num(in.pHMin, 2) + " to " + num(in.pHMax, 2) +
                            " in steps of " + num(step, 4) + " pH units.");
    if (!haveSalt) {
        r.warnings.push_back("No salt-limited plateau is drawn: pHmax is " + r.pHmax.source +
                             ", and a kink whose position is unknown is not drawn at all.");
    }
    return r;
}

// ----------------------------------------------------------------- dissolution

DissolutionReport dissolutionTimeCourse(const DissolutionInput& in) {
    DissolutionReport r;

    auto missing = [&](const char* what) {
        r.timeTo85Pct = notComputed(what);
        r.warnings.push_back(std::string("Dissolution not simulated: ") + what +
                             " is required and was not supplied.");
        return r;
    };
    if (in.doseMg <= 0.0) return missing("dose");
    if (in.molWeight <= 0.0) return missing("molecular weight");
    if (in.initialRadiusUm <= 0.0) return missing("initial particle radius");
    if (in.densityGPerCm3 <= 0.0) return missing("true density");
    if (in.diffusivityCm2PerS <= 0.0) return missing("aqueous diffusivity");
    if (in.diffusionLayerUm <= 0.0) return missing("diffusion layer thickness");
    if (in.solubilityMolar <= 0.0) return missing("solubility");
    if (in.volumeL <= 0.0) return missing("dissolution volume");
    if (in.precipitation && !in.hasKppt)
        return missing("precipitation rate constant kppt (supplied or fitted, never predicted)");

    const double rhoMg = in.densityGPerCm3 * 1000.0;          // mg/cm3
    const double volCm3 = in.volumeL * 1000.0;                // cm3
    const double r0 = in.initialRadiusUm * 1.0e-4;            // cm
    const double film = in.diffusionLayerUm * 1.0e-4;         // cm
    const double csMg = in.solubilityMolar * in.molWeight;    // mol/L * g/mol == mg/cm3
    const double spMg = in.precipSolubilityMolar * in.molWeight;
    const double particles = in.doseMg / (rhoMg * (4.0 / 3.0) * kPi * r0 * r0 * r0);

    // y = {dissolved mg, solid mg, precipitated mg}. Solid is a state and the radius
    // is derived from it as r = r0 (m/m0)^(1/3) - which IS Hixson-Crowell - so the
    // three derivatives sum to zero identically and RK4 conserves the total exactly.
    auto radiusOf = [&](double solidMg) {
        return solidMg <= 0.0 ? 0.0 : r0 * std::cbrt(solidMg / in.doseMg);
    };

    const auto f = [&](double, const std::vector<double>& y, std::vector<double>& dy) {
        const double dissolved = y[0], solid = y[1];
        const double c = dissolved / volCm3;               // mg/cm3
        const double rad = radiusOf(solid);
        double flux = 0.0;
        if (solid > 0.0 && c < csMg) {
            flux = particles * 4.0 * kPi * rad * rad * (in.diffusivityCm2PerS / film) *
                   (csMg - c);
        }
        double precip = 0.0;
        if (in.precipitation && in.hasKppt && c > spMg) {
            precip = in.kpptPerS * (c - spMg) * volCm3;
        }
        dy[0] = flux - precip;
        dy[1] = -flux;
        dy[2] = precip;
    };

    std::vector<double> y{0.0, in.doseMg, 0.0};
    const double step = in.stepS > 0.0 ? in.stepS : 0.5;
    double worst = 0.0;
    const auto observe = [&](double t, const std::vector<double>& s) {
        const double imbalance = std::abs(s[0] + s[1] + s[2] - in.doseMg);
        worst = std::max(worst, imbalance);
        // mg -> mol/L: (mg / 1000 -> g) / (g/mol) / L.
        const double molar = (s[0] / 1000.0 / in.molWeight) / in.volumeL;
        r.points.push_back({t, molar, s[1], s[2], radiusOf(s[1]) * 1.0e4});
    };
    numeric::rk4Integrate(0.0, in.horizonS, step, y, f, observe);
    r.maxMassImbalance = worst;

    // timeTo85Pct is interpolated between the two straddling observations rather than
    // snapped to the grid, so it does not silently inherit the step size.
    const double target = 0.85 * in.doseMg;
    r.timeTo85Pct = notComputed("85% never reached in the simulated horizon");
    for (std::size_t i = 1; i < r.points.size(); ++i) {
        const double prev = in.doseMg - r.points[i - 1].solidMg - r.points[i - 1].precipitatedMg;
        const double cur = in.doseMg - r.points[i].solidMg - r.points[i].precipitatedMg;
        if (cur >= target && prev < target) {
            const double frac = (target - prev) / (cur - prev);
            const double t = r.points[i - 1].timeS +
                             frac * (r.points[i].timeS - r.points[i - 1].timeS);
            r.timeTo85Pct = makeQuantity(t, "s", 0.0, Provenance::Model,
                                         "Noyes-Whitney dissolution of Hixson-Crowell "
                                         "shrinking spheres, RK4 step " + num(step, 4) + " s");
            break;
        }
    }

    r.assumptions.push_back("Monodisperse spheres: " + num(particles, 1) +
                            " particles of initial radius " + num(in.initialRadiusUm, 3) +
                            " um at density " + num(in.densityGPerCm3, 4) + " g/cm3.");
    r.assumptions.push_back("Noyes-Whitney film thickness h = " + num(in.diffusionLayerUm, 3) +
                            " um, diffusivity D = " + num(in.diffusivityCm2PerS, 8) +
                            " cm2/s, volume " + num(in.volumeL, 4) + " L, Cs = " +
                            num(in.solubilityMolar, 8) + " mol/L.");
    r.assumptions.push_back("RK4 fixed step " + num(step, 4) + " s over " +
                            num(in.horizonS, 1) + " s; worst mass imbalance " +
                            num(worst, 15) + " mg.");
    if (in.precipitation) {
        r.assumptions.push_back("pH-shift precipitation dC/dt = -kppt (C - S) with kppt = " +
                                num(in.kpptPerS, 8) + " 1/s (supplied or fitted, never "
                                "predicted) relaxing to S = " +
                                num(in.precipSolubilityMolar, 8) + " mol/L.");
    }
    return r;
}

Quantity fitPrecipitationRate(const std::vector<double>& timeS,
                              const std::vector<double>& concentrationMolar,
                              double solubilityMolar) {
    if (timeS.size() != concentrationMolar.size())
        return notComputed("time and concentration series of equal length");

    // ln(C - S) is linear in t with slope -kppt, so the fit is an exact linear
    // least squares on the supersaturated points; points at or below S carry no
    // information about kppt and are excluded rather than clamped.
    double n = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (std::size_t i = 0; i < timeS.size(); ++i) {
        const double excess = concentrationMolar[i] - solubilityMolar;
        if (excess <= 0.0) continue;
        const double x = timeS[i], y = std::log(excess);
        ++n; sx += x; sy += y; sxx += x * x; sxy += x * y;
    }
    if (n < 2.0) return notComputed("at least two supersaturated observations");
    const double denom = n * sxx - sx * sx;
    if (denom <= 0.0) return notComputed("observations at two distinct times");
    const double slope = (n * sxy - sx * sy) / denom;
    if (slope >= 0.0)
        return notComputed("a decaying supersaturation (the series does not decay)");

    // Standard error of the slope from the residual variance; with n == 2 the fit is
    // exact and the error bar is honestly 0 residual, not a fabricated interval.
    const double intercept = (sy - slope * sx) / n;
    double ss = 0.0;
    for (std::size_t i = 0; i < timeS.size(); ++i) {
        const double excess = concentrationMolar[i] - solubilityMolar;
        if (excess <= 0.0) continue;
        const double resid = std::log(excess) - (intercept + slope * timeS[i]);
        ss += resid * resid;
    }
    const double se = n > 2.0 ? std::sqrt((ss / (n - 2.0)) * n / denom) : 0.0;
    return makeQuantity(-slope, "1/s", se, Provenance::Model,
                        "kppt fitted to the supplied concentration decay against "
                        "dC/dt = -kppt (C - S) over " + num(n, 0) + " supersaturated points");
}

}  // namespace biocad::chem
