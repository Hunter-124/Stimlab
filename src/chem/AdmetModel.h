// chem/AdmetModel.h - structure-derived ADMET prediction model (header-only).
//
// ONE source of truth for the oral-bioavailability and stability-range
// predictions used across the app: callers include this and feed it the
// liabilities they perceive. Everything here is a
// pure function of scalar molecular descriptors plus a small set of detected
// structural liabilities - there is deliberately NO per-molecule hardcoding and
// no fixed default value: each number is computed from the inputs.
//
// The methods are mechanistic, literature-style heuristics:
//   * Absorption  - Veber/Egan-type passive permeability gated by polar surface
//                   area, H-bond donors and size.
//   * Bioavailability - F = f_abs x F_H, where F_H is the hepatic availability of
//                   the well-stirred model, F_H = Q_H / (Q_H + fu.CLint). Q_H is
//                   adult human hepatic blood flow, read from
//                   assets/packs/physiology.json (97 L/h, the same value the
//                   drug-interaction model uses) and stated as an assumption in
//                   the UI. fu.CLint is NOT predictable from structure: it is
//                   an explicit parameter, defaulting to a liability-derived
//                   ASSUMED value. That default is why catecholamines (catechol
//                   -> COMT/MAO) read as near-zero F while metabolically robust
//                   stimulants read high. The resulting score is rank-ordering
//                   only unless the caller supplies a measured CLint.
//   * Stability   - degradation-chemistry rules give an actual predicted pH window
//                   and temperature ceiling rather than a vague "narrow/broad".
//
// All functions are inline, so this adds no source file; it does take one link
// dependency, on biocad_core, for the physiological constants. That is the price
// of not hard-coding hepatic blood flow in a second place, and it is worth it: the
// literal 90 L/h that used to sit below already disagreed with the pack's 97.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include "core/Physiology.h"

namespace biocad::chem {

// Structural liabilities the PK / stability heuristics key on. Each backend fills
// these from its own perception (graph-based, from chem::detectGroups); the math
// below is identical regardless of where the flags came from.
struct PkLiabilities {
    bool catechol = false;        // two adjacent ring -OH: COMT + MAO, near-total first pass
    bool phenol = false;          // single ring -OH: phase-II (glucuronide/sulfate) first pass
    bool ester = false;           // carboxylic ester: carboxylesterase presystemic hydrolysis
    bool amide = false;           // amide: generally stable; slow hydrolysis only at pH extremes
    bool arylKetone = false;      // aryl / beta-keto carbonyl: reduction + base-catalysed condensation
    bool methylenedioxy = false;  // -O-CH2-O- arene bridge: CYP demethylenation (minor for F)
    bool maoLabileAmine = false;  // unsubstituted (no alpha-methyl) phenethylamine: MAO deamination
};

inline double admetClamp(double v, double lo, double hi) { return std::max(lo, std::min(hi, v)); }

// Small fixed-precision formatters so the rationale text is clean (no
// std::to_string trailing-zero noise) and identical across backends.
inline std::string fmt1(double v) { char b[32]; std::snprintf(b, sizeof b, "%.1f", v); return b; }
inline std::string fmt0(double v) { char b[32]; std::snprintf(b, sizeof b, "%.0f", v); return b; }

// ---------------------------------------------------------------------------
// Absorption + oral bioavailability
// ---------------------------------------------------------------------------
// Hepatic disposition assumptions. Nothing here is derivable from a structure, so
// each field is either a stated assumption or a measurement supplied by the user.
struct HepaticAssumptions {
    // Q_H: adult human hepatic blood flow, L/h, from the physiology pack. A
    // population average, not a patient. Zero when the pack is missing, which makes
    // the prediction below report the missing pack instead of inventing a flow.
    double hepaticBloodFlowLPerH = core::physiology().hepaticBloodFlowLPerH;
    // fu.CLint, L/h: unbound fraction x intrinsic clearance. Negative means
    // "derive an assumed value from the perceived structural liabilities".
    double unboundIntrinsicClearanceLPerH = -1.0;
    // True only when unboundIntrinsicClearanceLPerH came from an experiment.
    bool clIntMeasured = false;
};

struct BioavailabilityPrediction {
    double hiaPct = 0;             // fraction absorbed across the gut wall, %
    double firstPassSurvival = 0;  // F_H, 0..1, hepatic availability (well-stirred model)
    double bioavailabilityPct = 0; // oral F, % = hiaPct * F_H
    double unboundIntrinsicClearanceLPerH = 0;  // fu.CLint actually used, L/h
    double hepaticClearanceLPerH = 0;           // CL_H = Q_H.fu.CLint / (Q_H + fu.CLint)
    bool   clIntMeasured = false;  // false => the whole result is rank-ordering only
    // The Q_H actually used, so the rationale quotes the pack's value rather than a
    // number retyped into a format string.
    double hepaticBloodFlowLPerH = 0;
    std::string limitingRoute;     // dominant first-pass route (for the rationale)
};

// Fraction absorbed across the gut wall (human intestinal absorption proxy, %).
// TPSA is the primary determinant (Veber/Egan); donor count and very high MW or
// very low logP add dissolution / partition penalties.
inline double absorbedFractionPct(double tpsa, int hbd, double logP, double mw) {
    double hia = 100.0
               - std::max(0.0, tpsa - 60.0) * 0.70    // polar surface area
               - std::max(0, hbd - 3) * 4.0           // excess H-bond donors
               - std::max(0.0, mw - 480.0) * 0.05;    // size / dissolution limit
    if (logP < -1.0) hia -= (-1.0 - logP) * 6.0;      // too hydrophilic to partition
    return admetClamp(hia, 5.0, 99.0);
}

// fu.CLint (L/h) assumed from perceived structural liabilities. These coefficients
// are ordinal, not measured: they exist to rank a catechol below an amide, and they
// are expressed as a multiple of hepatic blood flow so the well-stirred model stays
// dimensionally honest. Any result built on them is Heuristic.
inline double assumedUnboundIntrinsicClearance(double logP, const PkLiabilities& L,
                                               double hepaticBloodFlowLPerH,
                                               std::string& limitingRoute) {
    double relative = 0.05;  // baseline: every drug loses a little to first pass
    limitingRoute = "minimal first-pass metabolism";
    if (L.catechol) {
        relative += 15.0;  // COMT O-methylation + MAO deamination, gut and liver
        limitingRoute = "catechol COMT/MAO first-pass metabolism";
    } else if (L.phenol) {
        relative += 0.30;  // intestinal/hepatic glucuronidation + sulfation
        limitingRoute = "phenol phase-II conjugation";
    }
    if (L.ester) {
        relative += 2.20;  // carboxylesterase hydrolysis in gut wall / liver / plasma
        if (!L.catechol) limitingRoute = "ester presystemic hydrolysis";
    }
    if (L.arylKetone) {
        relative += 0.90;  // carbonyl reduction / oxidative metabolism
        if (!L.catechol && !L.ester) limitingRoute = "carbonyl reductive metabolism";
    }
    if (L.maoLabileAmine) {
        relative += 2.50;  // monoamine oxidase deamination of the primary/secondary amine
        if (!L.catechol && !L.ester) limitingRoute = "MAO presystemic deamination";
    }
    relative += 0.08 * std::max(0.0, logP - 2.0);  // lipophilicity-driven CYP oxidation
    return relative * hepaticBloodFlowLPerH;
}

// Oral availability under the well-stirred hepatic model:
//   F_H  = Q_H / (Q_H + fu.CLint)
//   CL_H = Q_H . fu.CLint / (Q_H + fu.CLint)
// With the assumed fu.CLint above this is algebraically identical to the old
// 1/(1+burden) form - the difference is that the assumption is now named, carries
// units, and can be replaced by a measurement.
inline BioavailabilityPrediction predictBioavailability(
        double mw, double logP, double tpsa, int hbd, const PkLiabilities& L,
        const HepaticAssumptions& assume = {}) {
    BioavailabilityPrediction p;
    p.hiaPct = absorbedFractionPct(tpsa, hbd, logP, mw);

    const double qH = assume.hepaticBloodFlowLPerH;
    p.hepaticBloodFlowLPerH = qH;
    if (qH <= 0.0) {
        // No hepatic blood flow means no well-stirred model. Reporting the missing
        // pack is the honest answer; falling back to a built-in number would put an
        // uncited constant back into the one place this model must not have one.
        p.limitingRoute = "hepatic blood flow unavailable: " + core::physiology().error;
        p.firstPassSurvival = 0.0;
        p.hepaticClearanceLPerH = 0.0;
        p.bioavailabilityPct = 0.0;
        return p;
    }
    if (assume.unboundIntrinsicClearanceLPerH >= 0.0) {
        p.unboundIntrinsicClearanceLPerH = assume.unboundIntrinsicClearanceLPerH;
        p.clIntMeasured = assume.clIntMeasured;
        p.limitingRoute = assume.clIntMeasured ? "measured unbound intrinsic clearance"
                                               : "supplied unbound intrinsic clearance";
    } else {
        p.unboundIntrinsicClearanceLPerH =
            assumedUnboundIntrinsicClearance(logP, L, qH, p.limitingRoute);
        p.clIntMeasured = false;
    }

    const double fuClint = p.unboundIntrinsicClearanceLPerH;
    p.firstPassSurvival = qH / (qH + fuClint);
    p.hepaticClearanceLPerH = qH * fuClint / (qH + fuClint);
    p.bioavailabilityPct = admetClamp(p.hiaPct * p.firstPassSurvival, 1.0, 99.0);
    return p;
}

inline std::string bioavailabilityRationale(const BioavailabilityPrediction& p) {
    std::string s = "Hepatic availability under the stated assumptions: absorbed fraction (" +
                    fmt0(p.hiaPct) + "%) x F_H (" + fmt0(p.firstPassSurvival * 100.0) +
                    "%), well-stirred model with Q_H = " + fmt0(p.hepaticBloodFlowLPerH) +
                    " L/h (assets/packs/physiology.json) and fu.CLint = " +
                    fmt1(p.unboundIntrinsicClearanceLPerH) + " L/h";
    s += p.clIntMeasured ? " (measured)." : " (ASSUMED from structural liabilities - "
                                            "rank ordering only, not a percentage to quote).";
    s += " Limited by " + p.limitingRoute + ".";
    return s;
}

// ---------------------------------------------------------------------------
// Binding thermodynamics and ligand efficiency
//
// DELIBERATELY NOT WIRED TO DOCKING SCORES. AutoDock Vina's reported standard
// error is 2.85 kcal/mol (Trott & Olson 2010, PMC3041641), which is a factor of
// ~123 in Kd: converting a docked score to a nanomolar affinity manufactures
// precision that does not exist. These take a MEASURED dG or pActivity.
// ---------------------------------------------------------------------------

// R in kcal/(mol.K) (CODATA gas constant expressed in kcal).
inline constexpr double kGasConstantKcal = 1.987204259e-3;

// Kd = exp(dG / RT), Kd in mol/L, dG in kcal/mol.
inline double kdFromDeltaG(double deltaGKcalPerMol, double temperatureK = 298.15) {
    return std::exp(deltaGKcalPerMol / (kGasConstantKcal * temperatureK));
}

// dG = RT.ln(Kd), Kd in mol/L, dG in kcal/mol.
inline double deltaGFromKd(double kdMolar, double temperatureK = 298.15) {
    return kGasConstantKcal * temperatureK * std::log(kdMolar);
}

// LE = -dG / HAC, kcal/mol per heavy atom. Returns 0 for a degenerate heavy-atom count.
inline double ligandEfficiency(double deltaGKcalPerMol, int heavyAtoms) {
    return heavyAtoms > 0 ? -deltaGKcalPerMol / static_cast<double>(heavyAtoms) : 0.0;
}

// LLE = pActivity - logP (lipophilic ligand efficiency).
inline double lipophilicEfficiency(double pActivity, double logP) { return pActivity - logP; }

// LELP = logP / LE (the lipophilicity price of the efficiency).
inline double leLipophilicityPrice(double le, double logP) { return le != 0.0 ? logP / le : 0.0; }

// BEI = pActivity / (MW / 1000).
inline double bindingEfficiencyIndex(double pActivity, double molecularWeight) {
    return molecularWeight > 0.0 ? pActivity / (molecularWeight / 1000.0) : 0.0;
}

// SEI = pActivity / (TPSA / 100).
inline double surfaceEfficiencyIndex(double pActivity, double tpsa) {
    return tpsa > 0.0 ? pActivity / (tpsa / 100.0) : 0.0;
}

// 1.37 kcal/mol is RT.ln(10) at 300 K, the conversion between a pActivity unit and
// a binding free energy (Murray et al. 2014, PMC4060940). Use it to turn a measured
// pIC50/pKi into dG, never to turn a docking score into a pActivity.
inline constexpr double kRtLn10At300K = 1.37;
inline double deltaGFromPActivity(double pActivity) { return -kRtLn10At300K * pActivity; }

// ---------------------------------------------------------------------------
// STABILITY WINDOWS WERE DELETED HERE - Phase 14.
//
// This file used to carry predictPhWindow(), predictThermalWindow() and a
// shelf-life bucket. They worked by adding invented interval bounds per perceived
// functional group ("ester: 3.0-6.0", "catechol: 2.0-5.0"), averaging five 0-100
// flag counts, and mapping the average onto one of four strings such as
// "~24 months @ 25C/60%RH". Every number in that chain was authored, not measured
// or derived, and the output looked exactly like a stability study.
//
// A pH window, a temperature ceiling and a shelf life are extrapolations of MEASURED
// degradation rate constants. There is no structure-only predictor of any of them.
// They now live in sim::arrhenius(), sim::eyring(), sim::phRate() and
// sim::shelfLife(), which take the user's rate data, refuse to extrapolate from
// fewer than three temperatures, and carry a prediction interval. With no data,
// StabilityReport::shelfLife is notComputed("...") and the panel says which
// measurement is missing.
//
// Do not reintroduce a group-count stability window here. If a caller needs one
// without data, the correct answer is that it does not exist.

}  // namespace biocad::chem
