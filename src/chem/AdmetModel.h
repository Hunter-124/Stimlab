// chem/AdmetModel.h - structure-derived ADMET prediction model (header-only).
//
// ONE source of truth for the oral-bioavailability and stability-range
// predictions used by every backend (the real chem engine and the fakes both
// include this and feed it the liabilities they perceive). Everything here is a
// pure function of scalar molecular descriptors plus a small set of detected
// structural liabilities - there is deliberately NO per-molecule hardcoding and
// no fixed default value: each number is computed from the inputs.
//
// The methods are mechanistic, literature-style heuristics:
//   * Absorption  - Veber/Egan-type passive permeability gated by polar surface
//                   area, H-bond donors and size.
//   * Bioavailability - F = f_abs x F_H, where F_H is the hepatic availability of
//                   the well-stirred model, F_H = Q_H / (Q_H + fu.CLint). Q_H is
//                   adult human hepatic blood flow (90 L/h, an assumption stated
//                   in the UI). fu.CLint is NOT predictable from structure: it is
//                   an explicit parameter, defaulting to a liability-derived
//                   ASSUMED value. That default is why catecholamines (catechol
//                   -> COMT/MAO) read as near-zero F while metabolically robust
//                   stimulants read high. The resulting score is rank-ordering
//                   only unless the caller supplies a measured CLint.
//   * Stability   - degradation-chemistry rules give an actual predicted pH window
//                   and temperature ceiling rather than a vague "narrow/broad".
//
// Header-only (all functions inline) so it needs no library target and adds no
// link dependency; it includes only the C++ standard library.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

namespace biocad::chem {

// Structural liabilities the PK / stability heuristics key on. Each backend fills
// these from its own perception (graph-based in the real engine, motif-based in
// the fakes); the math below is identical regardless of where the flags came from.
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
    // Q_H: adult human hepatic blood flow, L/h. A population average, not a patient.
    double hepaticBloodFlowLPerH = 90.0;
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
                    "%), well-stirred model with Q_H = 90 L/h and fu.CLint = " +
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
// pH stability window
// ---------------------------------------------------------------------------
struct PhWindow {
    double low = 2.0;
    double high = 10.0;
    double optimal = 6.0;
    std::string mechanism = "no acid/base- or oxidation-labile groups; broad stable window.";
};

inline PhWindow predictPhWindow(const PkLiabilities& L) {
    PhWindow w;
    std::string why;
    bool labile = false;
    auto add = [&](double lo, double hi, const char* reason) {
        w.low = std::max(w.low, lo);
        w.high = std::min(w.high, hi);
        if (!why.empty()) why += " ";
        why += reason;
        labile = true;
    };
    if (L.ester)        add(3.0, 6.0, "ester hydrolysis is acid- and base-catalysed (slowest mildly acidic).");
    if (L.catechol)     add(2.0, 5.0, "catechol auto-oxidation accelerates above ~pH 5 (keep acidic, exclude air/light).");
    else if (L.phenol)  add(3.0, 8.0, "phenol oxidation/ionisation accelerates above ~pH 8.");
    if (L.arylKetone)   add(3.0, 7.0, "alpha-keto enolisation/condensation is base-catalysed.");
    if (L.amide && !labile) add(3.0, 9.0, "amide hydrolyses only at pH extremes.");

    if (w.high < w.low + 1.0) w.high = w.low + 1.0;  // keep a sane minimum width
    w.optimal = labile ? (w.low + 0.30 * (w.high - w.low)) : 0.5 * (w.low + w.high);
    if (labile) w.mechanism = why;
    return w;
}

inline std::string phWindowText(const PhWindow& w) {
    return "Predicted stable window pH " + fmt1(w.low) + "-" + fmt1(w.high) +
           " (optimal ~pH " + fmt1(w.optimal) + "); " + w.mechanism;
}

// ---------------------------------------------------------------------------
// Temperature stability window
// ---------------------------------------------------------------------------
struct ThermalWindow {
    double stableToC = 40.0;    // accelerated-stability ceiling (~ICH 40C/75%RH)
    double storeBelowC = 25.0;  // recommended storage temperature
    bool refrigerate = false;
    std::string mechanism = "thermally robust scaffold; no low-barrier degradation route.";
};

inline ThermalWindow predictThermalWindow(const PkLiabilities& L) {
    ThermalWindow t;
    std::string why;
    bool any = false;
    auto cap = [&](double toC, double storeC, bool fridge, const char* reason) {
        t.stableToC = std::min(t.stableToC, toC);
        t.storeBelowC = std::min(t.storeBelowC, storeC);
        if (fridge) t.refrigerate = true;
        if (!why.empty()) why += " ";
        why += reason;
        any = true;
    };
    if (L.ester)        cap(30.0, 25.0, false, "ester hydrolysis accelerates with heat and humidity.");
    if (L.arylKetone)   cap(30.0, 25.0, false, "beta-keto condensation is thermally driven.");
    if (L.catechol)     cap(25.0,  8.0, true,  "catechol auto-oxidation is strongly temperature-dependent.");
    else if (L.phenol)  cap(35.0, 25.0, false, "phenol oxidation is mildly temperature-sensitive.");

    if (any) t.mechanism = why;
    return t;
}

inline std::string thermalWindowText(const ThermalWindow& t) {
    std::string s = "Predicted stable to ~" + fmt0(t.stableToC) + "C; store below ~" +
                    fmt0(t.storeBelowC) + "C";
    if (t.refrigerate) s += " (refrigerate)";
    s += "; " + t.mechanism;
    return s;
}

}  // namespace biocad::chem
