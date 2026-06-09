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
//   * Bioavailability - F = f_abs x F_fp, where F_fp is the fraction surviving
//                   gut + hepatic first pass, modelled as a saturable extraction
//                   1/(1+L) from a summed presystemic intrinsic-clearance burden
//                   L (well-stirred-model form). This is what makes catecholamines
//                   (catechol -> COMT/MAO) read as near-zero F while metabolically
//                   robust stimulants read high - instead of everything pinning at
//                   the old 95% ceiling.
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

namespace stimlab::chem {

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
struct BioavailabilityPrediction {
    double hiaPct = 0;             // fraction absorbed across the gut wall, %
    double firstPassSurvival = 0;  // 0..1, fraction escaping gut + hepatic first pass
    double bioavailabilityPct = 0; // oral F, % = hiaPct * firstPassSurvival
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

inline BioavailabilityPrediction predictBioavailability(
        double mw, double logP, double tpsa, int hbd, const PkLiabilities& L) {
    BioavailabilityPrediction p;
    p.hiaPct = absorbedFractionPct(tpsa, hbd, logP, mw);

    // Summed presystemic intrinsic-clearance burden (dimensionless, additive on a
    // log-clearance scale). Larger burden -> larger first-pass extraction.
    double burden = 0.05;  // baseline: every drug loses a little to first pass
    p.limitingRoute = "minimal first-pass metabolism";
    if (L.catechol) {
        burden += 15.0;  // COMT O-methylation + MAO deamination, gut and liver
        p.limitingRoute = "catechol COMT/MAO first-pass metabolism";
    } else if (L.phenol) {
        burden += 0.30;  // intestinal/hepatic glucuronidation + sulfation
        p.limitingRoute = "phenol phase-II conjugation";
    }
    if (L.ester) {
        burden += 2.20;  // carboxylesterase hydrolysis in gut wall / liver / plasma
        if (!L.catechol) p.limitingRoute = "ester presystemic hydrolysis";
    }
    if (L.arylKetone) {
        burden += 0.90;  // carbonyl reduction / oxidative metabolism
        if (!L.catechol && !L.ester) p.limitingRoute = "carbonyl reductive metabolism";
    }
    if (L.maoLabileAmine) {
        burden += 2.50;  // monoamine oxidase deamination of the primary/secondary amine
        if (!L.catechol && !L.ester) p.limitingRoute = "MAO presystemic deamination";
    }
    burden += 0.08 * std::max(0.0, logP - 2.0);  // lipophilicity-driven CYP oxidation

    p.firstPassSurvival = 1.0 / (1.0 + burden);
    p.bioavailabilityPct = admetClamp(p.hiaPct * p.firstPassSurvival, 1.0, 99.0);
    return p;
}

inline std::string bioavailabilityRationale(const BioavailabilityPrediction& p) {
    return "Oral F = absorbed fraction (" + fmt0(p.hiaPct) + "%) x first-pass survival (" +
           fmt0(p.firstPassSurvival * 100.0) + "%); limited by " + p.limitingRoute + ".";
}

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

}  // namespace stimlab::chem
