#include "fakes/FakeBackend.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <string_view>
#include <vector>

#include "chem/AdmetModel.h"

namespace biocad {
namespace {

bool contains(std::string_view hay, std::string_view needle) {
    return hay.find(needle) != std::string_view::npos;
}

double clampd(double v, double lo, double hi) { return std::max(lo, std::min(hi, v)); }

// ------------------------------------------------------------------ defaults
// Curated reference library. Property values are approximate public-data values,
// sufficient for the heuristic models below. Extend freely.
std::vector<Molecule> buildLibrary() {
    return {
        {"caffeine", "Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "C8H10N4O2",
         194.19, -0.07, 58.4, 0, 6, 0, "Xanthine stimulant", "Unscheduled (US)",
         "Adenosine receptor antagonist; common reference stimulant."},
        {"amphetamine", "Amphetamine", "CC(N)Cc1ccccc1", "C9H13N",
         135.21, 1.76, 26.0, 1, 1, 2, "Phenethylamine stimulant", "Schedule II (US)",
         "Releaser/reuptake inhibitor at DAT/NET."},
        {"methamphetamine", "Methamphetamine", "CNC(C)Cc1ccccc1", "C10H15N",
         149.23, 2.07, 12.0, 1, 1, 2, "Phenethylamine stimulant", "Schedule II (US)",
         "N-methyl amphetamine; higher CNS penetration."},
        {"mdma", "MDMA", "CNC(C)Cc1ccc2OCOc2c1", "C11H15NO2",
         193.24, 1.06, 30.0, 1, 3, 3, "Entactogen (phenethylamine)", "Schedule I (US)",
         "3,4-methylenedioxy substitution; serotonergic releaser."},
        {"methylphenidate", "Methylphenidate", "COC(=O)C(c1ccccc1)C1CCCCN1", "C14H19NO2",
         233.31, 2.30, 38.3, 1, 3, 4, "Piperidine stimulant", "Schedule II (US)",
         "DAT/NET reuptake inhibitor; methyl ester (hydrolysis-labile)."},
        {"modafinil", "Modafinil", "NC(=O)CS(=O)C(c1ccccc1)c1ccccc1", "C15H15NO2S",
         273.35, 1.40, 79.1, 1, 3, 4, "Eugeroic (wakefulness)", "Schedule IV (US)",
         "Atypical DAT inhibitor; sulfinyl + primary amide."},
        {"cocaine", "Cocaine", "COC(=O)C1C(OC(=O)c2ccccc2)CC3CCC1N3C", "C17H21NO4",
         303.35, 2.30, 55.8, 0, 5, 5, "Tropane stimulant", "Schedule II (US)",
         "Two ester groups -> hydrolysis to benzoylecgonine/ecgonine."},
        {"pseudoephedrine", "Pseudoephedrine", "CNC(C)C(O)c1ccccc1", "C10H15NO",
         165.23, 0.90, 32.3, 2, 2, 3, "Phenethylamine decongestant", "OTC/Schedule V (US)",
         "Beta-hydroxy amphetamine analog; sympathomimetic."},
        {"ephedrine", "Ephedrine", "CNC(C)C(O)c1ccccc1", "C10H15NO",
         165.23, 1.00, 32.3, 2, 2, 3, "Phenethylamine sympathomimetic", "Schedule (varies)",
         "Diastereomer of pseudoephedrine."},
        {"cathinone", "Cathinone", "CC(N)C(=O)c1ccccc1", "C9H11NO",
         149.19, 0.60, 43.1, 1, 2, 2, "Beta-keto phenethylamine", "Schedule I (US)",
         "Beta-keto group reduces stability (oxidation/condensation)."},
        {"methcathinone", "Methcathinone", "CNC(C)C(=O)c1ccccc1", "C10H13NO",
         163.22, 0.90, 29.1, 1, 2, 2, "Cathinone stimulant", "Schedule I (US)",
         "N-methyl cathinone."},
        {"mephedrone", "Mephedrone (4-MMC)", "CNC(C)C(=O)c1ccc(C)cc1", "C11H15NO",
         177.24, 1.60, 29.1, 1, 2, 3, "Substituted cathinone", "Schedule I (US)",
         "4-methyl methcathinone; 'bath salts' class."},
        {"mda", "MDA", "CC(N)Cc1ccc2OCOc2c1", "C10H13NO2",
         179.22, 1.00, 47.6, 1, 3, 2, "Entactogen (phenethylamine)", "Schedule I (US)",
         "3,4-methylenedioxyamphetamine; MDMA's primary metabolite."},
        {"methylone", "Methylone (bk-MDMA)", "CNC(C)C(=O)c1ccc2OCOc2c1", "C11H13NO3",
         207.23, 0.70, 56.4, 1, 4, 3, "Substituted cathinone", "Schedule I (US)",
         "Beta-keto MDMA analog."},
        {"mdpv", "MDPV", "O=C(C(CCC)N1CCCC1)c1ccc2OCOc2c1", "C16H21NO3",
         275.34, 2.80, 47.6, 0, 4, 4, "Substituted cathinone", "Schedule I (US)",
         "Potent DAT/NET inhibitor; pyrrolidinophenone."},
        {"alpha-pvp", "alpha-PVP", "O=C(C(CCC)N1CCCC1)c1ccccc1", "C15H21NO",
         231.33, 3.20, 20.3, 0, 2, 4, "Pyrrolidinophenone stimulant", "Schedule I (US)",
         "'Flakka'; MDPV des-methylenedioxy analog."},
        {"4-fa", "4-Fluoroamphetamine", "CC(N)Cc1ccc(F)cc1", "C9H12FN",
         153.20, 1.80, 26.0, 1, 2, 2, "Phenethylamine stimulant", "Schedule (varies)",
         "Para-fluoro amphetamine; releaser."},
        {"dmaa", "Methylhexanamine (DMAA)", "CCCC(C)CC(C)N", "C7H17N",
         115.22, 2.00, 26.0, 1, 1, 4, "Aliphatic amine stimulant", "Banned supplement (US)",
         "Non-aromatic sympathomimetic; no arene for CYP2D6/MAO motifs."},
        {"bupropion", "Bupropion", "CC(NC(C)(C)C)C(=O)c1cccc(Cl)c1", "C13H18ClNO",
         239.74, 3.20, 29.1, 1, 2, 3, "Aminoketone (NDRI)", "Rx (US)",
         "Antidepressant/cessation aid; cathinone-like aminoketone."},
        {"phenethylamine", "Phenethylamine (parent)", "NCCc1ccccc1", "C8H11N",
         121.18, 1.40, 26.0, 1, 1, 2, "Phenethylamine (trace amine)", "Endogenous",
         "Structural parent of the amphetamine class; TAAR1 agonist."},
        {"tyramine", "Tyramine", "NCCc1ccc(O)cc1", "C8H11NO",
         137.18, 0.90, 46.2, 2, 2, 2, "Phenethylamine (trace amine)", "Endogenous",
         "Dietary; MAO substrate (MAOI 'cheese' interaction)."},
        {"theobromine", "Theobromine", "Cn1cnc2c1c(=O)[nH]c(=O)n2C", "C7H8N4O2",
         180.16, -0.80, 67.2, 1, 6, 0, "Xanthine stimulant", "Unscheduled (US)",
         "Cocoa alkaloid; milder than caffeine."},
        {"theophylline", "Theophylline", "Cn1c(=O)c2[nH]cnc2n(C)c1=O", "C7H8N4O2",
         180.16, -0.02, 69.3, 1, 6, 0, "Xanthine bronchodilator", "Rx (US)",
         "Narrow therapeutic index; CYP1A2 substrate."},
        {"norepinephrine", "Norepinephrine", "NCC(O)c1ccc(O)c(O)c1", "C8H11NO3",
         169.18, -1.20, 86.7, 4, 4, 2, "Catecholamine neurotransmitter", "Endogenous",
         "Catechol; COMT/MAO substrate; poor oral/CNS absorption."},
        {"epinephrine", "Epinephrine", "CNCC(O)c1ccc(O)c(O)c1", "C9H13NO3",
         183.20, -1.40, 75.7, 3, 4, 2, "Catecholamine hormone", "Rx (US)",
         "Adrenaline; catechol; not orally bioavailable."},
        {"lisdexamfetamine", "Lisdexamfetamine", "CC(Cc1ccccc1)NC(=O)C(N)CCCCN", "C15H25N3O",
         263.38, 0.80, 95.0, 3, 4, 8, "Amphetamine prodrug", "Schedule II (US)",
         "L-lysine amide prodrug of dexamphetamine; amide hydrolysis activates."},
        {"atomoxetine", "Atomoxetine (reference)", "CNCCC(Oc1ccccc1C)c1ccccc1", "C17H21NO",
         255.35, 3.90, 21.3, 1, 2, 5, "Selective NRI (non-stimulant)", "Rx (US)",
         "Non-controlled ADHD reference; NET-selective."},
        {"dopamine", "Dopamine", "NCCc1ccc(O)c(O)c1", "C8H11NO2",
         153.18, -0.98, 66.5, 3, 3, 2, "Catecholamine neurotransmitter", "Endogenous",
         "Catechol -> rapid autoxidation; COMT/MAO substrate."},
        {"nicotine", "Nicotine", "CN1CCCC1c1cccnc1", "C10H14N2",
         162.23, 1.17, 16.1, 0, 2, 1, "Pyridine alkaloid stimulant", "Unscheduled (US)",
         "nAChR agonist."},
        {"phentermine", "Phentermine", "CC(C)(N)Cc1ccccc1", "C10H15N",
         149.23, 1.90, 26.0, 1, 1, 1, "Phenethylamine anorectic", "Schedule IV (US)",
         "Alpha,alpha-dimethyl amphetamine; resists MAO."},
        {"acetaminophen", "Acetaminophen (reference)", "CC(=O)Nc1ccc(O)cc1", "C8H9NO2",
         151.16, 0.46, 49.3, 2, 2, 1, "Analgesic (non-stimulant reference)", "OTC (US)",
         "Metabolism reference: bioactivation to reactive NAPQI."},
    };
}

// ----------------------------------------------------------- feature sniffing
struct Features {
    bool ester;
    bool amide;
    bool arylKetone;
    bool catechol;
    bool phenol;
    bool phenethylamine;
    bool methylenedioxy;
};

Features sniff(const Molecule& m) {
    Features f{};
    f.ester          = contains(m.smiles, "OC(=O)") || contains(m.smiles, "C(=O)OC");
    f.amide          = contains(m.smiles, "NC(=O)") || contains(m.smiles, "C(=O)N");
    f.arylKetone     = contains(m.smiles, "C(=O)c") || contains(m.smiles, "C(=O)C1CC");
    f.catechol       = contains(m.smiles, "c(O)c(O)");
    f.phenol         = contains(m.smiles, "c(O)") || contains(m.smiles, "(O)c");
    f.methylenedioxy = contains(m.smiles, "OCOc");
    f.phenethylamine = contains(m.drugClass, "henethylamine") ||
                       contains(m.drugClass, "atecholamine") ||
                       contains(m.drugClass, "athinone") ||
                       contains(m.drugClass, "ntactogen");
    return f;
}

// Map perceived motifs onto the shared PK / stability model's liability flags.
// maoLabileAmine: a phenethylamine whose amine sits on an unsubstituted -CH2-CH2-
// chain (the "NCC" motif) is a monoamine-oxidase substrate; an alpha-methylated
// amine (amphetamine-type, "C(N)" / "NC(C)") is MAO-resistant and is not flagged.
chem::PkLiabilities liabilitiesOf(const Molecule& m, const Features& f) {
    chem::PkLiabilities L;
    L.catechol = f.catechol;
    L.phenol = f.phenol;
    L.ester = f.ester;
    L.amide = f.amide;
    L.arylKetone = f.arylKetone;
    L.methylenedioxy = f.methylenedioxy;
    L.maoLabileAmine = f.phenethylamine && contains(m.smiles, "NCC");
    return L;
}

// ----------------------------------------------------------------- stability
StabilityReport computeStability(const Molecule& m) {
    const Features f = sniff(m);
    StabilityReport r;
    r.moleculeId = m.id;

    const double hydrolysis = f.ester ? 34.0 : (f.amide ? 70.0 : 92.0);
    const double oxidation  = f.catechol ? 22.0 : (f.phenol ? 58.0 : (m.logP > 2.0 ? 76.0 : 85.0));
    const double photolysis = (f.catechol || f.arylKetone) ? 54.0 : 82.0;
    const double thermal    = f.arylKetone ? 66.0 : (f.ester ? 72.0 : 88.0);
    const double pHsens     = (f.ester || f.arylKetone || f.catechol) ? 58.0 : 80.0;

    // Predicted, numeric stability ranges from the shared degradation-chemistry
    // model: an actual pH window and temperature ceiling instead of "narrow/broad".
    const chem::PkLiabilities sl = liabilitiesOf(m, f);
    const auto ph = chem::predictPhWindow(sl);
    const auto th = chem::predictThermalWindow(sl);

    r.factors = {
        {"Hydrolysis resistance", hydrolysis,
         f.ester ? "Ester present - susceptible to acid/base/enzymatic hydrolysis."
                 : (f.amide ? "Amide present - slow hydrolysis under stress."
                            : "No hydrolyzable carbonyl ester; resistant.")},
        {"Oxidation resistance", oxidation,
         f.catechol ? "Catechol - rapid autoxidation to ortho-quinone."
                    : (f.phenol ? "Phenol - moderate oxidative liability."
                                : "Limited oxidation-prone motifs.")},
        {"Photostability", photolysis,
         (f.catechol || f.arylKetone) ? "Chromophore/carbonyl - photodegradation risk."
                                      : "No strong photolabile chromophore."},
        {"Thermal stability (range)", thermal, chem::thermalWindowText(th)},
        {"pH stability (range)", pHsens, chem::phWindowText(ph)},
    };

    double sum = 0;
    for (const auto& fac : r.factors) sum += fac.score;
    r.overallScore = sum / static_cast<double>(r.factors.size());

    if (f.ester) r.degradants.push_back({"Hydrolyzed acid + alcohol", "Ester hydrolysis",
                                         "Primary degradation route under humidity."});
    if (f.amide) r.degradants.push_back({"Parent amine + carboxylic acid", "Amide hydrolysis",
                                         "Slow; accelerated at extreme pH."});
    if (f.catechol) r.degradants.push_back({"ortho-Quinone / melanin-like polymers", "Autoxidation",
                                            "Air/light accelerated; color change."});
    if (f.arylKetone) r.degradants.push_back({"Reduced alcohol / aldol products", "Carbonyl reactivity",
                                              "Condensation under heat."});
    if (r.degradants.empty()) r.degradants.push_back({"N-oxide (minor)", "N-oxidation",
                                                      "Trace; not shelf-life limiting."});

    const std::string store = th.refrigerate
        ? ("<=" + chem::fmt0(th.storeBelowC) + "C, dark")
        : (chem::fmt0(th.storeBelowC) + "C/60%RH");
    const std::string months = r.overallScore >= 85 ? "~36 months"
                             : r.overallScore >= 70 ? "~24 months"
                             : r.overallScore >= 55 ? "~12 months"
                                                    : "~6 months";
    r.shelfLifeEstimate = months + " @ " + store;

    r.summary = "Overall stability " + std::to_string(static_cast<int>(r.overallScore)) +
                "/100. Limiting factor: " +
                (f.ester ? "ester hydrolysis."
                         : (f.catechol ? "catechol oxidation."
                                       : (f.arylKetone ? "carbonyl reactivity." : "none major."))) +
                " Predicted stable pH " + chem::fmt1(ph.low) + "-" + chem::fmt1(ph.high) +
                ", to ~" + chem::fmt0(th.stableToC) + "C (store " + store + ").";
    return r;
}

// ---------------------------------------------------------------- absorption
AbsorptionReport computeAbsorption(const Molecule& m) {
    AbsorptionReport r;
    r.moleculeId = m.id;
    const double tpsa = m.tpsa, logP = m.logP, mw = m.molWeight;

    // Oral F is computed mechanistically (absorbed fraction x first-pass survival),
    // not pinned to a 95% ceiling: structural metabolic liabilities set the
    // presystemic extraction, so catecholamines read near-zero F and metabolically
    // robust stimulants read high.
    const Features f = sniff(m);
    const chem::PkLiabilities L = liabilitiesOf(m, f);
    const auto bio = chem::predictBioavailability(mw, logP, tpsa, m.hbd, L);

    r.hiaPct           = bio.hiaPct;
    r.bioavailabilityPct = bio.bioavailabilityPct;
    r.caco2LogPapp     = clampd(-4.80 + 0.30 * logP - 0.008 * tpsa, -7.0, -4.0);
    r.logBB            = clampd(0.15 + 0.17 * logP - 0.011 * tpsa, -2.0, 1.2);
    r.logS             = clampd(0.80 - 0.0100 * mw - 0.55 * logP, -6.5, 1.0);
    r.pgpSubstrate     = (mw > 400.0 && m.hba >= 6) || tpsa > 120.0;
    r.cnsPenetrant     = (r.logBB > -0.30) && (tpsa < 90.0) && (mw < 450.0);

    auto band = [](double v, double good, double ok) {
        return v >= good ? Verdict::Good : (v >= ok ? Verdict::Info : Verdict::Warn);
    };

    r.metrics = {
        {"Human intestinal absorption", r.hiaPct, "%", band(r.hiaPct, 80, 50),
         "Fraction crossing the gut wall; TPSA / H-bond-donor limited (passive permeability)."},
        {"Hepatic availability (F, assumed CLint)", r.bioavailabilityPct, "%",
         band(r.bioavailabilityPct, 70, 40),
         chem::bioavailabilityRationale(bio)},
        {"Caco-2 permeability", r.caco2LogPapp, "log(cm/s)", band(r.caco2LogPapp, -5.0, -6.0),
         "Higher (less negative) = more permeable monolayer flux."},
        {"BBB partition", r.logBB, "logBB", r.cnsPenetrant ? Verdict::Good : Verdict::Info,
         r.cnsPenetrant ? "CNS-penetrant - relevant for central stimulant action."
                        : "Limited CNS penetration predicted."},
        {"Aqueous solubility", r.logS, "logS", band(r.logS, -4.0, -5.0),
         "Low solubility can dissolution-limit absorption."},
        {"P-gp efflux substrate", r.pgpSubstrate ? 1.0 : 0.0, "bool",
         r.pgpSubstrate ? Verdict::Warn : Verdict::Good,
         r.pgpSubstrate ? "Efflux may reduce net absorption / CNS exposure."
                        : "Not a likely P-gp substrate."},
    };

    r.summary = "Hepatic availability under the stated assumptions ~" +
                chem::fmt0(r.bioavailabilityPct) + "% (HIA ~" + chem::fmt0(r.hiaPct) +
                "%, F_H ~" + chem::fmt0(bio.firstPassSurvival * 100.0) +
                "% from an ASSUMED fu.CLint; " + bio.limitingRoute + "). " +
                (r.cnsPenetrant ? "CNS-penetrant." : "Low CNS penetration.");
    return r;
}

// --------------------------------------------------------------------- admet
AdmetReport computeAdmet(const Molecule& m) {
    const Features f = sniff(m);
    AdmetReport r;
    r.moleculeId = m.id;

    auto add = [&](std::string name, Verdict v, std::string detail) {
        r.endpoints.push_back({std::move(name), v, std::move(detail)});
    };

    if (f.phenethylamine && !f.catechol) {
        add("MAO-A/B substrate", Verdict::Warn,
            "Phenethylamine core - monoamine oxidase metabolism; tyramine/MAOI DDI risk.");
        add("CYP2D6 substrate", Verdict::Warn,
            "Basic amine beta to an arene - CYP2D6-mediated metabolism / polymorphic clearance.");
    }
    if (f.catechol) {
        add("COMT substrate", Verdict::Warn, "Catechol - O-methylation by COMT.");
        add("MAO substrate", Verdict::Warn, "Catecholamine - oxidative deamination by MAO.");
        add("Catechol autoxidation", Verdict::Warn, "Forms reactive ortho-quinones.");
    }
    if (m.id == "acetaminophen") {
        add("Bioactivation -> NAPQI", Verdict::Danger,
            "CYP2E1/3A oxidation yields reactive NAPQI; glutathione-depleting hepatotoxin.");
    }
    if (m.id == "caffeine") {
        add("CYP1A2 substrate", Verdict::Info, "Primary metabolic route; smoking induces clearance.");
        add("Adenosine A1/A2A antagonism", Verdict::Info, "Mechanism of stimulant action.");
    }
    if (f.ester) {
        add("Esterase hydrolysis", Verdict::Info, "Carboxylesterase/plasma esterase clearance.");
    }
    // hERG-style cardiac liability screen.
    if (m.logP > 2.5 && m.molWeight > 250.0) {
        add("hERG liability (possible)", Verdict::Warn, "Lipophilic + sizeable - screen QT risk.");
    } else {
        add("hERG liability", Verdict::Info, "Low predicted hERG affinity.");
    }
    add("Plasma protein binding", Verdict::Info,
        m.logP > 2.0 ? "Moderate-high binding expected." : "Low-moderate binding expected.");

    Verdict overall = Verdict::Info;
    for (const auto& e : r.endpoints) overall = std::max(overall, e.verdict);
    r.overall = overall;
    r.summary = std::string("Overall ADMET: ") + verdictLabel(overall) + ". " +
                std::to_string(r.endpoints.size()) + " endpoints screened.";
    return r;
}

// ---------------------------------------------------------------- similarity
double pseudoSimilarity(const Molecule& a, const Molecule& b) {
    const double d =
        std::abs(a.logP - b.logP) * 0.6 +
        std::abs(a.tpsa - b.tpsa) * 0.02 +
        std::abs(a.molWeight - b.molWeight) * 0.006 +
        std::abs(a.hbd - b.hbd) * 0.3 +
        std::abs(a.hba - b.hba) * 0.2 +
        std::abs(a.rotatableBonds - b.rotatableBonds) * 0.15;
    return clampd(std::exp(-d * 0.55), 0.0, 1.0);
}

SimilarityReport computeSimilarity(const Molecule& m, const std::vector<Molecule>& lib) {
    SimilarityReport r;
    r.moleculeId = m.id;
    for (const auto& other : lib) {
        if (other.id == m.id) continue;
        const double t = pseudoSimilarity(m, other);
        const double p = clampd(t * 0.85 + 0.10, 0.0, 1.0);
        r.hits.push_back({other.name, other.drugClass, other.legalStatus, t, p});
    }
    std::sort(r.hits.begin(), r.hits.end(),
              [](const SimilarityHit& x, const SimilarityHit& y) { return x.tanimoto > y.tanimoto; });
    if (r.hits.size() > 8) r.hits.resize(8);
    if (!r.hits.empty()) {
        r.nearestName = r.hits.front().referenceName;
        r.nearestScore = r.hits.front().tanimoto;
    }
    r.summary = r.hits.empty() ? "No references." :
                ("Nearest: " + r.nearestName + " (Tanimoto " +
                 std::to_string(r.nearestScore).substr(0, 4) + ").");
    return r;
}

// --------------------------------------------------------------------- legal
LegalScorecard computeLegal(const Molecule& m, const std::vector<Molecule>& lib) {
    LegalScorecard r;
    r.moleculeId = m.id;
    r.jurisdiction = "US Federal Analogue Act (illustrative)";

    double bestControlled = 0.0;
    std::string nearestControlled;
    for (const auto& other : lib) {
        if (other.id == m.id) continue;
        const bool controlled = contains(other.legalStatus, "Schedule I") ||
                                contains(other.legalStatus, "Schedule II");
        if (!controlled) continue;
        const double t = pseudoSimilarity(m, other);
        if (t > bestControlled) { bestControlled = t; nearestControlled = other.name; }
    }
    r.substantialSimilarity = bestControlled * 100.0;

    if (r.substantialSimilarity >= 75)      r.classification = "Likely controlled analog";
    else if (r.substantialSimilarity >= 50) r.classification = "Possible analog (fact-specific)";
    else                                    r.classification = "Unlikely analog";

    r.rationale.push_back("Nearest controlled reference: " +
                          (nearestControlled.empty() ? "none" : nearestControlled));
    r.rationale.push_back("Structural-similarity component: " +
                          std::to_string(static_cast<int>(r.substantialSimilarity)) + "/100");
    r.rationale.push_back("NOTE: heuristic only - not legal advice; analog status is fact- and "
                          "intent-dependent.");
    return r;
}

// ------------------------------------------------------------------- docking
DockingResult computeDock(const Molecule& m, const std::string& target) {
    DockingResult r;
    r.moleculeId = m.id;
    r.targetName = target;
    const double base = -(5.0 + m.logP * 0.8 + m.molWeight * 0.004);
    for (int i = 0; i < 6; ++i) {
        Pose p;
        p.rank = i + 1;
        p.affinityKcalPerMol = base + static_cast<double>(i) * 0.35;
        p.rmsd = static_cast<double>(i) * 0.9;
        r.poses.push_back(p);
    }
    r.bestAffinity = r.poses.front().affinityKcalPerMol;
    r.summary = "Best predicted affinity " + std::to_string(r.bestAffinity).substr(0, 6) +
                " kcal/mol at " + target + " (binding/pharmacology only).";
    return r;
}

// ============================================================ fake module impls
class FakeLibrary final : public ILibrary {
public:
    explicit FakeLibrary(std::vector<Molecule> lib) : lib_(std::move(lib)) {}
    std::vector<Molecule> all() const override { return lib_; }
    std::optional<Molecule> byId(const std::string& id) const override {
        for (const auto& m : lib_) if (m.id == id) return m;
        return std::nullopt;
    }
    std::size_t count() const override { return lib_.size(); }
    const std::vector<Molecule>& ref() const { return lib_; }
private:
    std::vector<Molecule> lib_;
};

class FakeStability final : public IStabilityModule {
public:
    StabilityReport analyze(const Molecule& m) const override { return computeStability(m); }
};
class FakeAdmet final : public IAdmetModule {
public:
    AdmetReport screen(const Molecule& m) const override { return computeAdmet(m); }
};
class FakeAbsorption final : public IAbsorptionModule {
public:
    AbsorptionReport predict(const Molecule& m) const override { return computeAbsorption(m); }
};
class FakeSimilarity final : public ISimilarityModule {
public:
    explicit FakeSimilarity(const std::vector<Molecule>& lib) : lib_(lib) {}
    SimilarityReport search(const Molecule& m) const override { return computeSimilarity(m, lib_); }
private:
    const std::vector<Molecule>& lib_;
};
class FakeLegal final : public ILegalModule {
public:
    explicit FakeLegal(const std::vector<Molecule>& lib) : lib_(lib) {}
    LegalScorecard score(const Molecule& m) const override { return computeLegal(m, lib_); }
private:
    const std::vector<Molecule>& lib_;
};
class FakeDocking final : public IDockingModule {
public:
    std::vector<std::string> targets() const override {
        return {"DAT (dopamine transporter)", "NET (norepinephrine transporter)",
                "SERT (serotonin transporter)", "TAAR1", "hERG (cardiac safety)"};
    }
    DockingResult dock(const Molecule& m, const std::string& target) const override {
        return computeDock(m, target);
    }
};
class FakeRuns final : public IRunStore {
public:
    std::vector<RunRecord> recent() const override {
        return {
            {"run-0007", "Absorption", "Methamphetamine", "complete", "2026-06-08 02:40",
             "F ~90%, CNS-penetrant."},
            {"run-0006", "Stability", "Cocaine", "complete", "2026-06-08 02:31",
             "62/100 - ester hydrolysis limiting."},
            {"run-0005", "ADMET", "Amphetamine", "complete", "2026-06-08 02:20",
             "WARN - MAO + CYP2D6."},
            {"run-0004", "Docking", "MDMA -> SERT", "complete", "2026-06-08 02:05",
             "Best -7.8 kcal/mol."},
            {"run-0003", "Similarity", "Mephedrone", "complete", "2026-06-08 01:55",
             "Nearest: Methcathinone (0.86)."},
        };
    }
};

}  // namespace

// =============================================================== FakeBackend
struct FakeBackend::Impl {
    FakeLibrary    library{buildLibrary()};
    FakeStability  stability;
    FakeAdmet      admet;
    FakeAbsorption absorption;
    FakeSimilarity similarity{library.ref()};
    FakeLegal      legal{library.ref()};
    FakeDocking    docking;
    FakeRuns       runs;
};

FakeBackend::FakeBackend() : impl_(std::make_unique<Impl>()) {}
FakeBackend::~FakeBackend() = default;

Services FakeBackend::services() {
    Services s;
    s.library    = &impl_->library;
    s.stability  = &impl_->stability;
    s.admet      = &impl_->admet;
    s.absorption = &impl_->absorption;
    s.similarity = &impl_->similarity;
    s.legal      = &impl_->legal;
    s.docking    = &impl_->docking;
    s.runs       = &impl_->runs;
    return s;
}

}  // namespace biocad
