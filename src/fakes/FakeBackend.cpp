#include "fakes/FakeBackend.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <string_view>
#include <vector>

#include "chem/AdmetModel.h"
#include "packs/Pack.h"
#include "pkpd/Pharmacodynamics.h"

namespace biocad {
namespace {

bool contains(std::string_view hay, std::string_view needle) {
    return hay.find(needle) != std::string_view::npos;
}

double clampd(double v, double lo, double hi) { return std::max(lo, std::min(hi, v)); }

// ------------------------------------------------------------------ defaults
// The fake's catalog is a data pack like any other - the same schema, parsed by
// the same packs::parse() - just embedded in the binary so tests stay hermetic and
// deterministic and never touch the filesystem. There is no second catalog format.
constexpr const char* kFixturePack = R"PACK({
  "schemaVersion": 1,
  "id": "fake-reference",
  "title": "Deterministic test fixture pack",
  "description": "Hermetic in-binary catalog used by FakeBackend so tests never touch the filesystem. Same schema, same parser, same code path as a shipped pack.",
  "compounds": [
    {
      "id": "caffeine",
      "name": "Caffeine",
      "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
      "drugClass": "Xanthine stimulant",
      "legalUs": "Unscheduled (US)",
      "notes": "Adenosine receptor antagonist; common reference stimulant.",
      "properties": {
        "formula": "C8H10N4O2",
        "molWeight": 194.19,
        "logP": -0.07,
        "tpsa": 58.4,
        "hbd": 0,
        "hba": 6,
        "rotatableBonds": 0
      }
    },
    {
      "id": "amphetamine",
      "name": "Amphetamine",
      "smiles": "CC(N)Cc1ccccc1",
      "drugClass": "Phenethylamine stimulant",
      "legalUs": "Schedule II (US)",
      "notes": "Releaser/reuptake inhibitor at DAT/NET.",
      "properties": {
        "formula": "C9H13N",
        "molWeight": 135.21,
        "logP": 1.76,
        "tpsa": 26.0,
        "hbd": 1,
        "hba": 1,
        "rotatableBonds": 2
      }
    },
    {
      "id": "methamphetamine",
      "name": "Methamphetamine",
      "smiles": "CNC(C)Cc1ccccc1",
      "drugClass": "Phenethylamine stimulant",
      "legalUs": "Schedule II (US)",
      "notes": "N-methyl amphetamine; higher CNS penetration.",
      "properties": {
        "formula": "C10H15N",
        "molWeight": 149.23,
        "logP": 2.07,
        "tpsa": 12.0,
        "hbd": 1,
        "hba": 1,
        "rotatableBonds": 2
      }
    },
    {
      "id": "mdma",
      "name": "MDMA",
      "smiles": "CNC(C)Cc1ccc2OCOc2c1",
      "drugClass": "Entactogen (phenethylamine)",
      "legalUs": "Schedule I (US)",
      "notes": "3,4-methylenedioxy substitution; serotonergic releaser.",
      "properties": {
        "formula": "C11H15NO2",
        "molWeight": 193.24,
        "logP": 1.06,
        "tpsa": 30.0,
        "hbd": 1,
        "hba": 3,
        "rotatableBonds": 3
      }
    },
    {
      "id": "methylphenidate",
      "name": "Methylphenidate",
      "smiles": "COC(=O)C(c1ccccc1)C1CCCCN1",
      "drugClass": "Piperidine stimulant",
      "legalUs": "Schedule II (US)",
      "notes": "DAT/NET reuptake inhibitor; methyl ester (hydrolysis-labile).",
      "properties": {
        "formula": "C14H19NO2",
        "molWeight": 233.31,
        "logP": 2.3,
        "tpsa": 38.3,
        "hbd": 1,
        "hba": 3,
        "rotatableBonds": 4
      }
    },
    {
      "id": "modafinil",
      "name": "Modafinil",
      "smiles": "NC(=O)CS(=O)C(c1ccccc1)c1ccccc1",
      "drugClass": "Eugeroic (wakefulness)",
      "legalUs": "Schedule IV (US)",
      "notes": "Atypical DAT inhibitor; sulfinyl + primary amide.",
      "properties": {
        "formula": "C15H15NO2S",
        "molWeight": 273.35,
        "logP": 1.4,
        "tpsa": 79.1,
        "hbd": 1,
        "hba": 3,
        "rotatableBonds": 4
      }
    },
    {
      "id": "cocaine",
      "name": "Cocaine",
      "smiles": "COC(=O)C1C(OC(=O)c2ccccc2)CC3CCC1N3C",
      "drugClass": "Tropane stimulant",
      "legalUs": "Schedule II (US)",
      "notes": "Two ester groups -> hydrolysis to benzoylecgonine/ecgonine.",
      "properties": {
        "formula": "C17H21NO4",
        "molWeight": 303.35,
        "logP": 2.3,
        "tpsa": 55.8,
        "hbd": 0,
        "hba": 5,
        "rotatableBonds": 5
      }
    },
    {
      "id": "pseudoephedrine",
      "name": "Pseudoephedrine",
      "smiles": "CNC(C)C(O)c1ccccc1",
      "drugClass": "Phenethylamine decongestant",
      "legalUs": "OTC/Schedule V (US)",
      "notes": "Beta-hydroxy amphetamine analog; sympathomimetic.",
      "properties": {
        "formula": "C10H15NO",
        "molWeight": 165.23,
        "logP": 0.9,
        "tpsa": 32.3,
        "hbd": 2,
        "hba": 2,
        "rotatableBonds": 3
      }
    },
    {
      "id": "ephedrine",
      "name": "Ephedrine",
      "smiles": "CNC(C)C(O)c1ccccc1",
      "drugClass": "Phenethylamine sympathomimetic",
      "legalUs": "Schedule (varies)",
      "notes": "Diastereomer of pseudoephedrine.",
      "properties": {
        "formula": "C10H15NO",
        "molWeight": 165.23,
        "logP": 1.0,
        "tpsa": 32.3,
        "hbd": 2,
        "hba": 2,
        "rotatableBonds": 3
      }
    },
    {
      "id": "cathinone",
      "name": "Cathinone",
      "smiles": "CC(N)C(=O)c1ccccc1",
      "drugClass": "Beta-keto phenethylamine",
      "legalUs": "Schedule I (US)",
      "notes": "Beta-keto group reduces stability (oxidation/condensation).",
      "properties": {
        "formula": "C9H11NO",
        "molWeight": 149.19,
        "logP": 0.6,
        "tpsa": 43.1,
        "hbd": 1,
        "hba": 2,
        "rotatableBonds": 2
      }
    },
    {
      "id": "methcathinone",
      "name": "Methcathinone",
      "smiles": "CNC(C)C(=O)c1ccccc1",
      "drugClass": "Cathinone stimulant",
      "legalUs": "Schedule I (US)",
      "notes": "N-methyl cathinone.",
      "properties": {
        "formula": "C10H13NO",
        "molWeight": 163.22,
        "logP": 0.9,
        "tpsa": 29.1,
        "hbd": 1,
        "hba": 2,
        "rotatableBonds": 2
      }
    },
    {
      "id": "mephedrone",
      "name": "Mephedrone (4-MMC)",
      "smiles": "CNC(C)C(=O)c1ccc(C)cc1",
      "drugClass": "Substituted cathinone",
      "legalUs": "Schedule I (US)",
      "notes": "4-methyl methcathinone; 'bath salts' class.",
      "properties": {
        "formula": "C11H15NO",
        "molWeight": 177.24,
        "logP": 1.6,
        "tpsa": 29.1,
        "hbd": 1,
        "hba": 2,
        "rotatableBonds": 3
      }
    },
    {
      "id": "mda",
      "name": "MDA",
      "smiles": "CC(N)Cc1ccc2OCOc2c1",
      "drugClass": "Entactogen (phenethylamine)",
      "legalUs": "Schedule I (US)",
      "notes": "3,4-methylenedioxyamphetamine; MDMA's primary metabolite.",
      "properties": {
        "formula": "C10H13NO2",
        "molWeight": 179.22,
        "logP": 1.0,
        "tpsa": 47.6,
        "hbd": 1,
        "hba": 3,
        "rotatableBonds": 2
      }
    },
    {
      "id": "methylone",
      "name": "Methylone (bk-MDMA)",
      "smiles": "CNC(C)C(=O)c1ccc2OCOc2c1",
      "drugClass": "Substituted cathinone",
      "legalUs": "Schedule I (US)",
      "notes": "Beta-keto MDMA analog.",
      "properties": {
        "formula": "C11H13NO3",
        "molWeight": 207.23,
        "logP": 0.7,
        "tpsa": 56.4,
        "hbd": 1,
        "hba": 4,
        "rotatableBonds": 3
      }
    },
    {
      "id": "mdpv",
      "name": "MDPV",
      "smiles": "O=C(C(CCC)N1CCCC1)c1ccc2OCOc2c1",
      "drugClass": "Substituted cathinone",
      "legalUs": "Schedule I (US)",
      "notes": "Potent DAT/NET inhibitor; pyrrolidinophenone.",
      "properties": {
        "formula": "C16H21NO3",
        "molWeight": 275.34,
        "logP": 2.8,
        "tpsa": 47.6,
        "hbd": 0,
        "hba": 4,
        "rotatableBonds": 4
      }
    },
    {
      "id": "alpha-pvp",
      "name": "alpha-PVP",
      "smiles": "O=C(C(CCC)N1CCCC1)c1ccccc1",
      "drugClass": "Pyrrolidinophenone stimulant",
      "legalUs": "Schedule I (US)",
      "notes": "'Flakka'; MDPV des-methylenedioxy analog.",
      "properties": {
        "formula": "C15H21NO",
        "molWeight": 231.33,
        "logP": 3.2,
        "tpsa": 20.3,
        "hbd": 0,
        "hba": 2,
        "rotatableBonds": 4
      }
    },
    {
      "id": "4-fa",
      "name": "4-Fluoroamphetamine",
      "smiles": "CC(N)Cc1ccc(F)cc1",
      "drugClass": "Phenethylamine stimulant",
      "legalUs": "Schedule (varies)",
      "notes": "Para-fluoro amphetamine; releaser.",
      "properties": {
        "formula": "C9H12FN",
        "molWeight": 153.2,
        "logP": 1.8,
        "tpsa": 26.0,
        "hbd": 1,
        "hba": 2,
        "rotatableBonds": 2
      }
    },
    {
      "id": "dmaa",
      "name": "Methylhexanamine (DMAA)",
      "smiles": "CCCC(C)CC(C)N",
      "drugClass": "Aliphatic amine stimulant",
      "legalUs": "Banned supplement (US)",
      "notes": "Non-aromatic sympathomimetic; no arene for CYP2D6/MAO motifs.",
      "properties": {
        "formula": "C7H17N",
        "molWeight": 115.22,
        "logP": 2.0,
        "tpsa": 26.0,
        "hbd": 1,
        "hba": 1,
        "rotatableBonds": 4
      }
    },
    {
      "id": "bupropion",
      "name": "Bupropion",
      "smiles": "CC(NC(C)(C)C)C(=O)c1cccc(Cl)c1",
      "drugClass": "Aminoketone (NDRI)",
      "legalUs": "Rx (US)",
      "notes": "Antidepressant/cessation aid; cathinone-like aminoketone.",
      "properties": {
        "formula": "C13H18ClNO",
        "molWeight": 239.74,
        "logP": 3.2,
        "tpsa": 29.1,
        "hbd": 1,
        "hba": 2,
        "rotatableBonds": 3
      }
    },
    {
      "id": "phenethylamine",
      "name": "Phenethylamine (parent)",
      "smiles": "NCCc1ccccc1",
      "drugClass": "Phenethylamine (trace amine)",
      "legalUs": "Endogenous",
      "notes": "Structural parent of the amphetamine class; TAAR1 agonist.",
      "properties": {
        "formula": "C8H11N",
        "molWeight": 121.18,
        "logP": 1.4,
        "tpsa": 26.0,
        "hbd": 1,
        "hba": 1,
        "rotatableBonds": 2
      }
    },
    {
      "id": "tyramine",
      "name": "Tyramine",
      "smiles": "NCCc1ccc(O)cc1",
      "drugClass": "Phenethylamine (trace amine)",
      "legalUs": "Endogenous",
      "notes": "Dietary; MAO substrate (MAOI 'cheese' interaction).",
      "properties": {
        "formula": "C8H11NO",
        "molWeight": 137.18,
        "logP": 0.9,
        "tpsa": 46.2,
        "hbd": 2,
        "hba": 2,
        "rotatableBonds": 2
      }
    },
    {
      "id": "theobromine",
      "name": "Theobromine",
      "smiles": "Cn1cnc2c1c(=O)[nH]c(=O)n2C",
      "drugClass": "Xanthine stimulant",
      "legalUs": "Unscheduled (US)",
      "notes": "Cocoa alkaloid; milder than caffeine.",
      "properties": {
        "formula": "C7H8N4O2",
        "molWeight": 180.16,
        "logP": -0.8,
        "tpsa": 67.2,
        "hbd": 1,
        "hba": 6,
        "rotatableBonds": 0
      }
    },
    {
      "id": "theophylline",
      "name": "Theophylline",
      "smiles": "Cn1c(=O)c2[nH]cnc2n(C)c1=O",
      "drugClass": "Xanthine bronchodilator",
      "legalUs": "Rx (US)",
      "notes": "Narrow therapeutic index; CYP1A2 substrate.",
      "properties": {
        "formula": "C7H8N4O2",
        "molWeight": 180.16,
        "logP": -0.02,
        "tpsa": 69.3,
        "hbd": 1,
        "hba": 6,
        "rotatableBonds": 0
      }
    },
    {
      "id": "norepinephrine",
      "name": "Norepinephrine",
      "smiles": "NCC(O)c1ccc(O)c(O)c1",
      "drugClass": "Catecholamine neurotransmitter",
      "legalUs": "Endogenous",
      "notes": "Catechol; COMT/MAO substrate; poor oral/CNS absorption.",
      "properties": {
        "formula": "C8H11NO3",
        "molWeight": 169.18,
        "logP": -1.2,
        "tpsa": 86.7,
        "hbd": 4,
        "hba": 4,
        "rotatableBonds": 2
      }
    },
    {
      "id": "epinephrine",
      "name": "Epinephrine",
      "smiles": "CNCC(O)c1ccc(O)c(O)c1",
      "drugClass": "Catecholamine hormone",
      "legalUs": "Rx (US)",
      "notes": "Adrenaline; catechol; not orally bioavailable.",
      "properties": {
        "formula": "C9H13NO3",
        "molWeight": 183.2,
        "logP": -1.4,
        "tpsa": 75.7,
        "hbd": 3,
        "hba": 4,
        "rotatableBonds": 2
      }
    },
    {
      "id": "lisdexamfetamine",
      "name": "Lisdexamfetamine",
      "smiles": "CC(Cc1ccccc1)NC(=O)C(N)CCCCN",
      "drugClass": "Amphetamine prodrug",
      "legalUs": "Schedule II (US)",
      "notes": "L-lysine amide prodrug of dexamphetamine; amide hydrolysis activates.",
      "properties": {
        "formula": "C15H25N3O",
        "molWeight": 263.38,
        "logP": 0.8,
        "tpsa": 95.0,
        "hbd": 3,
        "hba": 4,
        "rotatableBonds": 8
      }
    },
    {
      "id": "atomoxetine",
      "name": "Atomoxetine (reference)",
      "smiles": "CNCCC(Oc1ccccc1C)c1ccccc1",
      "drugClass": "Selective NRI (non-stimulant)",
      "legalUs": "Rx (US)",
      "notes": "Non-controlled ADHD reference; NET-selective.",
      "properties": {
        "formula": "C17H21NO",
        "molWeight": 255.35,
        "logP": 3.9,
        "tpsa": 21.3,
        "hbd": 1,
        "hba": 2,
        "rotatableBonds": 5
      }
    },
    {
      "id": "dopamine",
      "name": "Dopamine",
      "smiles": "NCCc1ccc(O)c(O)c1",
      "drugClass": "Catecholamine neurotransmitter",
      "legalUs": "Endogenous",
      "notes": "Catechol -> rapid autoxidation; COMT/MAO substrate.",
      "properties": {
        "formula": "C8H11NO2",
        "molWeight": 153.18,
        "logP": -0.98,
        "tpsa": 66.5,
        "hbd": 3,
        "hba": 3,
        "rotatableBonds": 2
      }
    },
    {
      "id": "nicotine",
      "name": "Nicotine",
      "smiles": "CN1CCCC1c1cccnc1",
      "drugClass": "Pyridine alkaloid stimulant",
      "legalUs": "Unscheduled (US)",
      "notes": "nAChR agonist.",
      "properties": {
        "formula": "C10H14N2",
        "molWeight": 162.23,
        "logP": 1.17,
        "tpsa": 16.1,
        "hbd": 0,
        "hba": 2,
        "rotatableBonds": 1
      }
    },
    {
      "id": "phentermine",
      "name": "Phentermine",
      "smiles": "CC(C)(N)Cc1ccccc1",
      "drugClass": "Phenethylamine anorectic",
      "legalUs": "Schedule IV (US)",
      "notes": "Alpha,alpha-dimethyl amphetamine; resists MAO.",
      "properties": {
        "formula": "C10H15N",
        "molWeight": 149.23,
        "logP": 1.9,
        "tpsa": 26.0,
        "hbd": 1,
        "hba": 1,
        "rotatableBonds": 1
      }
    },
    {
      "id": "acetaminophen",
      "name": "Acetaminophen (reference)",
      "smiles": "CC(=O)Nc1ccc(O)cc1",
      "drugClass": "Analgesic (non-stimulant reference)",
      "legalUs": "OTC (US)",
      "notes": "Metabolism reference: bioactivation to reactive NAPQI.",
      "properties": {
        "formula": "C8H9NO2",
        "molWeight": 151.16,
        "logP": 0.46,
        "tpsa": 49.3,
        "hbd": 2,
        "hba": 2,
        "rotatableBonds": 1
      }
    }
  ],
  "targets": [
    {
      "id": "DAT",
      "name": "DAT (dopamine transporter)",
      "pdb": "4M48",
      "box": {
        "cx": -2.1,
        "cy": 12.4,
        "cz": -6.8,
        "sx": 20.0,
        "sy": 20.0,
        "sz": 20.0
      },
      "headline": true
    },
    {
      "id": "NET",
      "name": "NET (norepinephrine transporter)",
      "pdb": "5I6X",
      "box": {
        "cx": 8.3,
        "cy": -4.2,
        "cz": 14.1,
        "sx": 20.0,
        "sy": 20.0,
        "sz": 20.0
      },
      "headline": true
    },
    {
      "id": "SERT",
      "name": "SERT (serotonin transporter)",
      "pdb": "5I73",
      "box": {
        "cx": -7.4,
        "cy": 16.2,
        "cz": 9.5,
        "sx": 22.0,
        "sy": 22.0,
        "sz": 22.0
      },
      "headline": true
    }
  ]
})PACK";

std::vector<Molecule> buildLibrary() {
    static const packs::Pack pack = packs::parseString(kFixturePack, "<embedded fixture pack>");
    std::vector<Molecule> v;
    v.reserve(pack.compounds.size());
    for (const auto& c : pack.compounds) v.push_back(c.molecule());
    return v;
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
         r.cnsPenetrant ? "Crosses the blood-brain barrier; relevant for a central site of action."
                        : "Limited blood-brain-barrier partition."},
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
                (r.cnsPenetrant ? "BBB-permeant." : "Low BBB partition.");
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
             "F ~90%, BBB-permeant."},
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

// The PK/PD fake is NOT a second set of numbers: the real implementation is a set
// of pure, deterministic functions (no clock, no randomness, no filesystem), so the
// fake delegates to it and keeps one code path under test.
class FakePharmacodynamics final : public IPharmacodynamicsModule {
public:
    CurveFit fitFourParameterLogistic(
        const std::vector<DoseResponsePoint>& pts) const override {
        return real_.fitFourParameterLogistic(pts);
    }
    Quantity kiFromIc50(const ChengPrusoffInput& in) const override {
        return real_.kiFromIc50(in);
    }
    SchildResult schild(const std::vector<SchildPoint>& pts) const override {
        return real_.schild(pts);
    }
    PkProfile simulate(const PkModelSpec& spec, const DoseRegimen& reg) const override {
        return real_.simulate(spec, reg);
    }
    OccupancyCurve occupancy(const PkProfile& p, const Quantity& kd) const override {
        return real_.occupancy(p, kd);
    }
private:
    pkpd::RealPharmacodynamics real_;
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
    FakePharmacodynamics pharmacodynamics;
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
    s.pharmacodynamics = &impl_->pharmacodynamics;
    return s;
}

}  // namespace biocad
