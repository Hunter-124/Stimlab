#include "modules/RealBackend.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <string>
#include <vector>

#include "chem/AdmetModel.h"
#include "chem/Analysis.h"
#include "chem/Descriptors.h"
#include "chem/Embed3D.h"
#include "chem/Smiles.h"
#include "modules/BioModules.h"
#include "modules/docking/Backends.h"
#include "modules/docking/Presets.h"
#include "modules/docking/ReceptorPrep.h"
#include "packs/Pack.h"
#include "pkpd/Pharmacodynamics.h"
#ifdef BIOCAD_HAVE_CUDA
#  include "modules/docking/CudaBackend.h"
#endif
#include "storage/RunStore.h"

namespace biocad {
namespace {

namespace chem = biocad::chem;

double clampd(double v, double lo, double hi) { return std::max(lo, std::min(hi, v)); }

// Library compounds come from the loaded data packs (assets/packs/*.json plus the
// %APPDATA%/BioCAD/packs overlay). There is no hard-coded catalog: the domain is
// whatever the packs say it is.

// The real engine always recomputes descriptors from the SMILES: any property a
// pack authored is metadata for the pack browser, never a substitute for running
// the in-house chem engine.
Molecule computeMolecule(const packs::PackCompound& c) {
    Molecule m = c.molecule();
    if (auto pm = chem::parseSmiles(c.smiles)) {
        m.formula = chem::molecularFormula(*pm);
        m.molWeight = chem::molecularWeight(*pm);
        m.logP = chem::crippenLogP(*pm);
        m.tpsa = chem::tpsa(*pm);
        m.hbd = chem::hbdCount(*pm);
        m.hba = chem::hbaCount(*pm);
        m.rotatableBonds = chem::rotatableBondCount(*pm);
    }
    return m;
}

chem::FunctionalGroups groupsOf(const Molecule& m) {
    if (auto pm = chem::parseSmiles(m.smiles)) return chem::detectGroups(*pm);
    return {};
}

bool hasAnilide(const Molecule& m) {
    auto pm = chem::parseSmiles(m.smiles);
    if (!pm) return false;
    for (size_t i = 0; i < pm->atoms.size(); ++i) {
        const auto& a = pm->atoms[i];
        if (a.z != 7) continue;
        bool amide = false, aryl = false;
        for (int nb : a.nbr) {
            const auto& x = pm->atoms[nb];
            if (x.z == 6 && x.aromatic) aryl = true;
            if (x.z == 6) {
                for (int b2 : x.bonds)
                    if (pm->bonds[b2].order == 2.0 && pm->atoms[pm->bonds[b2].other(nb)].z == 8)
                        amide = true;
            }
        }
        if (amide && aryl) return true;
    }
    return false;
}

// A flexible-chain primary/secondary amine whose alpha carbon is an unsubstituted
// methylene (-CH2-) is readily deaminated by monoamine oxidase, so it carries a
// heavy presystemic-metabolism burden. The strict tests below pick out exactly the
// MAO-vulnerable beta-phenethylamine class while excluding the look-alikes:
//   * alpha-METHYL amines (amphetamine-type): the alpha carbon is a CH (1 H), not
//     a CH2, so it is MAO-resistant and not flagged.
//   * N-METHYL substituents (e.g. methamphetamine's N-CH3): a terminal CH3 has 3 H,
//     not 2, so the methyl neighbour does not count as the alpha methylene.
//   * RING amines (piperidine/pyrrolidine/tropane in methylphenidate, MDPV, cocaine,
//     nicotine): the nitrogen and its alpha carbons are in a ring and are excluded.
// Combined with the phenethylamine flag at the call site so only Ar-C-C-N amines
// reach here.
bool hasMaoLabileAmine(const Molecule& m) {
    auto pm = chem::parseSmiles(m.smiles);
    if (!pm) return false;
    for (const auto& a : pm->atoms) {
        if (a.z != 7 || a.aromatic || a.inRing) continue;
        bool amideN = false;  // skip nitrogens bonded to a carbonyl carbon
        for (int nb : a.nbr) {
            const auto& x = pm->atoms[nb];
            if (x.z != 6) continue;
            for (int b2 : x.bonds)
                if (pm->bonds[b2].order == 2.0 && pm->atoms[pm->bonds[b2].other(nb)].z == 8) amideN = true;
        }
        if (amideN) continue;
        const int h = a.totalH(), deg = a.degree();
        const bool primaryOrSecondary = (deg == 1 && h >= 2) || (deg == 2 && h >= 1);
        if (!primaryOrSecondary) continue;
        for (int nb : a.nbr) {
            const auto& c = pm->atoms[nb];
            // alpha carbon must be an acyclic methylene (-CH2-), i.e. exactly 2 H
            if (c.z == 6 && !c.aromatic && !c.inRing && c.totalH() == 2) return true;
        }
    }
    return false;
}

// ----------------------------------------------------------------- stability
StabilityReport realStability(const Molecule& m) {
    const auto f = groupsOf(m);
    StabilityReport r;
    r.moleculeId = m.id;

    const double hydrolysis = f.ester ? 34.0 : (f.amide ? 70.0 : 92.0);
    const double oxidation = f.catechol ? 22.0 : (f.phenol ? 58.0 : (m.logP > 2.0 ? 76.0 : 85.0));
    const double photolysis = (f.catechol || f.arylKetone) ? 54.0 : 82.0;
    const double thermal = f.arylKetone ? 66.0 : (f.ester ? 72.0 : 88.0);
    const double phStab = (f.ester || f.arylKetone || f.catechol) ? 58.0 : 80.0;

    // Predicted, numeric stability ranges (not just a "narrow/broad" word): the
    // degradation chemistry of the perceived groups sets an actual pH window and
    // temperature ceiling, which we surface in the rationale text below.
    chem::PkLiabilities sl;
    sl.catechol = f.catechol;
    sl.phenol = f.phenol;
    sl.ester = f.ester;
    sl.amide = f.amide;
    sl.arylKetone = f.arylKetone;
    sl.methylenedioxy = f.methylenedioxy;
    const auto ph = chem::predictPhWindow(sl);
    const auto th = chem::predictThermalWindow(sl);

    r.factors = {
        {"Hydrolysis resistance", hydrolysis,
         f.ester ? "Ester detected - susceptible to acid/base/enzymatic hydrolysis."
                 : (f.amide ? "Amide present - slow hydrolysis under stress."
                            : "No hydrolyzable ester; resistant.")},
        {"Oxidation resistance", oxidation,
         f.catechol ? "Catechol - rapid autoxidation to ortho-quinone."
                    : (f.phenol ? "Phenol - moderate oxidative liability." : "Few oxidation-prone motifs.")},
        {"Photostability", photolysis,
         (f.catechol || f.arylKetone) ? "Chromophore/carbonyl - photodegradation risk."
                                      : "No strong photolabile chromophore."},
        {"Thermal stability (range)", thermal, chem::thermalWindowText(th)},
        {"pH stability (range)", phStab, chem::phWindowText(ph)},
    };
    double sum = 0;
    for (const auto& x : r.factors) sum += x.score;
    r.overallScore = sum / static_cast<double>(r.factors.size());

    if (f.ester) r.degradants.push_back({"Hydrolyzed acid + alcohol", "Ester hydrolysis",
                                         "Primary route under humidity."});
    if (f.amide) r.degradants.push_back({"Amine + carboxylic acid", "Amide hydrolysis",
                                         "Slow; extreme pH accelerated."});
    if (f.catechol) r.degradants.push_back({"ortho-Quinone / melanin polymers", "Autoxidation",
                                            "Air/light accelerated."});
    if (f.arylKetone) r.degradants.push_back({"Aldol / reduced products", "Carbonyl reactivity",
                                              "Condensation under heat."});
    if (r.degradants.empty()) r.degradants.push_back({"N-oxide (minor)", "N-oxidation", "Trace."});

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
                (f.ester ? "ester hydrolysis." : (f.catechol ? "catechol oxidation."
                          : (f.arylKetone ? "carbonyl reactivity." : "none major."))) +
                " Predicted stable pH " + chem::fmt1(ph.low) + "-" + chem::fmt1(ph.high) +
                ", to ~" + chem::fmt0(th.stableToC) + "C (store " + store + ").";
    return r;
}

// ---------------------------------------------------------------- absorption
AbsorptionReport realAbsorption(const Molecule& m) {
    AbsorptionReport r;
    r.moleculeId = m.id;
    const double tpsa = m.tpsa, logP = m.logP, mw = m.molWeight;

    // Oral F is computed mechanistically as F = absorbed fraction x first-pass
    // survival, NOT pinned to a fixed 95% ceiling. The presystemic metabolic
    // liabilities perceived from the structure drive the first-pass extraction,
    // so e.g. catecholamines (catechol -> COMT/MAO) correctly read near-zero F
    // while metabolically robust stimulants read high.
    const auto g = groupsOf(m);
    chem::PkLiabilities L;
    L.catechol = g.catechol;
    L.phenol = g.phenol;
    L.ester = g.ester;
    L.amide = g.amide;
    L.arylKetone = g.arylKetone;
    L.methylenedioxy = g.methylenedioxy;
    L.maoLabileAmine = g.phenethylamine && hasMaoLabileAmine(m);
    const auto bio = chem::predictBioavailability(mw, logP, tpsa, m.hbd, L);

    r.hiaPct = bio.hiaPct;
    r.bioavailabilityPct = bio.bioavailabilityPct;
    r.caco2LogPapp = clampd(-4.80 + 0.30 * logP - 0.008 * tpsa, -7.0, -4.0);
    r.logBB = clampd(0.15 + 0.17 * logP - 0.011 * tpsa, -2.0, 1.2);
    r.logS = clampd(0.80 - 0.0100 * mw - 0.55 * logP, -6.5, 1.0);
    r.pgpSubstrate = (mw > 400.0 && m.hba >= 6) || tpsa > 120.0;
    r.cnsPenetrant = (r.logBB > -0.30) && (tpsa < 90.0) && (mw < 450.0);

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
         "Higher (less negative) = more permeable."},
        {"BBB partition", r.logBB, "logBB", r.cnsPenetrant ? Verdict::Good : Verdict::Info,
         r.cnsPenetrant ? "Crosses the blood-brain barrier; central site of action plausible." : "Limited blood-brain-barrier partition."},
        {"Aqueous solubility", r.logS, "logS", band(r.logS, -4.0, -5.0),
         "Low solubility can dissolution-limit absorption."},
        {"P-gp efflux substrate", r.pgpSubstrate ? 1.0 : 0.0, "bool",
         r.pgpSubstrate ? Verdict::Warn : Verdict::Good,
         r.pgpSubstrate ? "Efflux may reduce net absorption." : "Not a likely P-gp substrate."},
    };
    r.summary = "Hepatic availability under the stated assumptions ~" +
                chem::fmt0(r.bioavailabilityPct) + "% (HIA ~" + chem::fmt0(r.hiaPct) +
                "%, F_H ~" + chem::fmt0(bio.firstPassSurvival * 100.0) +
                "% from an ASSUMED fu.CLint; " + bio.limitingRoute + "). " +
                (r.cnsPenetrant ? "BBB-permeant." : "Low BBB partition.");
    return r;
}

// --------------------------------------------------------------------- admet
AdmetReport realAdmet(const Molecule& m) {
    const auto f = groupsOf(m);
    AdmetReport r;
    r.moleculeId = m.id;
    auto add = [&](std::string name, Verdict v, std::string detail) {
        r.endpoints.push_back({std::move(name), v, std::move(detail)});
    };

    if (f.phenethylamine && !f.catechol) {
        add("MAO-A/B substrate", Verdict::Warn,
            "Phenethylamine core - monoamine oxidase metabolism; tyramine/MAOI DDI risk.");
        add("CYP2D6 substrate", Verdict::Warn,
            "Basic amine beta to an arene - CYP2D6 metabolism / polymorphic clearance.");
    }
    if (f.catechol) {
        add("COMT substrate", Verdict::Warn, "Catechol - O-methylation by COMT.");
        add("MAO substrate", Verdict::Warn, "Catecholamine - oxidative deamination by MAO.");
        add("Catechol autoxidation", Verdict::Warn, "Forms reactive ortho-quinones.");
    }
    if (f.phenol && hasAnilide(m)) {
        add("Bioactivation -> quinone-imine (NAPQI-type)", Verdict::Danger,
            "Para-aminophenol/anilide - CYP oxidation to a reactive electrophile; glutathione-depleting.");
    }
    if (f.ester) add("Esterase hydrolysis", Verdict::Info, "Carboxylesterase/plasma esterase clearance.");
    if (f.nitrile) add("Nitrile metabolism", Verdict::Warn, "Possible cyanide release on oxidation.");
    if (m.logP > 2.5 && m.molWeight > 250.0)
        add("hERG liability (possible)", Verdict::Warn, "Lipophilic + sizeable - screen QT risk.");
    else
        add("hERG liability", Verdict::Info, "Low predicted hERG affinity.");
    add("Plasma protein binding", Verdict::Info,
        m.logP > 2.0 ? "Moderate-high binding expected." : "Low-moderate binding expected.");

    Verdict overall = Verdict::Info;
    for (const auto& e : r.endpoints) overall = std::max(overall, e.verdict);
    r.overall = overall;
    r.summary = std::string("Overall ADMET: ") + verdictLabel(overall) + ". " +
                std::to_string(r.endpoints.size()) + " endpoints screened (structure-derived).";
    return r;
}

// ============================================================ module classes
class RealLibrary final : public ILibrary {
public:
    explicit RealLibrary(std::vector<Molecule> lib) : lib_(std::move(lib)) {}
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

class RealStability final : public IStabilityModule {
public:
    StabilityReport analyze(const Molecule& m) const override { return realStability(m); }
};
class RealAdmet final : public IAdmetModule {
public:
    AdmetReport screen(const Molecule& m) const override { return realAdmet(m); }
};
class RealAbsorption final : public IAbsorptionModule {
public:
    AbsorptionReport predict(const Molecule& m) const override { return realAbsorption(m); }
};

class RealSimilarity final : public ISimilarityModule {
public:
    explicit RealSimilarity(const std::vector<Molecule>& lib) : lib_(lib) {
        for (const auto& m : lib_) {
            std::vector<std::uint32_t> fp, fp1;
            if (auto pm = chem::parseSmiles(m.smiles)) {
                fp = chem::morganFingerprint(*pm, 2);
                fp1 = chem::morganFingerprint(*pm, 1);
            }
            fps_.push_back(std::move(fp));
            fps1_.push_back(std::move(fp1));
        }
    }
    SimilarityReport search(const Molecule& m) const override {
        SimilarityReport r;
        r.moleculeId = m.id;
        std::vector<std::uint32_t> q, q1;
        if (auto pm = chem::parseSmiles(m.smiles)) {
            q = chem::morganFingerprint(*pm, 2);
            q1 = chem::morganFingerprint(*pm, 1);
        }
        for (std::size_t i = 0; i < lib_.size(); ++i) {
            if (lib_[i].id == m.id) continue;
            const double t = chem::tanimoto(q, fps_[i]);
            const double p = chem::tanimoto(q1, fps1_[i]);
            r.hits.push_back({lib_[i].name, lib_[i].drugClass, lib_[i].legalStatus, t, p});
        }
        std::sort(r.hits.begin(), r.hits.end(),
                  [](const SimilarityHit& a, const SimilarityHit& b) { return a.tanimoto > b.tanimoto; });
        if (r.hits.size() > 8) r.hits.resize(8);
        if (!r.hits.empty()) { r.nearestName = r.hits.front().referenceName; r.nearestScore = r.hits.front().tanimoto; }
        r.summary = r.hits.empty() ? "No references." :
                    ("Nearest: " + r.nearestName + " (Tanimoto " + std::to_string(r.nearestScore).substr(0, 4) +
                     ", ECFP4).");
        return r;
    }
private:
    const std::vector<Molecule>& lib_;
    std::vector<std::vector<std::uint32_t>> fps_, fps1_;
};

class RealLegal final : public ILegalModule {
public:
    RealLegal(const std::vector<Molecule>& lib, const RealSimilarity& sim) : lib_(lib), sim_(sim) {}
    LegalScorecard score(const Molecule& m) const override {
        LegalScorecard r;
        r.moleculeId = m.id;
        r.jurisdiction = "US Federal Analogue Act (illustrative)";
        const auto sim = sim_.search(m);
        double best = 0.0; std::string nearest;
        for (const auto& h : sim.hits) {
            const bool controlled = h.legalStatus.find("Schedule I") != std::string::npos ||
                                    h.legalStatus.find("Schedule II") != std::string::npos;
            if (controlled && h.tanimoto > best) { best = h.tanimoto; nearest = h.referenceName; }
        }
        r.substantialSimilarity = best * 100.0;
        if (r.substantialSimilarity >= 75) r.classification = "Likely controlled analog";
        else if (r.substantialSimilarity >= 50) r.classification = "Possible analog (fact-specific)";
        else r.classification = "Unlikely analog";
        r.rationale.push_back("Nearest controlled reference: " + (nearest.empty() ? "none" : nearest));
        r.rationale.push_back("Structural-similarity (ECFP4 Tanimoto) component: " +
                              std::to_string(static_cast<int>(r.substantialSimilarity)) + "/100");
        r.rationale.push_back("NOTE: heuristic only - not legal advice; analog status is fact-dependent.");
        return r;
    }
private:
    const std::vector<Molecule>& lib_;
    const RealSimilarity& sim_;
};

// Docking now runs behind the frozen IDockingBackend seam: pick the best AVAILABLE
// real engine (Vina, then smina) for the resolved CNS receptor preset; if none can
// actually dock (no binary / no prepared receptor), fall back to the clearly-labeled
// descriptor estimate. Both the legacy DockingResult and the 3D-bearing DockJobResult
// are derived from one path so the two UI views always agree.
class RealDocking final : public IDockingModule {
public:
    std::vector<std::string> targets() const override { return docking::presetNames(); }

    std::vector<ReceptorTarget> presets() const override { return docking::targetPresets(); }

    DockJobResult dockDetailed(const Molecule& m, const std::string& target) const override {
        ReceptorTarget tgt;
        if (const ReceptorTarget* p = docking::findPreset(target)) tgt = *p;
        else { tgt.id = target; tgt.name = target; }

        // Bind a prepared receptor (cache-only; no network on the hot path). When one
        // exists, an available engine can produce a REAL dock; otherwise the backends
        // report "no prepared receptor" and we fall back to the labeled estimate. The
        // prepared receptor also carries a box center in ITS OWN coordinate frame
        // (from the co-crystal ligand); use it so the search box actually overlaps the
        // receptor - the preset's literature-style center may be in a different frame.
        if (tgt.receptorPath.empty()) {
            const auto rec = docking::locatePreparedReceptor(tgt.id);
            if (rec.ready) {
                tgt.receptorPath = rec.path;
                if (rec.hasBox) { tgt.box.cx = rec.cx; tgt.box.cy = rec.cy; tgt.box.cz = rec.cz; }
            }
        }

        auto graph = chem::parseSmiles(m.smiles);
        if (!graph) {
            DockJobResult r;
            r.engine = "none"; r.provenance = Provenance::Heuristic; r.targetId = tgt.id;
            r.log = "Ligand SMILES could not be parsed.";
            return r;
        }
        const chem::Conformer conf = chem::embed3D(*graph);

        for (const IDockingBackend* b : realEngines()) {
            if (!b->available()) continue;
            DockJobResult r = b->dock(*graph, conf, tgt);
            if (r.fromEngine() && !r.poses.empty()) { r.targetId = tgt.id; return r; }
        }
        DockJobResult r = estimate_.dock(*graph, conf, tgt);
        r.targetId = tgt.id;
        return r;
    }

    DockingResult dock(const Molecule& m, const std::string& target) const override {
        const DockJobResult d = dockDetailed(m, target);
        DockingResult r;
        r.moleculeId = m.id;
        r.targetName = target;
        for (const auto& p : d.poses) r.poses.push_back({p.rank, p.affinityKcalPerMol, p.rmsdLb});
        r.bestAffinity = d.bestAffinity();
        const std::string aff = std::to_string(r.bestAffinity).substr(0, 6);
        r.summary = d.fromEngine()
            ? ("Best affinity " + aff + " kcal/mol at " + target + " (docked with " + d.engine + ").")
            : ("Estimated affinity " + aff + " kcal/mol at " + target +
               " (" + d.engine + " - structure-descriptor model, not a docked score).");
        return r;
    }

private:
    using BackendList = std::vector<const IDockingBackend*>;
    // Engines to TRY, in order, for the current compute mode. The backends are static
    // singletons (stateless); only the small ordered pointer list is built per call
    // (mode is a runtime user setting, so this can't be one fixed static list).
    //   Auto - Vina/smina (the accurate CPU path) first; then the GPU engines as a
    //          fallback when no CPU engine binary is present. Among GPU engines, Vina-GPU
    //          (the real Vina Monte-Carlo + BFGS search, OpenCL) before the first-party
    //          CUDA engine (a fixed rotation x translation grid).
    //   Cpu  - Vina/smina only (skip the GPU engines).
    //   Gpu  - GPU engines only: Vina-GPU first, then CUDA.
    // Vina-GPU is a subprocess engine so it is compiled in EVERY build (no CUDA toolkit
    // needed); the first-party CUDA engine only exists under BIOCAD_HAVE_CUDA.
    // dockDetailed() falls back to the labeled descriptor estimate when none of these
    // produce a real dock, so every mode degrades honestly.
    BackendList realEngines() const {
        static const docking::VinaBackend vina{docking::Engine::Vina};
        static const docking::VinaBackend smina{docking::Engine::Smina};
        static const docking::VinaGpuBackend vinaGpu;
#ifdef BIOCAD_HAVE_CUDA
        static const docking::CudaBackend cuda;
#endif
        const docking::ComputeMode mode = docking::computeMode();
        BackendList list;
        if (mode == docking::ComputeMode::Gpu) {
            list.push_back(&vinaGpu);  // GPU forced: Vina-GPU (real Vina search) first...
#ifdef BIOCAD_HAVE_CUDA
            list.push_back(&cuda);     // ...then the first-party CUDA grid engine
#endif
        } else if (mode == docking::ComputeMode::Cpu) {
            list.push_back(&vina);
            list.push_back(&smina);  // CPU only: skip the GPU engines
        } else {                     // Auto: accurate CPU first, GPU engines as a fallback
            list.push_back(&vina);
            list.push_back(&smina);
            list.push_back(&vinaGpu);
#ifdef BIOCAD_HAVE_CUDA
            list.push_back(&cuda);
#endif
        }
        return list;
    }
    docking::EstimateBackend estimate_;
};

}  // namespace

struct RealBackend::Impl {
    // The loaded packs are the catalog. A pack that failed to load is reported
    // here, surfaced by the Presets panel, and never silently dropped.
    packs::LoadReport packReport = packs::loadBuiltin();
    RealLibrary library{[this] {
        std::vector<Molecule> v;
        for (const auto& c : packReport.compounds()) v.push_back(computeMolecule(c));
        return v;
    }()};
    RealStability stability;
    RealAdmet admet;
    RealAbsorption absorption;
    RealSimilarity similarity{library.ref()};
    RealLegal legal{library.ref(), similarity};
    RealDocking docking;
    SqliteRunStore runs;
    pkpd::RealPharmacodynamics pharmacodynamics;
    RealSequence sequence;
    RealStructure structure;
};

RealBackend::RealBackend() : impl_(std::make_unique<Impl>()) {}
RealBackend::~RealBackend() = default;

Services RealBackend::services() {
    Services s;
    s.library = &impl_->library;
    s.stability = &impl_->stability;
    s.admet = &impl_->admet;
    s.absorption = &impl_->absorption;
    s.similarity = &impl_->similarity;
    s.legal = &impl_->legal;
    s.docking = &impl_->docking;
    s.runs = &impl_->runs;
    s.pharmacodynamics = &impl_->pharmacodynamics;
    s.sequence = &impl_->sequence;
    s.structure = &impl_->structure;
    return s;
}

}  // namespace biocad
