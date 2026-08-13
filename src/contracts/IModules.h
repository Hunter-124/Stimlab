// contracts/IModules.h - frozen analysis-module interfaces.
// Every analysis capability is a pure-virtual interface; real impls (RDKit-backed)
// implement these. The UI codes ONLY against these contracts, never against an
// implementation - which is what lets a module be replaced without touching a panel.
//
// SAFETY SCOPE: interfaces describe identity, pharmacology (binding affinity),
// stability, absorption/PK, similarity and legal status. There is intentionally
// no synthesis/route/manufacturability interface.
#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "bio/Structure.h"
#include "contracts/IDockingBackend.h"
#include "data/Domain.h"

namespace biocad {

// Curated compound library (defaults + user imports).
class ILibrary {
public:
    virtual ~ILibrary() = default;
    virtual std::vector<Molecule> all() const = 0;
    virtual std::optional<Molecule> byId(const std::string& id) const = 0;
    virtual std::size_t count() const = 0;
};

// Molecular stability (replaces the out-of-scope manufacturability score).
class IStabilityModule {
public:
    virtual ~IStabilityModule() = default;
    virtual StabilityReport analyze(const Molecule& m) const = 0;
};

// ADMET / metabolism (D, M, E, T) - harmful metabolites, DDIs, safety flags.
class IAdmetModule {
public:
    virtual ~IAdmetModule() = default;
    virtual AdmetReport screen(const Molecule& m) const = 0;
};

// Structural alerts for metabolic bioactivation.
//
// SAFETY SCOPE: this interface emits LIABILITY FLAGS. A flag says a substructure
// has been associated with reactive-metabolite formation in the literature; it
// does not say the compound is toxic, and there is deliberately no toxicity
// verdict, no composite risk score and no Verdict::Danger anywhere in
// AlertReport. An empty flag list means "no alert in the pack matched", which is
// not a safety claim about the compound.
class IAlertsModule {
public:
    virtual ~IAlertsModule() = default;
    virtual AlertReport screen(const Molecule& m) const = 0;
};

// Curated, cited biotransformations for a library compound.
//
// SAFETY SCOPE: this interface RETRIEVES facts; it never enumerates hypotheses.
// Rule-based metabolite prediction was measured at 1.1-29% precision and
// 14.7-28.3% sensitivity across every major tool (Boyce et al. 2022, Comput
// Toxicol 21:100208), so hypothesis generation is a separate surface and there is
// no entry point here that would produce one. An empty MetabolismReport::known is
// a statement about BioCAD's curation, not about the compound, which is why
// coverageNote is a required field of the result.
class IMetabolismFactsModule {
public:
    virtual ~IMetabolismFactsModule() = default;
    virtual MetabolismReport known(const Molecule& m) const = 0;
};

// Absorption / pharmacokinetics (the "A" of ADMET) - permeability, F%, BBB.
class IAbsorptionModule {
public:
    virtual ~IAbsorptionModule() = default;
    virtual AbsorptionReport predict(const Molecule& m) const = 0;
};

// Structural + pharmacophore similarity vs the known-substance reference set.
class ISimilarityModule {
public:
    virtual ~ISimilarityModule() = default;
    virtual SimilarityReport search(const Molecule& m) const = 0;
};

// Legal-analog "substantially similar" scorecard.
class ILegalModule {
public:
    virtual ~ILegalModule() = default;
    virtual LegalScorecard score(const Molecule& m) const = 0;
};

// Docking = ligand->protein binding affinity (pharmacology/activity).
class IDockingModule {
public:
    virtual ~IDockingModule() = default;
    virtual std::vector<std::string> targets() const = 0;
    virtual DockingResult dock(const Molecule& m, const std::string& target) const = 0;

    // Phase D: richer result carrying 3D poses (for the molecular viewer) and the
    // engine/fallback provenance. Default = no 3D detail, so a legacy impl keeps
    // compiling; the real backend-backed impl overrides both.
    virtual DockJobResult dockDetailed(const Molecule& m, const std::string& target) const {
        (void)m; (void)target; return {};
    }
    // CNS receptor presets (PDB ref + binding-site box). Default = none.
    virtual std::vector<ReceptorTarget> presets() const { return {}; }
};

// Run history.
class IRunStore {
public:
    virtual ~IRunStore() = default;
    virtual std::vector<RunRecord> recent() const = 0;
    // Phase D: persist a run. Default no-op so a legacy store keeps compiling;
    // the SQLite-backed store overrides it.
    virtual void record(const RunRecord&) {}
};

// Pharmacodynamics and pharmacokinetics.
//
// SAFETY SCOPE: this interface emits EXPOSURE SCENARIOS under stated assumptions.
// It has no dose-recommendation entry point, and adding one is out of scope by
// design. Every returned number carries a Provenance; a missing prerequisite is a
// NotComputed Quantity naming what is missing, never a silently assumed default.
class IPharmacodynamicsModule {
public:
    virtual ~IPharmacodynamicsModule() = default;

    // Four-parameter logistic fit in log-concentration space.
    virtual CurveFit fitFourParameterLogistic(const std::vector<DoseResponsePoint>&) const = 0;

    // Cheng-Prusoff Ki from an IC50. Returns NotComputed naming the missing field
    // when the modality's required inputs are absent.
    virtual Quantity kiFromIc50(const ChengPrusoffInput&) const = 0;

    // Schild regression: pA2, slope with a confidence interval, and KB only when
    // the slope CI includes 1.
    virtual SchildResult schild(const std::vector<SchildPoint>&) const = 0;

    // Integrate an exposure profile for a dosing regimen.
    virtual PkProfile simulate(const PkModelSpec&, const DoseRegimen&) const = 0;

    // Fractional target occupancy from a free-concentration profile and a Kd.
    virtual OccupancyCurve occupancy(const PkProfile&, const Quantity& kd) const = 0;
};

// Pairwise sequence alignment (Gotoh affine gaps, BLOSUM62 from a pack).
//
// SAFETY SCOPE: alignGlobal() returns a NotComputed E-value. Karlin-Altschul
// statistics are defined for local alignments only; offering an E-value on a
// global alignment would attach a significance claim to a number that has none.
class ISequenceModule {
public:
    virtual ~ISequenceModule() = default;
    virtual SequenceAlignment alignGlobal(const std::string& a, const std::string& b) const = 0;
    virtual SequenceAlignment alignLocal(const std::string& a, const std::string& b) const = 0;
};

// Protein structure I/O and comparison.
//
// There is deliberately no entry point taking a chem::Conformer: that type is a
// small-molecule distance-geometry embedding, and comparing it as a protein
// would yield an lDDT/TM-score that looks meaningful and is not. bio::Structure
// has no such constructor, so the misuse is a compile error.
class IStructureModule {
public:
    virtual ~IStructureModule() = default;

    // Reads a local .pdb/.cif file (including one previously downloaded into the
    // cache). std::nullopt when the file is unreadable or the format is unknown;
    // recoverable problems arrive in Structure::warnings instead.
    virtual std::optional<bio::Structure> load(const std::filesystem::path& file) const = 0;

    virtual StructureComparison compare(const bio::Structure& ref,
                                        const bio::Structure& model) const = 0;

    // Solvent-accessible surface area. `source` carries the full parameter
    // string (algorithm, probe radius, point count, radii set, hydrogen policy)
    // because a SASA without them is not reproducible.
    virtual Quantity sasa(const bio::Structure& s) const = 0;
};

}  // namespace biocad
