// contracts/IModules.h - frozen analysis-module interfaces.
// Every analysis capability is a pure-virtual interface; real impls (RDKit-backed)
// and thick fakes both implement these. The UI codes ONLY against these contracts.
//
// SAFETY SCOPE: interfaces describe identity, pharmacology (binding affinity),
// stability, absorption/PK, similarity and legal status. There is intentionally
// no synthesis/route/manufacturability interface.
#pragma once

#include <optional>
#include <string>
#include <vector>

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
    // engine/fallback provenance. Default = no 3D detail, so legacy/fake impls keep
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
    // Phase D: persist a run. Default no-op so fakes/legacy stores keep compiling;
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

}  // namespace biocad
