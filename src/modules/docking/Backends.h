// modules/docking/Backends.h - IDockingBackend implementations (Phase D, WP-3).
//
// EstimateBackend  - always-available descriptor fallback (real=false). Produces
//                    finite, ranked affinities from logP/MW and attaches the
//                    embedded ligand conformer to every pose so the 3D viewer has
//                    geometry even with no engine present.
// VinaBackend      - AutoDock Vina / smina subprocess (real=true when it runs):
//                    rigid-PDBQT ligand prep -> subprocess -> parse REMARK VINA
//                    RESULT poses. Degrades to a real=false result (never throws)
//                    when its binary or a prepared receptor is missing, so the
//                    caller can fall back to the estimate.
//
// SAFETY SCOPE: these produce ligand->target BINDING AFFINITY (pharmacology). No
// synthesis/route/manufacturability content is produced anywhere.
#pragma once

#include "contracts/IDockingBackend.h"
#include "modules/docking/EngineLocator.h"

namespace stimlab::docking {

class EstimateBackend final : public IDockingBackend {
public:
    std::string id() const override { return "estimate"; }
    std::string displayName() const override { return "descriptor-estimate"; }
    bool        available() const override { return true; }
    DockJobResult dock(const chem::Molecule& graph, const chem::Conformer& ligand3d,
                       const ReceptorTarget& target) const override;
};

class VinaBackend final : public IDockingBackend {
public:
    explicit VinaBackend(Engine engine) : engine_(engine) {}
    std::string id() const override;
    std::string displayName() const override;
    bool        available() const override;   // binary located (receptor checked in dock)
    DockJobResult dock(const chem::Molecule& graph, const chem::Conformer& ligand3d,
                       const ReceptorTarget& target) const override;

private:
    Engine engine_;
};

}  // namespace stimlab::docking
