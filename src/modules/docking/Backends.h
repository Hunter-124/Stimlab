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
// VinaGpuBackend  - Vina-GPU (OpenCL) docking on the GPU (real=true when it runs).
//                    Mirrors VinaBackend exactly - same rigid-PDBQT ligand prep, same
//                    REMARK VINA RESULT pose parsing - but shells the GPU executable in
//                    its own folder (so it finds its compiled Kernel2_Opt.bin) and uses
//                    Vina-GPU's --thread lanes instead of Vina's --exhaustiveness. It
//                    runs the REAL AutoDock Vina search (Monte-Carlo + BFGS local
//                    optimisation) on the GPU, so it is a more thorough search than the
//                    first-party CUDA engine's fixed rotation x translation grid.
//                    Degrades to real=false (never throws) when the engine, its compiled
//                    kernel, or a prepared receptor is missing.
//
// SAFETY SCOPE: these produce ligand->target BINDING AFFINITY (pharmacology). No
// synthesis/route/manufacturability content is produced anywhere.
#pragma once

#include "contracts/IDockingBackend.h"
#include "modules/docking/EngineLocator.h"

namespace biocad::docking {

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

class VinaGpuBackend final : public IDockingBackend {
public:
    std::string id() const override { return "vina-gpu"; }
    std::string displayName() const override { return "GPU (OpenCL, Vina-GPU)"; }
    bool        available() const override;   // exe + a compiled kernel located
    DockJobResult dock(const chem::Molecule& graph, const chem::Conformer& ligand3d,
                       const ReceptorTarget& target) const override;
};

}  // namespace biocad::docking
