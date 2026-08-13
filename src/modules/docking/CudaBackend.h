// modules/docking/CudaBackend.h - first-party CUDA GPU docking backend.
//
// A genuine GPU docking engine for NVIDIA hardware, compiled ONLY with
// BIOCAD_ENABLE_CUDA (BIOCAD_HAVE_CUDA). It runs a RIGID-BODY pose search of the
// embedded ligand conformer inside the receptor box and scores every sampled pose on
// the GPU with the AutoDock Vina inter-molecular scoring function (Trott & Olson 2010).
//
// HONEST SCOPE: this is a rigid-body search - it does NOT vary ligand torsions, so
// AutoDock Vina (flexible, CPU) stays the accurate path; this engine is a fast GPU
// pre-screen / fallback, labeled as such. Like every IDockingBackend it degrades
// gracefully (returns real=false, never throws) when no CUDA device or no prepared
// receptor is present.
//
// SAFETY SCOPE: produces ligand->target BINDING AFFINITY (pharmacology / target
// engagement). No synthesis, route, precursor, or manufacturability content.
#pragma once

#include "contracts/IDockingBackend.h"

namespace biocad::docking {

class CudaBackend final : public IDockingBackend {
public:
    std::string   id() const override { return "cuda-gpu"; }
    std::string   displayName() const override { return "CUDA GPU (rigid, Vina scoring)"; }
    bool          available() const override;  // a usable CUDA device is present
    DockJobResult dock(const chem::Molecule& graph, const chem::Conformer& ligand3d,
                       const ReceptorTarget& target) const override;
};

}  // namespace biocad::docking
