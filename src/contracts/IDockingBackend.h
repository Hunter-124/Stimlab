// contracts/IDockingBackend.h - FROZEN docking-engine seam (Phase D).
//
// An IDockingBackend turns a 3D ligand + a receptor box into ranked, 3D-bearing
// binding poses. Real engines (AutoDock Vina, smina) and a labeled descriptor
// estimate all implement it; the UI's docking module (IDockingModule) selects a
// backend and exposes the result. The 3D ligand geometry comes from the shared
// chem::Embed3D conformer seam (WP-1).
//
// SAFETY SCOPE: a DockPose carries ligand->protein BINDING AFFINITY, a
// pharmacology / target-engagement signal. It is NEVER a synthesis or
// make-it signal. There is no synthesis/route/manufacturability data here.
#pragma once

#include <string>
#include <vector>

#include "chem/Embed3D.h"
#include "chem/Molecule.h"

namespace biocad {

// Search box for docking, in Angstroms (Vina center + size convention).
struct DockBox {
    double cx = 0.0, cy = 0.0, cz = 0.0;     // center
    double sx = 22.0, sy = 22.0, sz = 22.0;  // size (edge length)
};

// A CNS receptor preset: identity + PDB reference + the binding-site box derived
// from the co-crystal ligand. An optional prepared-receptor path (PDBQT) is set
// once a receptor has been provisioned/prepared; empty means "box only".
struct ReceptorTarget {
    std::string id;             // stable key, e.g. "DAT"
    std::string name;           // display, e.g. "DAT (dopamine transporter)"
    std::string pdb;            // PDB id, e.g. "4M48"
    std::string receptorPath;   // prepared receptor file (PDBQT); may be empty
    DockBox     box;
};

// One scored binding pose. `ligand` holds the docked ligand coordinates (heavy +
// H, indexed like a chem::Conformer); it may be empty if the backend produced a
// score without 3D coordinates (e.g. the descriptor estimate).
struct DockPose {
    int    rank = 0;
    double affinityKcalPerMol = 0.0;   // more negative = stronger predicted binding
    double rmsdLb = 0.0;               // RMSD lower bound vs the top pose (Vina)
    double rmsdUb = 0.0;               // RMSD upper bound vs the top pose (Vina)
    chem::Conformer ligand;            // docked 3D ligand (may be empty)
};

// The outcome of one docking job. `real` distinguishes a true engine run from the
// labeled descriptor fallback so the UI can be honest about what produced a score.
struct DockJobResult {
    std::string engine;          // "AutoDock Vina 1.2.5", "smina", "descriptor-estimate"
    bool        real = false;    // true iff a real docking engine produced the poses
    std::string targetId;        // ReceptorTarget::id used
    std::vector<DockPose> poses; // ranked best (most negative) first
    std::string log;             // status / notes / fallback reason

    // ---- search confidence (set by the convergent real-engine search) ---------
    // A real dock now repeats an independent flexible-ligand search with fresh seeds
    // and escalating exhaustiveness until the best affinity stabilises (or a budget
    // is hit). These let the UI report how trustworthy the number is.
    int    searchRuns = 1;        // independent search runs aggregated into this result
    double affinitySpread = 0.0;  // kcal/mol spread of the per-run best affinity
    bool   converged = false;     // best affinity stabilised within tolerance
    int    torsions = 0;          // active rotatable bonds the ligand was docked with

    [[nodiscard]] bool empty() const { return poses.empty(); }
    [[nodiscard]] double bestAffinity() const {
        return poses.empty() ? 0.0 : poses.front().affinityKcalPerMol;
    }
};

// A swappable docking engine. Implementations provision/locate their binary and
// run a single ligand->target job; they must degrade gracefully (return a result
// with real=false rather than throwing) when their engine is unavailable.
class IDockingBackend {
public:
    virtual ~IDockingBackend() = default;

    virtual std::string id() const = 0;            // "vina" / "smina" / "estimate"
    virtual std::string displayName() const = 0;   // human label
    virtual bool        available() const = 0;     // engine binary present?

    virtual DockJobResult dock(const chem::Molecule& ligandGraph,
                               const chem::Conformer& ligand3d,
                               const ReceptorTarget& target) const = 0;
};

}  // namespace biocad
