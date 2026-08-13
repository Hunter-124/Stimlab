// modules/docking/ReceptorPrep.h - prepare a rigid receptor PDBQT (WP-F).
//
// Turns a Protein Data Bank structure into a Vina-ready RIGID receptor PDBQT, with
// NO RDKit and NO hard OpenBabel linkage. Two paths, both producing the same
// cached artifact under %APPDATA%/BioCAD/runtime/receptors/<id>.pdbqt:
//
//   1. Built-in writer (always available): pdbToRigidReceptor() parses PDB
//      ATOM/HETATM records, strips waters + crystallization additives + the
//      co-crystal ligand, infers elements, assigns AutoDock4 atom types, writes a
//      heavy-atom rigid receptor (no torsion tree - receptors are rigid by
//      definition). Vina's empirical scoring is partial-charge-INDEPENDENT, so a
//      correctly-typed heavy-atom receptor docks; we still emit approximate
//      charges for AutoDock4-style consumers. Fully unit-testable from a string.
//   2. obabel.exe (preferred when provisioned): higher-quality prep with added
//      polar H (obabel <pdb> -xr -O <pdbqt>); used iff an obabel binary is located.
//
// ensureReceptor() is the orchestrator: returns a cached prepared receptor if one
// exists; otherwise (allowDownload) best-effort fetches the preset's PDB from RCSB,
// prepares it, caches it, and returns the path. It NEVER throws and NEVER blocks the
// caller on a hard failure - it degrades to "not ready" with a human note so the
// docking module can fall back to the clearly-labeled descriptor estimate.
//
// SAFETY SCOPE: receptor prep defines WHERE/AGAINST WHAT a ligand binds for a
// binding-affinity (target-engagement / pharmacology) prediction. There is no
// synthesis, route, precursor, or manufacturability content anywhere here.
#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "contracts/IDockingBackend.h"

namespace biocad::docking {

// The directory prepared receptor PDBQTs (and their fetched source PDBs) live in:
// %APPDATA%/BioCAD/runtime/receptors.
std::filesystem::path receptorsDir();

// Result of the built-in PDB -> rigid receptor PDBQT conversion.
struct ReceptorPdbqt {
    std::string text;            // full receptor PDBQT (ATOM records, no torsion tree)
    int         atomCount = 0;   // receptor ATOM lines written (heavy atoms)
    int         keptHetero = 0;  // non-water HETATM atoms kept (cofactors in keepHetero)
    int         keptModified = 0;   // modified-residue HETATM atoms kept as protein (e.g. MSE)
    int         droppedHetero = 0;  // waters/additives/ligand atoms removed
    std::string note;            // human-readable summary of what was kept/dropped

    // Docking box center derived from the structure's OWN coordinate frame (the
    // co-crystal ligand centroid when present, else the receptor centroid). This is
    // what makes a fetched receptor actually overlap the search box - the preset's
    // literature-style center may be in a different frame. boxSource explains which.
    double cx = 0.0, cy = 0.0, cz = 0.0;
    bool   hasBox = false;
    std::string boxSource;       // "co-crystal ligand <RES>" | "receptor centroid"

    [[nodiscard]] bool empty() const { return atomCount == 0; }
};

// Derive a docking box CENTER from a PDB structure's own coordinates: the centroid
// of the largest non-water HETATM residue (the co-crystal ligand) if one with >= 8
// atoms exists, else the centroid of all protein atoms. Returns false only for an
// atom-less structure. `source` is filled with a human-readable origin.
bool receptorBoxFromPdb(const std::string& pdbText, double& cx, double& cy, double& cz,
                        std::string& source);

// Convert PDB text to a rigid receptor PDBQT (built-in path; no external tools).
// Keeps standard protein ATOM records (first model, altLoc A/blank). Drops all
// HETATM by default; pass residue names in `keepHetero` (e.g. {"FAD"}) to retain a
// biologically relevant cofactor. Deterministic; never produces NaN coordinates.
ReceptorPdbqt pdbToRigidReceptor(const std::string& pdbText,
                                 const std::vector<std::string>& keepHetero = {});

// Outcome of locating / preparing a receptor for a target.
struct ReceptorPrepResult {
    bool        ready = false;   // a usable prepared PDBQT exists at `path`
    std::string path;            // runtime/receptors/<id>.pdbqt (when ready)
    std::string note;            // status / failure reason
    bool        usedObabel = false;
    int         atomCount = 0;

    // Structure-frame box center recovered from the prepared receptor's BIOCAD_BOX
    // remark (see receptorBoxFromPdb). When hasBox is true the docking module
    // overrides the preset center with this so the box overlaps the real receptor.
    double cx = 0.0, cy = 0.0, cz = 0.0;
    bool   hasBox = false;
};

// Ensure a prepared receptor PDBQT exists for `target`. If a cached one is present,
// returns it immediately (no network). Otherwise, iff `allowDownload`, fetches the
// preset's PDB from RCSB, prepares it (obabel if located, else the built-in writer),
// caches both, and returns the path. Best-effort; never throws.
ReceptorPrepResult ensureReceptor(const ReceptorTarget& target, bool allowDownload = false);

// Cache-only lookup: is a prepared receptor already on disk for this target id?
// (No network, no preparation - used on the hot docking path.)
ReceptorPrepResult locatePreparedReceptor(const std::string& targetId);

}  // namespace biocad::docking
