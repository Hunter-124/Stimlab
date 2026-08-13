// chem/Embed3D.h - generate 3D coordinates from the molecular graph (no RDKit).
//
// Pragmatic distance-geometry: build a distance matrix from bonded (1-2) and
// angle (1-3) ideal geometry (covalent radii + VSEPR angles), embed it with a
// metric-matrix (classical MDS) eigendecomposition via Eigen, then relax with a
// short steepest-descent on a simple bonded + angle + van-der-Waals force field.
// Explicit hydrogen positions are appended so the 3D viewer can draw them.
//
// This is the geometry seam shared by the 3D viewer (render/MolViewport) and the
// docking ligand-prep (a Conformer -> PDBQT). The API below is FROZEN.
#pragma once

#include <utility>
#include <vector>

#include "chem/Molecule.h"

namespace biocad::chem {

struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;
};

// A generated 3D conformation of a molecule. Heavy atoms occupy indices
// [0, heavyCount) in graph order (so pos[i] is the coordinate of Molecule::atoms[i]);
// appended explicit hydrogens occupy [heavyCount, size). `bonds` indexes into pos
// and includes the heavy-heavy graph bonds plus the generated heavy-H bonds.
struct Conformer {
    std::vector<Vec3> pos;                      // coordinates in Angstrom
    std::vector<int>  z;                        // atomic number per position
    std::vector<std::pair<int, int>> bonds;     // bonded index pairs (incl. C-H)
    int heavyCount = 0;                         // pos[0..heavyCount) are heavy atoms

    [[nodiscard]] bool empty() const { return pos.empty(); }
    [[nodiscard]] int  size() const { return static_cast<int>(pos.size()); }
};

// Element geometry helpers (Angstrom). Both fall back to a carbon-like default
// for elements outside the parsed set, never returning 0/NaN.
double covalentRadius(int z);   // single-bond covalent radius (bond-length geometry)
double vdwRadius(int z);        // van-der-Waals radius (clash term + spacefill display)

// Generate a 3D conformation for a molecular graph. Deterministic given `seed`
// (the seed only drives tiny symmetry-breaking jitter). Returns an empty
// Conformer for an empty molecule; never produces NaN/Inf coordinates.
Conformer embed3D(const Molecule& m, unsigned seed = 0x5715ABu);

}  // namespace biocad::chem
