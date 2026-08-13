// modules/docking/PdbqtWriter.h - minimal RIGID AutoDock PDBQT ligand writer.
//
// Serializes a chem::Conformer (+ its chem::Molecule graph, used for atom typing)
// into a single-ROOT, zero-torsion AutoDock PDBQT string. This deliberately emits
// a clean rigid ROOT/ENDROOT block (TORSDOF 0) rather than a Meeko-style flexible
// tree: a flexible PDBQT trips AutoDock Vina's tree.h parser, while
// a rigid ROOT always docks. When a real OpenBabel (obabel.exe) is provisioned the
// caller may instead shell out for a torsionally-flexible ligand; this writer is
// the always-available fallback and is fully unit-testable from its string output.
//
// Atom types follow the AutoDock4 set (C, A=aromatic C, N, NA, OA, SA, HD=polar H,
// F, Cl, Br, I, P, S, ...). Partial charges are an APPROXIMATE per-element /
// electronegativity scheme (documented below), not a full Gasteiger iteration;
// Vina is tolerant of approximate charges for ranking. Only polar hydrogens
// (bonded to N/O/S) are emitted, matching AutoDock's united-atom convention.
//
// SAFETY SCOPE: this is ligand geometry/typing for a binding-affinity prediction.
// There is no synthesis/route content anywhere in this file.
#pragma once

#include <string>
#include <vector>

#include "chem/Embed3D.h"
#include "chem/Molecule.h"

namespace biocad::docking {

// AutoDock atom type for a single position, given its element and bonding context.
// `aromatic` selects A (aromatic carbon) / NA / SA variants; `polarH` marks an HD.
std::string autodockAtomType(int z, bool aromatic, bool polarH);

// Approximate partial charge (electron units) for one atom from its element and an
// electronegativity-difference smear over its bonds. NOT a full Gasteiger solve.
double approxPartialCharge(const chem::Molecule& graph, const chem::Conformer& conf,
                           int posIndex);

// Result of writing a ligand: the PDBQT text plus the count of ATOM records emitted
// (heavy atoms + polar hydrogens), so callers/tests can verify without re-parsing.
struct PdbqtLigand {
    std::string text;       // full PDBQT (ROOT ... ENDROOT ... TORSDOF n)
    int         atomCount;  // number of ATOM lines written
    int         heavyCount; // heavy-atom subset of atomCount
    int         polarH;     // polar-hydrogen subset of atomCount
    int         torsions = 0;  // active rotatable bonds (TORSDOF; 0 for the rigid writer)
    // For each emitted ATOM (in file order), the source conformer index. The engine
    // preserves atom order in its output, so this maps a docked pose's coordinates
    // back onto the original conformer topology (the flexible writer reorders atoms
    // into the torsion tree, so an index-only mapping would scramble elements/bonds).
    std::vector<int> serialToConf;
};

// Write a rigid PDBQT for `conf` (typed via `graph`). `resName` is the 3-letter
// residue label (default "LIG"). Coordinates are written as-is (Angstrom). Returns
// an empty .text with atomCount==0 for an empty conformer.
PdbqtLigand writeRigidPdbqt(const chem::Molecule& graph, const chem::Conformer& conf,
                            const std::string& resName = "LIG");

// Write a TORSIONALLY FLEXIBLE PDBQT so the docking engine can BEND and ROTATE the
// ligand. Rotatable bonds (acyclic single bonds between two non-terminal heavy atoms,
// excluding amide C-N) are detected from `graph`; the ligand is split into rigid
// fragments on them and serialised as an AutoDock ROOT/BRANCH/ENDBRANCH torsion tree
// with TORSDOF = #rotatable bonds. Atom serials are assigned in strict file order so
// every BRANCH references already-emitted atoms (the correctness the engine's tree
// parser requires). Degrades to a single rigid ROOT when no rotatable bond exists.
PdbqtLigand writeFlexiblePdbqt(const chem::Molecule& graph, const chem::Conformer& conf,
                               const std::string& resName = "LIG");

}  // namespace biocad::docking
