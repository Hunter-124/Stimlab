// modules/docking/PdbqtWriter.h - minimal RIGID AutoDock PDBQT ligand writer.
//
// Serializes a chem::Conformer (+ its chem::Molecule graph, used for atom typing)
// into a single-ROOT, zero-torsion AutoDock PDBQT string. This deliberately emits
// a clean rigid ROOT/ENDROOT block (TORSDOF 0) rather than a Meeko-style flexible
// tree: a flexible PDBQT trips AutoDock Vina's tree.h parser (plan.md sec 9), while
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

#include "chem/Embed3D.h"
#include "chem/Molecule.h"

namespace stimlab::docking {

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
    std::string text;       // full PDBQT (ROOT ... ENDROOT ... TORSDOF 0)
    int         atomCount;  // number of ATOM lines written
    int         heavyCount; // heavy-atom subset of atomCount
    int         polarH;     // polar-hydrogen subset of atomCount
};

// Write a rigid PDBQT for `conf` (typed via `graph`). `resName` is the 3-letter
// residue label (default "LIG"). Coordinates are written as-is (Angstrom). Returns
// an empty .text with atomCount==0 for an empty conformer.
PdbqtLigand writeRigidPdbqt(const chem::Molecule& graph, const chem::Conformer& conf,
                            const std::string& resName = "LIG");

}  // namespace stimlab::docking
