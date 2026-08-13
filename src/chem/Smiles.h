// chem/Smiles.h - SMILES -> Molecule graph parser.
#pragma once

#include <optional>
#include <string_view>

#include "chem/Molecule.h"

namespace biocad::chem {

// Parse a SMILES string into a molecular graph. Returns nullopt on malformed input.
// Supports: organic subset (B C N O P S F Cl Br I) + aromatic lowercase, bracket
// atoms with isotopes/charge/explicit-H/chirality, branches, ring-bond closures
// (single digit and %NN), and the bond symbols - = # : / \ and disconnection '.'.
std::optional<Molecule> parseSmiles(std::string_view smiles);

// Fill Atom::implicitH for every non-bracket atom from the standard organic-subset
// valences, given the bonds already present.
//
// This is public because THREE places need the identical rule, and a second copy of
// it is a round-trip bug waiting to happen: the parser fills it while reading, the
// canonical writer must predict exactly what the parser will re-infer or the round
// trip silently loses hydrogens, and a graph built by hand (the structure sketcher)
// has no parser run to fill it in at all. One function, one rule.
void assignImplicitHydrogens(Molecule& m);

}  // namespace biocad::chem
