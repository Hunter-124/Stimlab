// chem/Smiles.h - SMILES -> Molecule graph parser.
#pragma once

#include <optional>
#include <string_view>

#include "chem/Molecule.h"

namespace stimlab::chem {

// Parse a SMILES string into a molecular graph. Returns nullopt on malformed input.
// Supports: organic subset (B C N O P S F Cl Br I) + aromatic lowercase, bracket
// atoms with isotopes/charge/explicit-H/chirality, branches, ring-bond closures
// (single digit and %NN), and the bond symbols - = # : / \ and disconnection '.'.
std::optional<Molecule> parseSmiles(std::string_view smiles);

}  // namespace stimlab::chem
