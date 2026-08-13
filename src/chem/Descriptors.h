// chem/Descriptors.h - molecular descriptors computed from the graph.
// Real, literature-based methods (Lipinski/Veber counts, Ertl TPSA,
// Crippen-style logP). Validated against known values in tests/test_chem.cpp.
#pragma once

#include <string>

#include "chem/Molecule.h"

namespace biocad::chem {

std::string molecularFormula(const Molecule& m);   // Hill system, includes H
double molecularWeight(const Molecule& m);          // average mass, g/mol
int    hbdCount(const Molecule& m);                 // H-bond donors (N/O bearing H)
int    hbaCount(const Molecule& m);                 // H-bond acceptors (Lipinski N+O)
int    rotatableBondCount(const Molecule& m);       // Veber definition
int    ringCount(const Molecule& m);                // SSSR size (cyclomatic number)
int    aromaticAtomCount(const Molecule& m);
double fractionCsp3(const Molecule& m);             // 0..1
double tpsa(const Molecule& m);                     // Ertl 2000 topological PSA (A^2)
double crippenLogP(const Molecule& m);              // Wildman-Crippen-style estimate
int    formalCharge(const Molecule& m);
int    heavyAtomCount(const Molecule& m);

}  // namespace biocad::chem
