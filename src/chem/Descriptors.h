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
// Ertl 2000 topological polar surface area (A^2), summed over N and O only.
//
// That is the convention every published TPSA threshold assumes (Veber <=140 A^2
// for oral bioavailability, the ~<=90 A^2 rule of thumb for CNS penetration), so
// it is the default. Ertl also tabulates sulfur and phosphorus; that sum is a
// DIFFERENT number, not a better one - the two conventions differ by 62 A^2 on
// famotidine, which straddles the Veber cutoff - so it gets its own function and
// any readout must state which one it used.
double tpsa(const Molecule& m);
double tpsaIncludingSulfurAndPhosphorus(const Molecule& m);
// Wildman-Crippen logP; delegates to chem::crippen() (chem/Crippen.h), which
// carries the citation and the method's ~0.67 log-unit RMS error. Returns 0.0
// if the descriptor pack is missing - call chem::crippen() to see that.
double crippenLogP(const Molecule& m);
int    formalCharge(const Molecule& m);
int    heavyAtomCount(const Molecule& m);

}  // namespace biocad::chem
