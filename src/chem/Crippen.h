// chem/Crippen.h - Wildman & Crippen (1999) atomic-contribution logP and molar
// refractivity, computed faithfully from the published class table.
//
// METHOD. S. A. Wildman and G. M. Crippen, "Prediction of Physicochemical
// Parameters by Atomic Contributions", J Chem Inf Comput Sci 1999;39(5):868-873
// (doi:10.1021/ci990307l). Every atom of the molecule - hydrogens included - is
// assigned to exactly one of the paper's classes (C1-C27, H1-H4, N1-N14, O1-O12,
// S1-S3, P, F, Cl, Br, I, Hal, Me1, Me2, plus the per-element catch-alls), and
//     logP = sum of the class logP contributions
//     MR   = sum of the class MR contributions
// Nothing is fitted here and nothing is invented: the classes, their SMARTS
// encodings, and both contribution columns are DATA, in
// assets/packs/descriptors/crippen.json, and the classification is first match
// wins in the array order that file ships. That is exactly how the paper's
// decision table works, so re-ordering the array changes the answer.
//
// HONESTY. The method's own published performance is an RMS error of about 0.67
// log units against experimental logP over its 9920-compound training set, with
// individual outliers far worse. Being faithful to the method is the goal of this
// file; agreeing with experiment is NOT claimed. Any Quantity built from this
// number must carry crippenCitation() as its source so the reader can see both
// the method and its error bar.
//
// HYDROGENS. chem::Molecule stores no explicit hydrogen atoms - H counts live on
// the heavy atom as Atom::totalH(). So the four hydrogen classes cannot be
// matched by SMARTS here and are instead derived from the attached heavy atom,
// which is the same decision the paper's own table encodes (its H patterns look
// only at the heavy neighbour, and at most one atom beyond it):
//     H1  0.1230  hydrocarbon H          H on carbon
//     H2 -0.2677  alcohol H              H on O whose neighbour is sp3 C or
//                                        aromatic c, or is neither C/N/O/S;
//                                        also water, and H on any heavy atom
//                                        that is not C, N or O (e.g. S-H, P-H)
//     H3  0.2142  amine H                H on nitrogen, and H on an O-N oxygen
//     H4  0.2980  acid H                 H on O whose carbon neighbour is doubly
//                                        bonded to C, N, O or S (carboxylic acid,
//                                        enol, oxime-type), and H on an O-O or
//                                        O-S oxygen
// The order above is the reference table's order and is likewise first match
// wins: a phenol OH is H2, not H4, because H2's "aromatic c" rule is tested
// first, and a carboxylic acid OH is H4 because its carbon is CX3, not CX4.
#pragma once

#include <string>
#include <vector>

#include "chem/Molecule.h"

namespace biocad::chem {

struct CrippenResult {
    double                   logP = 0.0;
    double                   molarRefractivity = 0.0;  // cm^3/mol
    std::vector<std::string> atomTypes;     // one class per heavy atom, for audit
    std::vector<std::string> unclassified;  // "<index>:<Z>" per unclassified heavy atom
    bool                     ok = false;    // false => logP/MR are meaningless, do not show them
    std::string              note;          // why ok is false, or the pack that was loaded
};

// Perceives rings and aromaticity on its own copy of `m` (the input is never
// mutated), so Kekule and aromatic SMILES of the same molecule give the same
// answer. A missing or malformed descriptor pack yields ok == false with `note`
// naming the file, never a silently wrong zero.
CrippenResult crippen(const Molecule& m);

// The exact string every Quantity built from crippen() must carry as its source.
const char* crippenCitation();

}  // namespace biocad::chem
