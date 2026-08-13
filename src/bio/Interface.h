// bio/Interface.h - protein-protein interfaces and the geometric alanine scan.
//
// WHAT IS EXACT HERE. Contacts, buried area, Levy's support/core/rim partition and
// the interaction inventory (hydrogen bonds, salt bridges, hydrophobic contacts,
// pi-pi, cation-pi, disulfides) are GEOMETRY over the coordinates supplied. They are
// Provenance::Measured, they carry the SASA parameter string that makes them
// reproducible, and none of them is an energy.
//
// WHY THE ALANINE SCAN IS UNIT-FREE. Truncating a side chain beyond C-beta and
// re-measuring how much interface disappears is a rank-ordering of which side chains
// hold the interface together. It is NOT a binding free energy: no solvation term, no
// entropy, no relaxation of the remaining structure, no electrostatic model. So the
// impact is a Provenance::Heuristic Quantity with an EMPTY unit, and makeQuantity()
// throws if anyone tries to attach kcal/mol to it. That throw is the feature.
//
// AND WHY THE BENCHMARK CORRELATION IS NotComputed. A geometric hotspot proxy is
// only worth its measured rank correlation, and this build has no measured
// antibody/antigen mutation set on disk. Hard-coding "Spearman 0.55" from memory
// would be the exact failure this file exists to avoid, so
// AlanineScanReport::benchmarkSpearman is notComputed("a measured benchmark subset")
// until a real set is shipped and the correlation is actually computed.
#pragma once

#include <string>
#include <vector>

#include "bio/Structure.h"
#include "data/Biologics.h"

namespace biocad::bio {

struct InterfaceOptions {
    double contactCutoff = 4.5;      // A, heavy atom to heavy atom
    double hydrogenBondCutoff = 3.5; // A, N/O donor to N/O acceptor, geometry-free
    double saltBridgeCutoff = 4.0;   // A, charged group to charged group
    double piStackingCutoff = 5.5;   // A, aromatic ring centroid to centroid
    double cationPiCutoff = 6.0;     // A, cation to ring centroid
    double disulfideCutoff = 2.5;    // A, SG to SG
    // Levy 2010 thresholds on relative accessibility.
    double burialThreshold = 0.25;
    // Antibody chains, so CDR contacts, epitope and paratope can be reported. Empty
    // means those fields stay empty rather than guessing which side is the antibody.
    std::vector<std::string> antibodyChains;
};

// Chain-id list parsing: "H,L" or "HL" both mean {H, L}. Multi-character mmCIF ids
// must use commas, and a chain that is not in the model is reported as a warning.
std::vector<std::string> parseChainList(const std::string& spec);

InterfaceReport interfaceOf(const Structure& complex, const std::string& chainsA,
                            const std::string& chainsB, const InterfaceOptions& opts = {});

AlanineScanReport alanineScan(const Structure& complex, const std::string& chainsA,
                              const std::string& chainsB, const InterfaceOptions& opts = {});

}  // namespace biocad::bio
