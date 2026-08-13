// chem/Rings.h - ring perception: Smallest Set of Smallest Rings (SSSR).
//
// WHY this exists: before this file the only ring information in the tree came
// from the SMILES parser's bridge detection (Smiles.cpp), which answers the
// single question "is this bond part of some cycle?" and nothing else. Every
// SMARTS ring primitive (R, R0, !R, rN), every aromaticity model, and every
// published biotransformation rule needs actual rings with actual sizes, so the
// real cycle basis is computed here once and cached in a RingInfo.
//
// HOW: the size of the cycle space is fixed by Euler's formula,
//
//     dim = |bonds| - |atoms| + |connected components|
//
// so we know up front exactly how many rings a correct answer contains (0 for
// any acyclic or disconnected-acyclic input). Candidate rings are generated as
// Horton's set - for every atom r and every bond (x,y), the cycle
// SP(r,x) + (x,y) + SP(y,r) when the two shortest paths are vertex-disjoint -
// which is proven to contain a minimum-weight cycle basis. Candidates are then
// sorted smallest-first and accepted greedily only when linearly independent of
// the already-accepted rings over GF(2) on their bond-incidence vectors
// (Gaussian elimination on a bitset). The loop stops as soon as the accepted
// count reaches `dim`, so it always terminates and never over-reports.
//
// HONESTY NOTE: an SSSR is NOT unique for fused/cage systems. Cubane is the
// standard example: it has six four-membered faces but a cycle space of
// dimension 12 - 8 + 1 = 5, so any SSSR must discard one face, and which face
// is discarded is arbitrary. This function returns ONE valid SSSR, chosen
// deterministically (candidates are ordered by ring size, then by their
// canonicalised atom sequence, so the same molecule always yields the same
// rings). That is exactly what ring-size queries and aromaticity models need;
// code that needs *every* small ring wants the relevant-cycles set instead,
// which this deliberately does not compute.
#pragma once

#include <vector>

#include "chem/Molecule.h"

namespace biocad::chem {

// One valid SSSR. atomRings[i] and bondRings[i] describe the same ring:
// atomRings[i] is the cycle walked in order (each consecutive pair, plus the
// wrap-around pair, is bonded); bondRings[i] holds the same ring's bond indices
// in the matching walk order. Rings are ordered smallest-first.
struct RingInfo {
    std::vector<std::vector<int>> atomRings;
    std::vector<std::vector<int>> bondRings;

    [[nodiscard]] std::size_t count() const { return atomRings.size(); }
};

// Compute the SSSR of `m`. Pure: does not modify the molecule.
RingInfo perceiveRings(const Molecule& m);

// Write the perception back onto the graph: Atom::inRing and Bond::inRing are
// cleared and then set from `info`. This replaces the parser's crude bridge
// flag with the real thing (the two agree for simple rings, but the parser's
// flag says nothing about ring size or membership counts).
void annotateRings(Molecule& m, const RingInfo& info);

// Size of the SMALLEST perceived ring containing `atomIndex`, or 0 if the atom
// is in no perceived ring. This is the SMARTS `rN` primitive's input.
int ringSizeOf(const RingInfo& info, int atomIndex);

// Is `atomIndex` a member of some perceived ring of exactly `size` atoms?
// Unlike ringSizeOf this looks at every ring the atom is in, so a fusion atom
// of an indole reports true for both 5 and 6.
bool inRingOfSize(const RingInfo& info, int atomIndex, int size);

// Number of perceived rings containing `atomIndex` (SMARTS `Rn`).
int ringCountOf(const RingInfo& info, int atomIndex);

}  // namespace biocad::chem
