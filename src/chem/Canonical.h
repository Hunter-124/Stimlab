// chem/Canonical.h - canonical ranking, canonical SMILES writer, and a graph hash.
//
// WHY this exists: without a canonical form, two identical molecules entered as
// different SMILES are unrecognisably different strings, so a library cannot be
// deduplicated, a cache cannot be keyed, and "are these the same structure?" has
// no answer. The only other writer in the tree (sketchToSmiles in the UI) walks
// the graph in input order and is therefore not canonical.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "chem/Molecule.h"

namespace biocad::chem {

// Morgan-style canonical ranking: a total order over the atoms.
//
// Initial invariants (element, degree, charge, total H, aromaticity, ring
// membership) are refined by iterated neighbour-rank hashing until the partition
// stops changing; remaining ties are broken one atom at a time, deterministically
// (lowest class first, then lowest atom index), refining after each break so the
// result is always a total order.
//
// HONESTY ABOUT THE STRENGTH OF THIS: breaking ties by atom index makes the
// ranking canonical for a GIVEN atom ordering. Full graph-automorphism
// canonicalisation - a ranking provably invariant under every relabelling of the
// input, as a proper refinement/backtracking canonicaliser (nauty, RDKit's
// canonical ranker) provides - is a strictly stronger property that this does NOT
// claim. What is made order-independent, and is what callers actually depend on,
// is canonicalSmiles(): it roots its walk at the canonical rank rather than at
// atom 0, so atoms that the refinement cannot separate are interchangeable and
// produce the same string. tests/test_chem_canonical.cpp verifies that
// experimentally over re-rooted inputs.
[[nodiscard]] std::vector<int> canonicalRanks(const Molecule& m);

// Canonical SMILES. DFS from the lowest-ranked atom, neighbours in canonical-rank
// order, ring-closure digits allocated in order of first need, disconnected
// components joined with '.' and ordered by their own canonical string.
//
// The contract is round-tripping: parseSmiles(canonicalSmiles(m)) yields a graph
// with the same canonical string. Aromaticity is taken from the molecule as
// annotated (see perceiveAromaticity in chem/Aromaticity.h) - this writer does
// not perceive it, so a Kekule graph is written Kekule.
[[nodiscard]] std::string canonicalSmiles(const Molecule& m);

// 64-bit hash of the canonical form, suitable as a cache key.
//
// It is a HASH, not an identity proof: a collision is a wrong cache hit, i.e. one
// molecule silently answering for another. Callers that cannot tolerate that must
// compare canonicalSmiles() as well. The hash therefore folds in the atom and
// bond counts alongside the canonical string, so the cheapest structural
// mismatches cannot collide at all.
[[nodiscard]] std::uint64_t graphHash(const Molecule& m);

}  // namespace biocad::chem
