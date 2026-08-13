// chem/Perceive.h - parse a SMILES into a graph that is ready to be computed on.
//
// parseSmiles() gives you a CONNECTIVITY graph. It sets Atom::aromatic only for
// lowercase input, because that is all the input told it. Ring membership from the
// parser is a crude bridge flag, not SSSR.
//
// Several descriptors and every SMARTS query branch on those flags, so computing
// straight off parseSmiles() makes the answer depend on how the molecule was
// SPELLED: caffeine's TPSA measures 56.22 A^2 from the Kekule spelling against
// 60.26 A^2 once aromaticity is perceived, and `c1ccccc1` and `C1=CC=CC=C1` would
// disagree about being aromatic at all.
//
// A descriptor that depends on the spelling of its input is not a descriptor. Use
// this instead of parseSmiles() anywhere a number or a substructure test follows.
#pragma once

#include <optional>
#include <string_view>

#include "chem/Aromaticity.h"
#include "chem/Molecule.h"
#include "chem/Rings.h"
#include "chem/Smiles.h"

namespace biocad::chem {

// Parse, then perceive SSSR rings and graph aromaticity. Returns nullopt only when
// the SMILES itself is malformed - perception never fails, it only annotates.
inline std::optional<Molecule> parsePerceived(std::string_view smiles) {
    auto m = parseSmiles(smiles);
    if (!m) return std::nullopt;
    const RingInfo rings = perceiveRings(*m);
    annotateRings(*m, rings);
    perceiveAromaticity(*m, rings);
    return m;
}

// Perceive in place, for a graph that was built by hand rather than parsed (the
// structure sketcher). Assumes implicit hydrogens are already assigned.
inline RingInfo perceiveInPlace(Molecule& m) {
    const RingInfo rings = perceiveRings(m);
    annotateRings(m, rings);
    perceiveAromaticity(m, rings);
    return rings;
}

}  // namespace biocad::chem
