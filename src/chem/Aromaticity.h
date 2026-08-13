// chem/Aromaticity.h - graph-based aromaticity perception (Daylight-style).
//
// WHY this exists: until now the engine had no aromaticity perception at all -
// it simply trusted the input SMILES. Lowercase `c1ccccc1` was aromatic and the
// exactly equivalent Kekule form `C1=CC=CC=C1` was not, which is a real bug:
// every SMARTS aromatic primitive (a, A, c, n, the `:` bond), every structural
// alert and every biotransformation rule would fire on one spelling of benzene
// and silently miss the other. Perception here reads the GRAPH, so both forms
// give the same answer.
//
// MODEL: per-ring Huckel counting over the SSSR from Rings.h.
//   1. every atom of the ring must be sp2-capable (see the table below); one
//      sp3 centre disqualifies the ring outright (cyclohexane, cyclopentadiene),
//   2. sum the pi-electron contributions of the ring atoms,
//   3. the ring is aromatic when the sum is 4n + 2 for some n >= 0 (2, 6, 10),
//   4. repeat to a fixed point so FUSED systems resolve (see below).
//
// PI-ELECTRON CONTRIBUTION TABLE (and the reason for each row):
//
//   atom in ring                                     e- | why
//   -------------------------------------------------+---+-------------------------
//   C with a ring double bond (order 2 or 1.5)        | 1 | one electron of the
//                                                     |   | shared C=C pi bond lies
//                                                     |   | in this ring
//   C with an exocyclic double bond to O/N/S          | 0 | the pi pair is pulled
//     (carbonyl / thiocarbonyl / imine carbon)        |   | out of the ring by the
//                                                     |   | electronegative atom;
//                                                     |   | the carbon stays sp2
//   C with an exocyclic double bond to a carbon that  | 1 | cross-ring conjugation
//     is itself perceived aromatic (fusion carbon)    |   | at a fused bond; this is
//                                                     |   | the fixed-point rule
//   C with an exocyclic double bond to a non-aromatic | 0 | fulvene-type exocyclic
//     carbon                                          |   | alkene donates nothing
//   C(-) carbanion, no ring double bond               | 2 | filled p orbital
//                                                     |   | (cyclopentadienyl anion)
//   C(+) carbocation, no ring double bond             | 0 | empty p orbital
//                                                     |   | (tropylium)
//   N/P with a ring double bond and no ring H         | 1 | pyridine-type: the lone
//                                                     |   | pair is in the sigma
//                                                     |   | plane, not the pi system
//   N(+)/P(+) in a pi ring bond                       | 1 | pyridinium / thiazolium /
//                                                     |   | isoquinolinium: the cation
//                                                     |   | has no lone pair left to
//                                                     |   | donate, so it is a
//                                                     |   | pyridine-type centre even
//                                                     |   | when it carries three
//                                                     |   | sigma bonds (N-methyl)
//   N/P with an H, or three sigma connections, and    | 2 | pyrrole-type: the lone
//     only single/aromatic ring bonds                 |   | pair completes the sextet
//   N(-) amide/azolate anion                          | 2 | filled p orbital
//   O/S/Se with two ring single (or aromatic) bonds   | 2 | furan/thiophene: one lone
//                                                     |   | pair enters the pi system
//   O(+)/S(+) with a ring double bond                 | 1 | pyrylium/thiopyrylium
//   B (three connections)                             | 0 | empty p orbital (borole,
//                                                     |   | borazine-type rings)
//   anything else                                     | - | not sp2: ring rejected
//
// FUSED SYSTEMS: after the first pass, rings are re-evaluated until nothing
// changes. This is required for Kekule input, because in `C1=CC=C2C=CC=CC2=C1`
// (naphthalene) the second ring's two fusion carbons have their double bonds
// pointing into the FIRST ring, so in isolation that ring counts only 4
// electrons. Once the first ring is known aromatic, the fusion carbons' partners
// are aromatic and row 3 of the table upgrades them to 1 each, giving 6.
// Naphthalene and indole therefore come out fully aromatic from either spelling.
//
// WHAT THIS MODEL GETS WRONG - stated plainly, because it is a heuristic RING
// model, not an electron-flow perception:
//   * 2-pyridone / cytosine / caffeine's pyrimidinedione ring are reported
//     AROMATIC. The exocyclic C=O contributes 0 and the amide N contributes 2,
//     which sums to 6. This matches RDKit's default model and ChEMBL-style
//     conventions, but a chemist may reasonably call these amides.
//   * Aromaticity delocalised over a PERIMETER rather than over an SSSR ring is
//     missed outright. Verified example: azulene (c1ccc2cccc2cc1) is reported
//     NOT aromatic, because its SSSR is a 5-ring and a 7-ring counting 5 and 7
//     electrons respectively, while the real 10-electron system runs around the
//     bicyclic perimeter, which is not an SSSR ring. Macrocyclic annulenes
//     (e.g. [18]annulene) fail the same way. Fixing this needs perimeter/
//     electron-flow perception, which this deliberately is not.
//   * No anti-aromaticity reasoning: a 4n ring is simply "not aromatic", never
//     flagged as destabilised (cyclobutadiene).
//   * Charge-separated and mesomeric-betaine rings, N-oxides, and metal-
//     coordinated rings are not modelled.
//   * The SSSR itself is not unique for cage systems (see Rings.h), so a cage
//     ring assignment is one valid choice, not the only one.
//
// This function sets PERCEIVED FLAGS ONLY (Atom::aromatic, Bond::aromatic). It
// never produces a data::Quantity, so no provenance tier is involved: a flag is
// a structural perception, not a measured or predicted number.
#pragma once

#include "chem/Molecule.h"
#include "chem/Rings.h"

namespace biocad::chem {

// Reset every Atom::aromatic / Bond::aromatic flag and set them from the graph
// using the model documented above. `info` must be the SSSR of `m`.
// A bond is flagged aromatic only when it lies inside a ring that was perceived
// aromatic (so the single bond joining the two rings of biphenyl stays plain).
void perceiveAromaticity(Molecule& m, const RingInfo& info);

// Convenience: perceive rings, annotate them, and perceive aromaticity.
// Returns the SSSR so callers can keep it instead of recomputing.
RingInfo perceiveRingsAndAromaticity(Molecule& m);

// Optional post-step for code that compares bond orders across input spellings
// (e.g. a canonical SMILES writer): set Bond::order to 1.5 on every perceived
// aromatic bond, so Kekule and lowercase input converge. Deliberately separate
// from perceiveAromaticity, because doing it unconditionally would destroy the
// Kekule structure that valence, H-count and carbonyl detection rely on.
void normalizeAromaticBondOrders(Molecule& m);

}  // namespace biocad::chem
