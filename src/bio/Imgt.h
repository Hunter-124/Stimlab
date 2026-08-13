// bio/Imgt.h - IMGT unique numbering of an immunoglobulin V-DOMAIN, and the
// table-driven views of it.
//
// WHY THERE IS NO PYTHON HERE. The reference implementation of antibody
// numbering (ANARCI) is BSD-3 but needs a Python runtime and HMMER. What ANARCI
// actually consumes is IMGT's own reference data, and IMGT publishes the V-REGION
// amino-acid sequences ALREADY GAPPED into the unique numbering: in
// assets/packs/biologics/imgt-reference.json, index i of a `gapped` string IS
// IMGT position i, and every functional entry carries Cys at index 104. So the
// position transfer is a table lookup through one alignment, not a re-derivation,
// and the alignment engine is Phase 5's Gotoh (bio/Align.h) - there is exactly one
// aligner in this tree.
//
// WHY A FAILED ANCHOR RETURNS NOTHING. IMGT 23 (1st-CYS), 41 (CONSERVED-TRP), 89
// (hydrophobic), 104 (2nd-CYS) and 118 (J-PHE/J-TRP) are the positions the whole
// numbering hangs off. If one of them does not land on the residue it must, the
// alignment found something that is not a V-DOMAIN in the register it thinks it
// is, and every downstream number - CDR lengths, liability positions, scheme
// conversion - is then confidently wrong. Wrong numbering is worse than none, so
// AbDomain::numbered stays false, AbDomain::residues stays EMPTY, and the failures
// are named. A T-cell receptor beta chain is rejected the same way.
//
// WHY "closestGermlineSet" IS NOT A SPECIES. It is the best-scoring entry in the
// pack, and the pack is human (plus what the retrieval could reach). A camelid VHH
// scores best against a human IGHV3 gene; that is a similarity result, not an
// identification, and nothing in this tree may render it as one.
#pragma once

#include <optional>
#include <string>
#include <vector>

#include "bio/Structure.h"
#include "data/Biologics.h"

namespace biocad::bio {

// One IMGT region and the closed position interval it owns.
struct ImgtRegion {
    std::string name;
    int         first = 0;
    int         last = 0;
};

// A conserved position and the residues that satisfy it.
struct ImgtAnchor {
    int         position = 0;
    std::string accept;    // set of acceptable one-letter codes
    std::string label;
};

// How a CDR shrinks and grows. `deletionOrder` is the order IMGT empties slots as
// the loop shortens (centre outward); `insertionAnchors` are the two central
// positions that carry insertion codes when the loop is LONGER than the region,
// left ascending and right descending, which is what keeps 110 and 113 fixed while
// a 20-residue CDR3 grows between them.
struct ImgtCdrLayout {
    std::string      region;
    std::vector<int> deletionOrder;
    std::vector<int> insertionAnchors;
};

struct ImgtVGene {
    std::string gene;        // "IGHV3-23*01"
    std::string geneClass;   // "IGHV" | "IGKV" | "IGLV" | "TRBV"
    std::string gapped;      // IMGT-gapped amino acids, '.' = gap slot
    std::string ungapped;    // gapped with the gaps removed, what alignment sees
    std::vector<int> column; // ungapped index -> IMGT position
};

struct ImgtJGene {
    std::string gene;
    std::string geneClass;
    std::string sequence;
    int         fr4Offset = -1;   // index of IMGT 118 inside `sequence`, -1 if unknown
};

// The loaded reference. `ok` false means nothing can be numbered; `errors` says why
// and every caller surfaces it rather than reporting "not an antibody".
struct ImgtReference {
    bool                       ok = false;
    std::string                source;
    std::string                licence;
    std::string                citation;
    std::vector<ImgtRegion>    regions;
    std::vector<ImgtAnchor>    anchors;
    std::vector<ImgtCdrLayout> cdrLayouts;
    std::vector<ImgtAnchor>    vhhHallmark;
    std::vector<ImgtVGene>     vGenes;
    std::vector<ImgtJGene>     jGenes;
    std::vector<std::string>   errors;
};

// Loaded once from assets/packs/biologics/imgt-reference.json.
const ImgtReference& imgtReference();

// Numbering options. The bit-score floor is the only defence against numbering a
// domain from a family the pack does not contain (the retrieval could not reach
// TRAV/TRAJ, so a TCR alpha chain is caught here rather than by name).
struct NumberingOptions {
    double minVBitScore = 40.0;
    // Residues past IMGT 128 are a constant domain, not the V-DOMAIN, and are
    // reported as a warning instead of being numbered.
    bool   warnOnTrailingResidues = true;
};

// IMGT-numbers one V-DOMAIN. Always returns an AbDomain: on refusal, `numbered` is
// false, `residues` is empty and `anchorFailures` names every reason.
AbDomain numberDomain(const std::string& sequence, const NumberingOptions& opts = {});

// Re-renders an IMGT-numbered domain in another scheme, purely by table lookup from
// assets/packs/biologics/scheme-maps.json. A scheme with no obtainable table is
// REFUSED: the returned domain has numbered == false and a warning naming the
// missing table, because a guessed correspondence renumbers residues plausibly and
// wrongly. Converting back to Imgt restores the IMGT positions exactly.
AbDomain convertScheme(const AbDomain& domain, NumberingScheme to);

// "111.3" for position 111 insertion 3, "111" when there is no insertion code.
std::string positionLabel(const NumberedResidue& r);

// Relative side-chain accessibility per numbered residue, or an empty vector when
// the chain is not in the structure. Index is parallel to `domain.residues`.
std::vector<double> relativeSasaFor(const AbDomain& domain, const Structure& s,
                                    const std::string& chainId);

}  // namespace biocad::bio
