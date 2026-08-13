#pragma once

// Records that describe a structure but are not part of the atom hierarchy: SEQRES, HELIX,
// SHEET, SSBOND and their mmCIF equivalents (_entity_poly, _struct_conf, _struct_sheet_range,
// _struct_conn). They live beside bio::Structure rather than inside it because the Structure
// contract is shared with the alignment, superposition and scoring slices and must stay a pure
// coordinate model; annotations are optional and format-specific.

#include <string>
#include <vector>

namespace biocad::bio {

// The deposited (full) sequence for a chain. It differs from the observed sequence whenever a
// loop is disordered, which is exactly why alignment must not be driven off coordinates alone.
struct SeqResChain {
    std::string chainId;
    std::vector<std::string> residueNames;   // three-letter codes, deposition order
    // mmCIF states the sequence one-letter (_entity_poly.pdbx_seq_one_letter_code) while PDB
    // states it three-letter; both readers fill both so the two formats can be compared.
    std::string oneLetterCode;
};

struct HelixRecord {
    std::string id;
    std::string chainId;
    int  startSeqId = 0;
    char startInsertionCode = ' ';
    int  endSeqId = 0;
    char endInsertionCode = ' ';
};

struct StrandRecord {
    std::string sheetId;
    std::string chainId;
    int  startSeqId = 0;
    char startInsertionCode = ' ';
    int  endSeqId = 0;
    char endInsertionCode = ' ';
};

// A disulfide is a covalent link between two residues that may sit in different chains, so it
// cannot be represented as a per-residue property.
struct DisulfideRecord {
    std::string chainId1;
    int  seqId1 = 0;
    char insertionCode1 = ' ';
    std::string chainId2;
    int  seqId2 = 0;
    char insertionCode2 = ' ';
};

struct Annotations {
    std::vector<SeqResChain>     seqres;
    std::vector<HelixRecord>     helices;
    std::vector<StrandRecord>    strands;
    std::vector<DisulfideRecord> disulfides;
};

}  // namespace biocad::bio
