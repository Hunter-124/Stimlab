#pragma once

// FASTA reading and writing.
//
// The header line is not free text in practice: UniProt and RCSB each impose a grammar, and the
// UniProt one cannot be parsed by splitting on spaces because the protein name itself contains
// spaces and GN= is simply absent when the gene is unknown. It is parsed right-to-left on the
// trailing KEY=VALUE tokens instead - see Fasta.cpp.

#include <filesystem>
#include <string>
#include <vector>

namespace biocad::bio {

enum class FastaHeaderKind {
    Plain,     // ">anything at all"
    UniProt,   // ">sp|P12345|NAME_HUMAN Protein name OS=... OX=... [GN=...] PE=... SV=..."
    Rcsb       // ">4M48_1|Chain A[auth B]|Protein name|Homo sapiens (9606)"
};

struct FastaRecord {
    std::string header;        // the raw header line without the leading '>'
    std::string sequence;      // residues, uppercased, gaps and stops retained
    FastaHeaderKind kind = FastaHeaderKind::Plain;

    // UniProt fields; empty when absent or when the header is not UniProt.
    std::string database;      // "sp" or "tr"
    std::string accession;     // "P12345"
    std::string entryName;     // "NAME_HUMAN"
    std::string proteinName;   // may contain spaces
    std::string organism;      // OS=
    std::string taxonId;       // OX=
    std::string geneName;      // GN=, often absent
    std::string evidence;      // PE=
    std::string version;       // SV=

    // RCSB fields; empty when absent.
    std::string entryId;       // "4M48_1"
    std::vector<std::string> chainIds;   // "A", "B" (auth ids when the header supplies them)
};

// Any residue letter accepted by the reader: the 20 standard residues plus B (Asx), Z (Glx),
// J (Leu/Ile), X (unknown), U (selenocysteine), O (pyrrolysine), '*' (stop) and '-' (gap).
bool isSequenceCharacter(char c);

std::vector<FastaRecord> readFasta(const std::string& text, std::vector<std::string>* warnings = nullptr);
std::vector<FastaRecord> readFastaFile(const std::filesystem::path& path,
                                       std::vector<std::string>* warnings = nullptr);

// Writes with a fixed line width; 60 is the de-facto standard and what UniProt emits.
std::string writeFasta(const std::vector<FastaRecord>& records, std::size_t lineWidth = 60);

}  // namespace biocad::bio
