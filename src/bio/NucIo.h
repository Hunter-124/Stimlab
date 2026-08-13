#pragma once

// FASTA and GenBank flat-file I/O for nucleic-acid records.
//
// The GenBank reader exists because a plasmid map is worth nothing without its
// features, and the feature table is the part of the format that is actually
// hard: locations are a small recursive language (join, order, complement,
// partial bounds), qualifier values are quoted and wrap across lines, and the
// continuation column is significant. A reader that handles only "n..m" opens a
// plasmid and silently loses every spliced CDS in it.
//
// Partial locations (<n, m>) are not silently squared off: the partiality is
// recorded as a qualifier on the feature and a warning on the record, because a
// CDS whose start is missing does not translate from its first stored base.

#include <string>
#include <vector>

#include "data/Nucleic.h"

namespace biocad::bio {

// ------------------------------------------------------------------- FASTA

// One record per '>' header. The id is the first whitespace-delimited token of
// the header, the rest is the description. Sequence characters are uppercased and
// validated as IUPAC; anything else is dropped with a warning on the record.
std::vector<NucRecord> readNucFasta(const std::string& text);
std::string writeNucFasta(const NucRecord& record, std::size_t lineWidth = 60);

// ----------------------------------------------------------------- GenBank

// Reads one or more records separated by "//". Understands LOCUS (including the
// circular/linear topology and the molecule type), DEFINITION, ACCESSION,
// VERSION, the FEATURES table and ORIGIN. Unrecognised header keywords are
// skipped rather than treated as errors: the format has open-ended sections and
// refusing to open a file because of a COMMENT block would be useless.
std::vector<NucRecord> readGenBank(const std::string& text);

// Writes a record that reads back with identical id, description, topology,
// sequence and features. It is not a byte-for-byte reproduction of an NCBI file -
// nothing outside the DTO can be reproduced - but it round-trips structurally,
// which is the property the tests assert.
std::string writeGenBank(const NucRecord& record);

// Sniffs FASTA vs GenBank on the first non-blank line so a caller with a file of
// unknown provenance has exactly one entry point.
std::vector<NucRecord> readNucleic(const std::string& text);

// ------------------------------------------------------- location expressions
//
// Exposed because the location language is independently testable, and because
// anything that writes a feature table needs the inverse.

struct ParsedLocation {
    std::vector<std::pair<int, int>> parts;      // [begin, end), 0-based, ascending
    Strand                           strand = Strand::Forward;
    bool                             ordered = false;   // order(...) not join(...)
    bool                             partial5 = false;  // <n
    bool                             partial3 = false;  // m>
    bool                             betweenBases = false;   // n^m
    std::string                      remoteAccession;       // ACC:1..10, unsupported
};

ParsedLocation parseLocation(const std::string& text, std::vector<std::string>* warnings = nullptr);
std::string formatLocation(const std::vector<std::pair<int, int>>& parts, Strand strand,
                           bool ordered = false);

}  // namespace biocad::bio
