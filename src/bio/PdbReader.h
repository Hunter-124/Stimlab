#pragma once

// Reader for the legacy fixed-column PDB format (v3.30).
//
// The format is FIXED COLUMN, not whitespace separated: several fields are legally adjacent
// with no separator (a five-digit resSeq beside a four-character atom name), so any
// tokenising reader corrupts real files. Every field below is sliced by column index.
//
// The reader is recoverable by design: a malformed or truncated line appends a message to
// Structure::warnings and is skipped. It never throws on bad input, because a half-downloaded
// or hand-edited file must still open with whatever it did contain.

#include <filesystem>
#include <string>

#include "bio/Annotations.h"
#include "bio/Structure.h"

namespace biocad::bio {

// `id` seeds Structure::id when the file carries no HEADER record. `annotations`, when non-null,
// receives SEQRES/HELIX/SHEET/SSBOND; passing null skips building them.
Structure readPdb(const std::string& text, std::string id = {}, std::string source = {},
                  Annotations* annotations = nullptr);
Structure readPdbFile(const std::filesystem::path& path, Annotations* annotations = nullptr);

}  // namespace biocad::bio
