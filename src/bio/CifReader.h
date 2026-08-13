#pragma once

// Reader for mmCIF (PDBx/mmCIF), via a minimal STAR/CIF tokeniser.
//
// This reader is MANDATORY, not a nicety. The legacy PDB format has been frozen at v3.30 since
// 21 November 2012 and the wwPDB issues NO PDB file at all for an entry with more than 62
// chains, more than 99999 atoms, a multi-character chain id, or a B-factor above 999.99. A
// PDB-only application therefore cannot open a large ribosome, a big virus capsid, or many
// cryo-EM depositions - and it fails by finding nothing to download rather than by reporting an
// unsupported structure, which is the worst possible failure mode.
//
// Categories read: _atom_site, _entity_poly, _struct_conf, _struct_sheet_range, _struct_conn.
// Columns are resolved BY TAG NAME. mmCIF loops carry no guaranteed column order, and files
// from different producers genuinely differ, so positional parsing works on the file you tested
// and corrupts the next one.

#include <filesystem>
#include <string>

#include "bio/Annotations.h"
#include "bio/Structure.h"

namespace biocad::bio {

Structure readCif(const std::string& text, std::string id = {}, std::string source = {},
                  Annotations* annotations = nullptr);
Structure readCifFile(const std::filesystem::path& path, Annotations* annotations = nullptr);

}  // namespace biocad::bio
