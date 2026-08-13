// packs/Pack.h - data packs: the compound and target catalogs, as data.
//
// A pack is one versioned JSON document describing compounds and receptor
// targets. Built-in packs ship beside the executable under assets/packs; user
// packs live under %APPDATA%/BioCAD/packs and override built-ins by pack id.
// Nothing in this file knows about CNS pharmacology, or about any other domain:
// the domain is whatever the loaded packs say it is.
//
// SAFETY SCOPE: a pack describes what a compound IS and where a ligand binds.
// There is deliberately no synthesis, route, precursor or manufacturability
// field, and a loader that encountered one would ignore it.
#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "contracts/IDockingBackend.h"
#include "data/Domain.h"

namespace biocad::packs {

// The only schema version this build understands. An unknown version is a load
// error surfaced in the UI, never a silent skip: silently ignoring a pack the
// user just wrote is indistinguishable from the app being broken.
inline constexpr int kSchemaVersion = 1;

// External database identifiers. Later phases join on these; nothing here
// dereferences them.
struct Xrefs {
    std::string chembl;      // e.g. "CHEMBL113"
    std::string inchikey;    // e.g. "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
    std::string uniprot;     // targets only
    long long   pubchemCid = 0;  // 0 = absent
};

// One catalog compound. A pack authors identity and metadata only: id, name,
// SMILES, drug class, legal status, notes and cross-references.
//
// A descriptor has exactly one source, the chem engine. Molecular weight, logP,
// TPSA, H-bond counts and rotatable bonds are computed from the SMILES by
// chem::Descriptors at the point of use, so a pack cannot author them: two
// sources of truth for a molecular weight would silently diverge, and the
// authored one would lose without saying so.
struct PackCompound {
    std::string id;
    std::string name;
    std::string smiles;
    std::string drugClass;
    std::string legalStatus;
    std::string notes;
    Xrefs       xrefs;

    // Identity and metadata only. Every numeric field of the Molecule is left at
    // zero for the caller to compute from `smiles`.
    [[nodiscard]] Molecule molecule() const;
};

// One receptor target. A target without a box is a coverage gap, not a docking
// target: it is listed so the user can see what is missing rather than being
// handed a fabricated binding site.
struct PackTarget {
    ReceptorTarget target;              // id, name, pdb, box
    bool           hasBox = false;      // false => box is meaningless, do not dock
    bool           headline = false;    // prepare this receptor during provisioning
    std::string    uniprot;
    std::string    chembl;
    std::vector<std::string> panels;    // e.g. {"safetyscreen44"}
};

struct Pack {
    int         schemaVersion = 0;
    std::string id;
    std::string title;
    std::string description;
    std::string sourcePath;             // file it was loaded from ("<embedded>" for literals)
    bool        builtin = false;        // shipped beside the exe rather than user-authored
    std::vector<PackCompound> compounds;
    std::vector<PackTarget>   targets;
};

// The merged result of loading a directory (or the built-in set). `packs` is in
// resolution order; `errors` is everything that failed, and is rendered in the
// Presets panel rather than swallowed.
struct LoadReport {
    std::vector<Pack>        packs;
    std::vector<std::string> errors;

    // Merged views, with later packs overriding earlier ones by pack id. A
    // duplicate compound/target id across packs is recorded in `errors` (which is
    // why these are non-const) and the later one is dropped, never silently
    // resolved behind the user's back.
    [[nodiscard]] std::vector<PackCompound> compounds();
    [[nodiscard]] std::vector<PackTarget>   targets();
    [[nodiscard]] std::optional<PackTarget> findTarget(const std::string& nameOrId) const;
};

// Parse one pack document. Throws Error (core/Error.h) when schemaVersion is
// absent, unknown, or the document is structurally wrong.
Pack parse(const nlohmann::json& j, std::string sourcePath = "<embedded>");

// Parse a pack from a JSON string, for embedded literals.
Pack parseString(const std::string& text, std::string sourcePath = "<embedded>");

// Load every *.json in `dir`, sorted by filename for determinism. A missing
// directory is not an error (it yields an empty report); an unparseable file is.
LoadReport loadFrom(const std::filesystem::path& dir, bool builtin = false);

// The shipping resolution order: assets/packs beside the executable, then the
// user overlay under AppPaths::instance().root()/"packs". Later wins by pack id.
LoadReport loadBuiltin();

// Directory the built-in packs are read from (exe-relative, with a source-tree
// fallback so a dev build works from the build directory).
std::filesystem::path builtinPackDir();

}  // namespace biocad::packs
