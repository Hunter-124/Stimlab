// bio/Connectivity.h - the bond graph a coordinate file does not contain.
//
// A PDB or mmCIF entry is positions and names. CONECT is optional and, when present, only
// covers non-standard groups; _struct_conn covers links, not intra-residue chemistry. So a
// renderer, an interface analysis or a cartoon has to obtain bonds from somewhere else, and
// that somewhere is assets/packs/structure/residue-templates.json: connectivity per standard
// residue keyed by ATOM NAME, which is the one piece of chemical identity the file does state.
//
// WHY NOT DISTANCE. A covalent-radius distance cutoff over a crystal structure invents bonds:
// it bridges the two halves of a disordered side chain, it fuses stacked bases, and at 2.0 A
// resolution it drops real bonds whose refined length is long. A template match, by contrast,
// is either present or reported missing. A distance pass exists (ConnectivityOptions::
// inferByDistance) because an unparameterised ligand has no other option, but it is OFF by
// default and every bond it creates is BondKind::DistanceInferred so a caller can never mistake
// a guess for a template.
//
// WHAT IS NEVER SILENT. An unrecognised residue name is reported with its occurrence count. A
// template bond whose atom is absent (a disordered side chain refined only to CB) is counted.
// A missing loop shows up as an explicit ChainGap with the measured C..N distance, and NO bond
// is created across it - a chain break rendered as a bond is a fabricated covalent link between
// residues that may be 30 A apart.
#pragma once

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

#include "bio/Structure.h"

namespace biocad::bio {

// Chemical bond order from the template. Aromatic is kept as its own value rather than being
// Kekule-alternated: the alternation is a drawing convention and choosing one would assert a
// tautomer the coordinate file never stated.
enum class BondOrder { Single, Double, Aromatic };

// Where a bond came from. This is provenance for geometry: a caller that wants only facts can
// filter out DistanceInferred, and the UI can label it.
enum class BondKind {
    Template,         // named pair from residue-templates.json
    PeptideLink,      // C(i)-N(i+1), inside the length window
    NucleicLink,      // O3'(i)-P(i+1)
    Disulfide,        // SG-SG by distance, reported separately from the templates
    DistanceInferred  // opt-in fallback only; never produced by default
};

// An atom addressed inside one Model. Indices, not pointers: the bond list outlives any
// reference into a vector that a caller may still be filling.
struct AtomRef {
    int chain = 0;
    int residue = 0;
    int atom = 0;
};

struct StructureBond {
    AtomRef   a;
    AtomRef   b;
    BondOrder order = BondOrder::Single;
    BondKind  kind = BondKind::Template;
};

// A break in a polymer chain: two consecutive residues in the file whose linking atoms are too
// far apart (or absent) to be bonded. `distance` is negative when an atom was missing entirely,
// which is a different failure from "present but far".
struct ChainGap {
    std::string chainId;
    int  fromSeqId = 0;
    char fromInsertionCode = ' ';
    int  toSeqId = 0;
    char toInsertionCode = ' ';
    double distance = -1.0;
    std::string reason;
};

struct UnknownResidue {
    std::string name;
    int count = 0;
    int atoms = 0;
};

struct ConnectivityDiagnostics {
    std::vector<UnknownResidue> unknownResidues;
    std::vector<ChainGap>       gaps;
    std::size_t unbondedAtoms = 0;         // atoms that ended up in no bond at all
    std::size_t missingTemplateAtoms = 0;  // template bonds skipped: an endpoint was absent
    std::size_t solventResidues = 0;       // HOH and friends: correctly bond-free, not unknown
    std::size_t monatomicResidues = 0;     // ions: one atom, no bond to make, not unknown
    std::size_t templateBonds = 0;
    std::size_t linkBonds = 0;
    std::size_t disulfides = 0;
    std::size_t inferredBonds = 0;
    std::vector<std::string> warnings;     // human-readable, ready for Structure::warnings
};

struct ConnectivityResult {
    std::vector<StructureBond>  bonds;
    ConnectivityDiagnostics     diagnostics;
};

struct ConnectivityOptions {
    // Opt-in distance fallback for residues with no template. Off by default; see the header
    // comment. `inferMaxAngstrom` is a plain cutoff, deliberately not element-aware, because an
    // element-aware cutoff would look like chemistry.
    bool   inferByDistance = false;
    double inferMaxAngstrom = 1.95;
    bool   findDisulfides = true;
};

// One residue's template: the heavy-atom bonds of the standard residue, by name.
struct TemplateBond {
    std::string a;
    std::string b;
    BondOrder   order = BondOrder::Single;
};

struct ResidueTemplate {
    std::string name;
    std::string kind;   // "protein" | "rna" | "dna"
    std::vector<TemplateBond> bonds;
};

struct LinkageRule {
    std::string from, to;
    double minAngstrom = 0.0;
    double maxAngstrom = 0.0;
};

struct ResidueTemplatePack {
    bool ok = false;
    int  schemaVersion = 0;
    std::vector<ResidueTemplate> residues;
    LinkageRule protein, nucleic, disulfide;
    std::vector<std::string> solvent;
    std::vector<std::string> errors;

    // Resolved through the alias table; null for an unrecognised name.
    [[nodiscard]] const ResidueTemplate* find(const std::string& residueName) const;
    [[nodiscard]] bool isSolvent(const std::string& residueName) const;

    std::vector<std::pair<std::string, std::string>> aliases;
};

// Loaded once from assets/packs/structure/residue-templates.json. A load failure leaves
// ok == false and errors populated; connect() then reports every residue as unknown rather
// than quietly bonding nothing.
const ResidueTemplatePack& residueTemplates();

// The same load from an explicit file, so a test can name the in-tree pack instead of depending
// on where the binary happens to live.
ResidueTemplatePack loadResidueTemplates(const std::filesystem::path& file);

ConnectivityResult connect(const Model& model, const ConnectivityOptions& options = {});

}  // namespace biocad::bio
