// modules/MechanismModule.h - retrieved mechanism of action, off-target panel
// coverage, pathway context, interaction flags and pharmacogenomic notes.
//
// Every method here RETRIEVES. Nothing in this file derives a mechanism from a
// structure, a fingerprint or a docking pose, and there is no entry point that
// would let it: the honest half of "what does this compound do" is a citation, and
// the dishonest half is a plausible sentence with no source behind it.
//
// LICENCE DISCIPLINE IS PART OF THE IMPLEMENTATION, not a note in the docs. The
// only sources this module queries or bundles are:
//   ChEMBL (CC BY-SA 3.0 - queried live, cached under cache/api, and any committed
//           ChEMBL-derived pack carries the share-alike notice)
//   Reactome (CC0), UniProt (CC BY 4.0), RCSB, AlphaFold DB (CC BY 4.0),
//   openFDA (US public domain), CPIC (CC0), and the FDA drug-interaction
//   labeling tables (US public domain).
// KEGG, STRING, DrugBank, PDBbind, BioLiP and Guide to PHARMACOLOGY are NOT
// queried and NOT bundled. KEGG appears only as a deep link to a kegg.jp pathway
// page, because linking is not redistribution. docs/mechanism.md carries the table
// with each verdict and its reason.
//
// With the science feature off there is no HTTP transport compiled in. Every
// network-backed method then returns its DTO with networkAvailable = false and the
// reason in `warnings` - never a build error, and never invented content.
#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "contracts/IModules.h"
#include "data/Mechanism.h"
#include "packs/Pack.h"

namespace biocad {

// The only schema version this build understands, for every mechanism pack kind.
// An unknown version is an error in `errors`, never a silent skip.
inline constexpr int kMechanismPackSchemaVersion = 1;

// ---------------------------------------------------------------- panel packs
// One row of a panel roster. `targetRef` is the receptor-target id to look for in
// the loaded target packs; empty means the panel names a target BioCAD has no
// receptor for at all, which is a coverage gap declared before docking is even
// attempted.
struct MechanismPanelTarget {
    std::string id;
    std::string name;
    std::string targetClass;
    std::string targetRef;
};

// `declaredSize` is the size the PANEL claims, which can exceed the number of rows
// this file enumerates (SafetyScreen87 declares 87 and enumerates 44). Coverage is
// always measured against the declared size, so an incomplete roster inflates the
// unscreened count instead of quietly shrinking the denominator.
struct MechanismPanelPack {
    std::string panelId;
    std::string title;
    int         declaredSize = 0;
    std::string licence;
    std::string sourceNote;
    std::string coverageNote;
    std::vector<MechanismPanelTarget> targets;
};

// ----------------------------------------------------- action-type vocabulary
struct ActionTypeValue {
    std::string value;   // the source's own token, e.g. "INHIBITOR"
    std::string gloss;   // BioCAD's editorial wording, never the source's text
};

struct ActionTypePack {
    std::string id;
    std::string title;
    int         declaredSize = 0;   // values the source's release documents
    std::string licence;
    std::string sourceNote;
    std::vector<ActionTypeValue> values;

    // Case-sensitive membership. An action_type that is not here is reported
    // VERBATIM and flagged, never mapped onto a neighbouring value.
    [[nodiscard]] bool recognises(const std::string& actionType) const;
};

// ------------------------------------------------------------ interaction pack
// A member's enzyme/transporter roles and mechanism classes. There is no severity
// and no numeric risk field: the mechanism is the information.
struct StackMember {
    std::string id;
    std::string name;
    std::string kind;                       // drug | supplement | food | lifestyle
    std::vector<std::string> substrateOf;   // CYP names
    std::vector<std::string> inhibitorOf;
    std::vector<std::string> inducerOf;
    std::vector<std::string> mechanisms;    // free-text mechanism classes
    std::string citation;
    std::string note;
};

// Pairs two mechanism classes into one flag, e.g. "MAO inhibition" x "monoamine
// release". Both orderings are tried; the left class names the perpetrator.
struct StackClassRule {
    std::string left;
    std::string right;
    std::string mechanism;
    std::string evidence;
    std::string citation;
};

struct InteractionPack {
    std::string id;
    std::string title;
    std::string licence;
    std::string boundaryNote;   // copied onto every flag this pack produces
    std::string coverageNote;
    std::vector<StackMember>    members;
    std::vector<StackClassRule> classRules;

    // Case-insensitive lookup by id or display name. nullptr when the member is
    // unknown, which is what puts it in StackReport::unknownMembers.
    [[nodiscard]] const StackMember* find(const std::string& idOrName) const;
};

// ------------------------------------------------------- pharmacogenomics pack
struct PgxPhenotype {
    std::string code;   // UM | RM | NM | IM | PM
    std::string term;   // the standardised term; "extensive metabolizer" is deprecated
};

struct PgxActivityBand {
    std::string gene;
    std::string phenotype;
    std::string band;   // e.g. "1.25 <= AS <= 2.25"
};

struct PgxNote {
    std::string compound;
    std::string gene;
    std::string phenotype;
    std::string implication;
};

struct PharmacogenomicsPack {
    std::string id;
    std::string title;
    std::string licence;
    std::string boundaryStatement;
    std::vector<PgxPhenotype>    phenotypes;
    std::vector<PgxActivityBand> bands;
    std::vector<PgxNote>         notes;
};

// Everything loaded from assets/packs/mechanism. `errors` is surfaced by the panels
// rather than swallowed: a pack that silently failed to parse is indistinguishable
// from a compound genuinely having no retrievable mechanism, which is exactly the
// confusion this module exists to prevent.
struct MechanismPacks {
    std::vector<MechanismPanelPack> panels;
    ActionTypePack                  actionTypes;
    InteractionPack                 interactions;
    PharmacogenomicsPack            pharmacogenomics;
    std::vector<std::string>        errors;
    std::string                     sourceDir;

    [[nodiscard]] const MechanismPanelPack* panel(const std::string& panelId) const;
};

// assets/packs/mechanism, resolved relative to the built-in pack root. Empty when
// the pack root itself could not be located.
std::filesystem::path defaultMechanismPackDir();

// Load every *.json in `dir`, dispatching on the document's "kind". A missing
// directory yields empty packs and one error; nothing here throws.
MechanismPacks loadMechanismPacks(const std::filesystem::path& dir);

// Parse one document into `into`, appending to `into.errors` on failure. Exposed so
// tests can drive the parser from a string without touching the filesystem.
void parseMechanismDocument(const std::string& text, const std::string& sourcePath,
                            MechanismPacks& into);

// KEGG is NOT a data source here: its REST API is licensed for academic use only, so
// this build never calls it and bundles none of its content. A DEEP LINK to a public
// kegg.jp pathway page is a different act - linking is not redistribution - so the
// Pathways panel offers one for the user to open themselves.
// `mapId` is a KEGG pathway map id such as "hsa04726"; the result is a page URL and
// is never fetched by BioCAD.
std::string keggPathwayUrl(const std::string& mapId);
const char* keggDeepLinkNote();

// True iff this build has an HTTP transport compiled in (the science feature).
bool mechanismNetworkAvailable();

// The one sentence every offline report carries, so the reason is identical
// everywhere instead of being re-worded per call site.
const char* mechanismOfflineReason();

// The measured hERG inputs. Both are USER-SUPPLIED measurements: a predicted hERG
// IC50 is prohibited, so there is no field for one and no code path that fills
// these in from structure. Units are mol/L.
struct HergInput {
    double      measuredIc50Molar = -1.0;   // < 0 = absent
    double      freeCmaxMolar     = -1.0;   // < 0 = absent
    std::string citation;                   // where the measured IC50 came from
};

// The one implementation of IMechanismModule.
class RealMechanism final : public IMechanismModule {
public:
    // `docking` is the EXISTING docking module - a panel screen runs it, it does not
    // reimplement scoring. `catalog` is the loaded compound catalog, used only to
    // resolve a compound id to its ChEMBL cross-reference.
    RealMechanism(const IDockingModule* docking, std::vector<packs::PackCompound> catalog);
    RealMechanism(const IDockingModule* docking, std::vector<packs::PackCompound> catalog,
                  MechanismPacks packs);

    MechanismReport       mechanisms(const std::string& compoundId) const override;
    PanelScreenReport     screenPanel(const Molecule& m, const std::string& panelId) const override;
    PathwayContext        pathways(const std::string& uniprotAccession) const override;
    StackReport           checkStack(const std::vector<std::string>& memberIds) const override;
    PharmacogenomicReport pharmacogenomics(const std::string& compoundId) const override;

    // The hERG margin is computed only from these. Setting an IC50 without a free
    // Cmax (or the reverse) leaves the margin NotComputed naming the missing one.
    void setHergInput(HergInput in) { herg_ = std::move(in); }
    [[nodiscard]] const HergInput& hergInput() const { return herg_; }

    [[nodiscard]] const MechanismPacks& packs() const { return packs_; }

private:
    const IDockingModule*            docking_ = nullptr;
    std::vector<packs::PackCompound> catalog_;
    MechanismPacks                   packs_;
    HergInput                        herg_;
};

}  // namespace biocad
