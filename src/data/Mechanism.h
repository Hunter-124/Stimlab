// data/Mechanism.h - JSON-serializable DTOs for mechanism of action, off-target
// panel coverage, pathway context, interaction flags and pharmacogenomic notes.
//
// HONESTY SCOPE, permanent:
//  - A mechanism is RETRIEVED with its reference, never inferred from a docking
//    pose or a fingerprint. There is no field for a predicted mechanism.
//  - The most prominent number in a panel-screen result is the count of targets
//    NOT screened, which is why `unscreened` exists and why `screened` alone is
//    never enough to render the view.
//  - There is no composite safety score and no cross-target score comparison:
//    different receptor preparations, box volumes and rotatable-bond penalties
//    make docking scores non-comparable across targets, so the DTO offers no
//    field in which to compare them.
//  - There is no pathway impact score. No database supports propagating a docking
//    score through a pathway graph, and such a number would be fabrication with a
//    scientific veneer.
//  - Interaction entries are a FLAG WITH A MECHANISM and a citation, never a
//    severity score. Pharmacogenomic entries are conditional notes, never a
//    genotype interpretation.
#pragma once

#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "data/Domain.h"

namespace biocad {

// One retrieved mechanism-of-action record. `actionType` comes from the source's
// own controlled vocabulary; `freeTextMechanism` is displayed AS TEXT and is never
// parsed into a vocabulary, because doing so invents structure the source does not
// have.
struct MechanismEntry {
    std::string        targetName;
    std::string        targetAccession;
    std::string        actionType;           // e.g. INHIBITOR, AGONIST
    std::string        freeTextMechanism;    // displayed verbatim, never parsed
    std::string        organism;
    std::vector<std::string> references;
    InhibitionModality modality = InhibitionModality::Unknown;   // BioCAD's own axis
    Provenance         provenance = Provenance::Measured;
    std::string        source;               // database name and release
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MechanismEntry, targetName, targetAccession, actionType,
                                   freeTextMechanism, organism, references, modality,
                                   provenance, source)

struct MechanismReport {
    std::string                 compoundId;
    std::vector<MechanismEntry> entries;
    bool                        retrievalAttempted = false;
    bool                        networkAvailable = false;
    // An empty entry list means "nothing was retrieved", which is a statement
    // about the query and the source, not about the compound.
    std::string                 coverageNote;
    std::vector<std::string>    warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MechanismReport, compoundId, entries, retrievalAttempted,
                                   networkAvailable, coverageNote, warnings)

// One target's result in a panel screen. There is no rank and no comparison field:
// see the header note on cross-target incomparability.
struct PanelTargetResult {
    std::string targetId;
    std::string targetName;
    bool        screened = false;
    std::string skipReason;        // why it was not screened, when it was not
    Quantity    affinity;          // whatever the docking module returned, with its tier
    std::string receptorPreparation;   // required: what was actually docked into
    std::string boxDefinition;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PanelTargetResult, targetId, targetName, screened,
                                   skipReason, affinity, receptorPreparation, boxDefinition)

// Panel coverage. `unscreened` is first in the struct because it is first in the
// view: a 44-target panel with 12 targets screened leaves the unknown fraction
// dominant, and the honest headline is the 32.
struct PanelScreenReport {
    std::string                    panelId;
    int                            unscreened = 0;
    int                            screened = 0;
    int                            panelSize = 0;
    std::vector<PanelTargetResult> results;
    // hERG margin only when the user supplied a MEASURED IC50; a predicted hERG
    // IC50 or any derived QT/TdP risk is prohibited, so there is no field for one.
    Quantity                       hergSafetyMargin;   // measured IC50 / free Cmax
    double                         hergMarginFlagBelow = 30.0;
    std::string                    coverageStatement;
    std::vector<std::string>       warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PanelScreenReport, panelId, unscreened, screened, panelSize,
                                   results, hergSafetyMargin, hergMarginFlagBelow,
                                   coverageStatement, warnings)

// Retrieved pathway membership, as a hierarchy. No score field, by design.
struct PathwayNode {
    std::string              stableId;
    std::string              name;
    std::string              species;
    std::vector<std::string> ancestorIds;
    std::string              url;         // deep link to the source
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PathwayNode, stableId, name, species, ancestorIds, url)

struct PathwayContext {
    std::string              accession;      // the UniProt accession queried
    std::vector<PathwayNode> pathways;
    std::string              source;         // database, licence and release
    bool                     networkAvailable = false;
    std::vector<std::string> warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PathwayContext, accession, pathways, source, networkAvailable,
                                   warnings)

// One interaction flag: a mechanism plus a citation. Deliberately no severity,
// no numeric risk and no recommendation - the mechanism IS the information.
struct InteractionFlag {
    std::string leftId;
    std::string rightId;
    std::string mechanism;        // "CYP1A2 inhibition", "MAO-A inhibition"
    std::string direction;        // which member is the perpetrator
    std::string evidence;         // what kind of evidence supports it
    std::string citation;
    std::string boundaryNote;     // the explicit not-a-recommendation sentence
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(InteractionFlag, leftId, rightId, mechanism, direction,
                                   evidence, citation, boundaryNote)

struct StackReport {
    std::vector<std::string>     members;      // compound and supplement ids entered
    std::vector<InteractionFlag> flags;
    std::vector<std::string>     unknownMembers;   // entered but not in any pack
    std::string                  coverageNote;
    std::vector<std::string>     warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(StackReport, members, flags, unknownMembers, coverageNote,
                                   warnings)

// A conditional pharmacogenomic note. Phenotype vocabulary only: UM, RM, NM, IM,
// PM - "extensive metabolizer" is deprecated and must not appear. Activity-score
// bands are carried as data so the panel does not restate them from memory. This
// is never a genotype interpretation and never a dosing statement.
struct PharmacogenomicNote {
    std::string gene;                 // CYP2D6, CYP2C19
    std::string phenotype;            // UM | RM | NM | IM | PM
    std::string activityScoreBand;    // e.g. "1.25 <= AS <= 2.25"
    std::string implication;          // the source's own conditional wording
    std::string citation;
    std::string source;               // database, licence and release
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PharmacogenomicNote, gene, phenotype, activityScoreBand,
                                   implication, citation, source)

struct PharmacogenomicReport {
    std::string                        compoundId;
    std::vector<PharmacogenomicNote>   notes;
    std::string                        boundaryStatement;
    bool                               networkAvailable = false;
    std::vector<std::string>           warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PharmacogenomicReport, compoundId, notes, boundaryStatement,
                                   networkAvailable, warnings)

}  // namespace biocad
