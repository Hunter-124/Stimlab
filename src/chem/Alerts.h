// chem/Alerts.h - structural alerts for metabolic bioactivation, as data.
//
// WHAT AN ALERT IS: a substructure that the medicinal-chemistry literature has
// associated with reactive-metabolite formation. A match is a LIABILITY FLAG and
// nothing more. It is not a toxicity verdict, not a prediction that
// bioactivation happens in vivo, and not a reason to reject a compound: whether
// a route is taken depends on the enzymes present, the competing clearance
// routes, the dose and the detoxication capacity, none of which a substructure
// knows. Several of the alerts below match widely used, marketed drugs - which
// is exactly why the flag is worded as a route, never as an outcome.
//
// WHY THE RULES ARE A PACK: every published alert set is expressed in SMARTS, so
// with a matcher in the tree the rules are data (assets/packs/rules/
// alerts-bioactivation.json) rather than hand-coded C++ conditions nobody can
// cite. Editing an alert, or adding a separately-licensed alert pack later, is
// then a data change.
//
// LICENSING: PAINS, Brenk and the ChEMBL-derived alert sets are deliberately not
// shipped in the built-in pack (BSD-3 / CC BY-SA 3.0 obligations that need their
// own pack and NOTICE entry). The built-in alerts are authored in-house from the
// bioactivation literature, each carrying its citation.
#pragma once

#include <string>
#include <vector>

#include "chem/Molecule.h"

namespace biocad::chem {

// One matched alert. `atoms` are molecule atom indices, so the UI can highlight
// exactly the substructure that raised the flag rather than the whole molecule.
struct AlertHit {
    std::string      key;        // stable pack id, e.g. "para-aminophenol"
    std::string      label;
    std::string      mechanism;  // the metabolic route, stated as a route
    std::string      citation;
    std::vector<int> atoms;
    bool             warn = false;  // pack severity: warn vs info. Never "danger".
};

// Screens `mol` against the built-in alert pack. Rings and aromaticity are
// perceived on a copy, so the caller's molecule is never mutated and uppercase
// SMILES input is still matched by the aromatic patterns.
//
// An empty result means "no alert in this pack matched". That is not a safety
// claim: the pack is a short, in-house list, and absence of a listed motif says
// nothing about the motifs it does not list.
std::vector<AlertHit> screenAlerts(const Molecule& mol);

// Problems found while loading the alert pack: a missing file, a malformed
// document, or a rule whose SMARTS does not parse (named, with the parser's
// message). A rule that cannot be parsed is dropped and reported here, never
// silently treated as a rule that matches nothing.
const std::vector<std::string>& alertPackErrors();

// Number of alerts the pack loaded. Zero with a non-empty alertPackErrors()
// means the screen is inoperative, which the UI must say out loud.
std::size_t alertRuleCount();

}  // namespace biocad::chem
