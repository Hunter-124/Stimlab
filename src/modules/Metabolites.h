// modules/Metabolites.h - the curated, cited metabolite fact pack.
//
// This is the honest and cheap half of metabolite prediction, and it deliberately
// ships before any rule-based enumeration exists. The independent EPA cross-tool
// benchmark of rule-based biotransformation predictors (Boyce et al. 2022,
// Computational Toxicology 21:100208) measured precision of 1.1-29% and sensitivity
// of 14.7-28.3% across SyGMa, Meteor, BioTransformer, TIMES, the OECD QSAR Toolbox
// and CTS. A single transformation that has actually been characterised, with the
// reference attached, therefore carries more information than a hundred enumerated
// candidates - and unlike them it can carry Provenance::Measured.
//
// SAFETY SCOPE: this module retrieves facts. It does not enumerate, score, rank or
// predict, and it emits no quantitative fraction of dose, clearance contribution or
// toxicity verdict. There is no entry point that would let it.
#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "contracts/IModules.h"
#include "data/Domain.h"

namespace biocad {

// One parsed fact document. `errors` is everything that failed to load and is
// surfaced rather than swallowed: a fact pack that silently failed to parse is
// indistinguishable from a compound genuinely having no curated metabolites, and
// that is exactly the confusion this whole module exists to prevent.
struct MetaboliteFactPack {
    int         schemaVersion = 0;
    std::string id;
    std::string title;
    std::string description;
    std::string sourcePath;
    std::vector<MetaboliteFact>  facts;
    std::vector<std::string>     errors;
};

// The only schema version this build understands. An unknown version is an error
// in `errors`, never a silent skip.
inline constexpr int kMetaboliteFactSchemaVersion = 1;

// assets/packs/rules/metabolism-facts.json, resolved relative to the pack root.
// Empty when the pack root itself could not be located.
std::filesystem::path defaultMetaboliteFactPath();

// Parse one fact document from disk. A missing or malformed file yields a pack with
// no facts and a populated `errors`; it never throws.
MetaboliteFactPack loadMetaboliteFacts(const std::filesystem::path& file);

// Parse from a JSON string, for tests and embedded literals.
MetaboliteFactPack parseMetaboliteFacts(const std::string& text,
                                        std::string sourcePath = "<embedded>");

// The coverage note. It is a function rather than a string literal at each call
// site because every report must carry the SAME sentence: an empty fact list means
// BioCAD has no curated fact, not that the compound has no metabolites.
std::string metaboliteCoverageNote(bool anyFacts);

// The fixed statement of why no hypothetical metabolites are enumerated here,
// with the Boyce 2022 numbers that justify it. Rendered by the panel verbatim.
const char* metaboliteNoEnumerationNote();

class RealMetabolismFacts final : public IMetabolismFactsModule {
public:
    RealMetabolismFacts();
    explicit RealMetabolismFacts(MetaboliteFactPack pack);

    MetabolismReport known(const Molecule& m) const override;

    [[nodiscard]] const MetaboliteFactPack& pack() const { return pack_; }

private:
    MetaboliteFactPack pack_;
};

}  // namespace biocad
