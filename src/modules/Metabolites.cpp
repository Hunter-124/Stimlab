#include "modules/Metabolites.h"

#include <fstream>

#include <nlohmann/json.hpp>

#include "packs/Pack.h"

namespace biocad {
namespace {

std::string requiredString(const nlohmann::json& j, const char* key) {
    if (!j.contains(key) || !j.at(key).is_string()) return {};
    return j.at(key).get<std::string>();
}

}  // namespace

std::filesystem::path defaultMetaboliteFactPath() {
    const auto root = packs::builtinPackDir();
    if (root.empty()) return {};
    return root / "rules" / "metabolism-facts.json";
}

MetaboliteFactPack parseMetaboliteFacts(const std::string& text, std::string sourcePath) {
    MetaboliteFactPack pack;
    pack.sourcePath = std::move(sourcePath);

    nlohmann::json j;
    try {
        j = nlohmann::json::parse(text);
    } catch (const std::exception& e) {
        pack.errors.push_back(pack.sourcePath + ": not valid JSON - " + e.what());
        return pack;
    }
    if (!j.is_object()) {
        pack.errors.push_back(pack.sourcePath + ": top level is not a JSON object.");
        return pack;
    }

    pack.schemaVersion = j.value("schemaVersion", 0);
    if (pack.schemaVersion != kMetaboliteFactSchemaVersion) {
        pack.errors.push_back(pack.sourcePath + ": schemaVersion " +
                              std::to_string(pack.schemaVersion) + " is not understood by this "
                              "build (expected " +
                              std::to_string(kMetaboliteFactSchemaVersion) + ").");
        return pack;
    }
    pack.id          = requiredString(j, "id");
    pack.title       = requiredString(j, "title");
    pack.description = requiredString(j, "description");

    if (!j.contains("facts") || !j.at("facts").is_array()) {
        pack.errors.push_back(pack.sourcePath + ": missing the \"facts\" array.");
        return pack;
    }

    std::size_t index = 0;
    for (const auto& e : j.at("facts")) {
        const std::string where = pack.sourcePath + " fact #" + std::to_string(index++);
        if (!e.is_object()) {
            pack.errors.push_back(where + ": not an object.");
            continue;
        }
        MetaboliteFact f;
        f.parentId         = requiredString(e, "parentId");
        f.metaboliteName   = requiredString(e, "metaboliteName");
        f.metaboliteSmiles = requiredString(e, "metaboliteSmiles");
        f.enzyme           = requiredString(e, "enzyme");
        f.reaction         = requiredString(e, "reaction");
        f.significance     = requiredString(e, "significance");
        f.citation         = requiredString(e, "citation");
        f.polymorphic      = e.value("polymorphic", false);

        // A fact without a citation is not a fact, and a transformation without the
        // enzyme that performs it cannot be reasoned about (a CYP route and an
        // esterase route have completely different interaction consequences). Both
        // are rejected at load rather than rendered as a blank cell.
        if (f.parentId.empty() || f.metaboliteName.empty()) {
            pack.errors.push_back(where + ": parentId and metaboliteName are required.");
            continue;
        }
        if (f.citation.empty()) {
            pack.errors.push_back(where + " (" + f.metaboliteName +
                                  "): no citation - an uncited transformation is not a fact.");
            continue;
        }
        if (f.enzyme.empty()) {
            pack.errors.push_back(where + " (" + f.metaboliteName +
                                  "): no enzyme - the responsible enzyme is required.");
            continue;
        }
        pack.facts.push_back(std::move(f));
    }
    return pack;
}

MetaboliteFactPack loadMetaboliteFacts(const std::filesystem::path& file) {
    MetaboliteFactPack pack;
    pack.sourcePath = file.string();
    if (file.empty()) {
        pack.errors.push_back("no pack root located, so assets/packs/rules/metabolism-facts.json "
                              "could not be read.");
        return pack;
    }
    std::error_code ec;
    if (!std::filesystem::is_regular_file(file, ec)) {
        pack.errors.push_back(pack.sourcePath + ": file not found.");
        return pack;
    }
    std::ifstream in(file, std::ios::binary);
    if (!in) {
        pack.errors.push_back(pack.sourcePath + ": could not be opened for reading.");
        return pack;
    }
    const std::string text{std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
    return parseMetaboliteFacts(text, file.string());
}

std::string metaboliteCoverageNote(bool anyFacts) {
    // One sentence, identical in both branches, because the honesty guarantee is
    // that the reader cannot mistake "we have nothing curated" for "there is
    // nothing to find". The count differs; the disclaimer does not.
    const std::string always =
        "An empty list means BioCAD has no curated metabolite fact for this compound. It does "
        "NOT mean the compound has no metabolites, and it is not evidence that the compound is "
        "metabolically stable: absence of a curated entry is absence of curation only.";
    if (anyFacts) {
        return "This list is the curated coverage BioCAD ships for this compound, not its "
               "complete metabolic fate; minor and unstudied routes are simply absent. " + always;
    }
    return always;
}

const char* metaboliteNoEnumerationNote() {
    return "BioCAD does not enumerate hypothetical metabolites on this surface. An independent "
           "EPA cross-tool benchmark (Boyce et al. 2022, Computational Toxicology 21:100208) "
           "measured rule-based biotransformation predictors - SyGMa, Meteor, BioTransformer, "
           "TIMES, the OECD QSAR Toolbox and CTS - at 1.1-29% precision and 14.7-28.3% "
           "sensitivity. Placing output at that accuracy beside a cited transformation would "
           "lend it credibility it has not earned, so only curated facts appear here.";
}

RealMetabolismFacts::RealMetabolismFacts()
    : pack_(loadMetaboliteFacts(defaultMetaboliteFactPath())) {}

RealMetabolismFacts::RealMetabolismFacts(MetaboliteFactPack pack) : pack_(std::move(pack)) {}

MetabolismReport RealMetabolismFacts::known(const Molecule& m) const {
    MetabolismReport r;
    r.moleculeId = m.id;
    for (const auto& f : pack_.facts) {
        if (f.parentId == m.id) r.known.push_back(f);
    }
    r.coverageNote = metaboliteCoverageNote(!r.known.empty());

    if (!pack_.errors.empty()) {
        // A load failure must not read as "no metabolites". Say it in the summary,
        // where the user is looking, before anything else.
        r.summary = "The curated metabolite fact pack did not load cleanly (" +
                    std::to_string(pack_.errors.size()) +
                    " error(s)); this list may be incomplete for that reason alone. ";
    }
    if (r.known.empty()) {
        r.summary += m.name.empty() ? "No curated metabolite fact is on file."
                                    : "No curated metabolite fact is on file for " + m.name + ".";
        return r;
    }

    std::size_t polymorphic = 0;
    for (const auto& f : r.known) polymorphic += f.polymorphic ? 1 : 0;
    r.summary += std::to_string(r.known.size()) +
                 (r.known.size() == 1 ? " curated transformation" : " curated transformations") +
                 " on file, each with a citation.";
    if (polymorphic > 0) {
        r.summary += " " + std::to_string(polymorphic) +
                     " of them run on a genetically polymorphic enzyme, so the exposure differs "
                     "between phenotypes rather than being a single number.";
    }
    return r;
}

}  // namespace biocad
