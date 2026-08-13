// tests/test_metabolites.cpp - the curated metabolite fact pack.
//
// The load-bearing assertion here is not that a lookup works: it is that a compound
// with NO curated fact comes back with an empty list AND a non-empty coverage note.
// That note is the only thing standing between this panel and the implication that
// BioCAD looked and found nothing - which would be a claim about the compound rather
// than about BioCAD's curation.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <filesystem>
#include <string>

#include "chem/Smiles.h"
#include "modules/Metabolites.h"
#include "packs/Pack.h"

using namespace biocad;

namespace {

std::filesystem::path factFile() {
    return std::filesystem::path(BIOCAD_ASSETS_DIR) / "packs" / "rules" / "metabolism-facts.json";
}

const MetaboliteFactPack& shippedPack() {
    static const MetaboliteFactPack p = loadMetaboliteFacts(factFile());
    return p;
}

Molecule mol(std::string id, std::string name = "") {
    Molecule m;
    m.id = std::move(id);
    m.name = name.empty() ? m.id : std::move(name);
    return m;
}

RealMetabolismFacts module_() { return RealMetabolismFacts(shippedPack()); }

bool hasMetabolite(const MetabolismReport& r, const std::string& needle) {
    return std::any_of(r.known.begin(), r.known.end(), [&](const MetaboliteFact& f) {
        return f.metaboliteName.find(needle) != std::string::npos;
    });
}

}  // namespace

TEST_CASE("the shipped metabolite fact pack loads without errors", "[metabolites]") {
    const auto& p = shippedPack();
    for (const auto& e : p.errors) WARN(e);
    REQUIRE(p.errors.empty());
    REQUIRE(p.schemaVersion == kMetaboliteFactSchemaVersion);
    REQUIRE(p.facts.size() >= 20);
}

TEST_CASE("every fact carries a citation and a responsible enzyme", "[metabolites]") {
    for (const auto& f : shippedPack().facts) {
        INFO(f.parentId << " -> " << f.metaboliteName);
        REQUIRE_FALSE(f.citation.empty());
        REQUIRE_FALSE(f.enzyme.empty());
        REQUIRE_FALSE(f.reaction.empty());
        REQUIRE_FALSE(f.significance.empty());
    }
}

TEST_CASE("every parentId resolves to a shipped library compound", "[metabolites]") {
    packs::LoadReport lib = packs::loadFrom(std::filesystem::path(BIOCAD_ASSETS_DIR) / "packs",
                                            /*builtin=*/true);
    REQUIRE(lib.errors.empty());
    const auto compounds = lib.compounds();
    REQUIRE_FALSE(compounds.empty());
    for (const auto& f : shippedPack().facts) {
        const bool found = std::any_of(compounds.begin(), compounds.end(),
                                       [&](const packs::PackCompound& c) {
                                           return c.id == f.parentId;
                                       });
        INFO("parentId '" << f.parentId << "' (metabolite " << f.metaboliteName
                          << ") is not a shipped compound id");
        REQUIRE(found);
    }
}

TEST_CASE("every authored metabolite SMILES parses", "[metabolites]") {
    for (const auto& f : shippedPack().facts) {
        if (f.metaboliteSmiles.empty()) continue;
        INFO(f.metaboliteName << ": " << f.metaboliteSmiles);
        REQUIRE(chem::parseSmiles(f.metaboliteSmiles).has_value());
    }
}

TEST_CASE("acetaminophen returns the NAPQI fact", "[metabolites]") {
    const auto r = module_().known(mol("acetaminophen", "Acetaminophen"));
    REQUIRE(hasMetabolite(r, "NAPQI"));
    const auto it = std::find_if(r.known.begin(), r.known.end(), [](const MetaboliteFact& f) {
        return f.metaboliteName.find("NAPQI") != std::string::npos;
    });
    REQUIRE(it != r.known.end());
    REQUIRE(it->enzyme.find("CYP2E1") != std::string::npos);
    REQUIRE(it->enzyme.find("CYP3A4") != std::string::npos);
    REQUIRE_FALSE(it->citation.empty());
}

TEST_CASE("N-demethylations of the shipped stimulants are on file", "[metabolites]") {
    REQUIRE(hasMetabolite(module_().known(mol("methamphetamine")), "Amphetamine"));
    REQUIRE(hasMetabolite(module_().known(mol("mdma")), "MDA"));
    REQUIRE(hasMetabolite(module_().known(mol("mdma")), "HMMA"));
}

TEST_CASE("a compound with no curated fact still carries the coverage note",
          "[metabolites]") {
    // Taurine is in the metabolic-and-nutrition pack and has no curated entry. The
    // empty list must arrive WITH the note, because the note is what stops the
    // absence being read as a finding.
    const auto r = module_().known(mol("taurine", "Taurine"));
    REQUIRE(r.known.empty());
    REQUIRE_FALSE(r.coverageNote.empty());
    REQUIRE(r.coverageNote.find("no curated") != std::string::npos);
    REQUIRE_FALSE(r.summary.empty());
}

TEST_CASE("a fact without a citation or an enzyme is rejected at load", "[metabolites]") {
    const auto p = parseMetaboliteFacts(R"({
      "schemaVersion": 1, "id": "t", "facts": [
        {"parentId":"x","metaboliteName":"Uncited","enzyme":"CYP3A4","reaction":"r",
         "significance":"s","citation":""},
        {"parentId":"x","metaboliteName":"NoEnzyme","enzyme":"","reaction":"r",
         "significance":"s","citation":"c"},
        {"parentId":"x","metaboliteName":"Good","enzyme":"CYP2D6","reaction":"r",
         "significance":"s","citation":"c","polymorphic":true}
      ]})");
    REQUIRE(p.facts.size() == 1);
    REQUIRE(p.facts[0].metaboliteName == "Good");
    REQUIRE(p.facts[0].polymorphic);
    REQUIRE(p.errors.size() == 2);
}

TEST_CASE("an unknown fact-pack schema version is an error, not a silent skip",
          "[metabolites]") {
    const auto p = parseMetaboliteFacts(R"({"schemaVersion": 99, "facts": []})");
    REQUIRE(p.facts.empty());
    REQUIRE(p.errors.size() == 1);
}
