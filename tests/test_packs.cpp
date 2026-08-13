// tests/test_packs.cpp - the catalog is data, and the loader is honest about it.
//
// The rules under test are the ones that keep a broken pack from looking like a
// working application: an unknown schema version is a named error, a duplicate id
// is reported rather than silently resolved, and a binding-site box may never
// exist without a structure to justify it.
#include <catch2/catch_test_macros.hpp>

#include <string>

#include "core/Error.h"
#include "packs/Pack.h"

using namespace biocad;

namespace {

// A minimal, valid pack. Tests mutate a copy of this text rather than sharing
// state through a fixture file.
std::string packText(const std::string& id, const std::string& compoundId,
                     const std::string& targetId, int schemaVersion = 1) {
    return std::string(R"PK({"schemaVersion": )PK") + std::to_string(schemaVersion) +
           R"PK(, "id": ")PK" + id +
           R"PK(", "title": "T", "description": "d",
                "compounds": [{"id": ")PK" + compoundId +
           R"PK(", "name": "Caffeine",
                              "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                              "drugClass": "Xanthine stimulant",
                              "legalUs": "Unscheduled",
                              "notes": "n",
                              "xrefs": {"chembl": "CHEMBL113", "pubchemCid": 2519}}],
                "targets": [{"id": ")PK" + targetId +
           R"PK(", "name": "DAT (dopamine transporter)",
                            "pdb": "4M48", "uniprot": "Q01959", "headline": true,
                            "panels": ["safetyscreen44"],
                            "box": {"cx": -2.1, "cy": 12.4, "cz": -6.8,
                                    "sx": 20.0, "sy": 20.0, "sz": 20.0}}]})PK";
}

}  // namespace

TEST_CASE("A pack round-trips through parse with its xrefs and metadata", "[packs]") {
    const auto p = packs::parseString(packText("a", "caffeine", "DAT"), "a.json");
    REQUIRE(p.schemaVersion == packs::kSchemaVersion);
    REQUIRE(p.id == "a");
    REQUIRE(p.sourcePath == "a.json");
    REQUIRE(p.compounds.size() == 1);
    REQUIRE(p.targets.size() == 1);

    const auto& c = p.compounds.front();
    REQUIRE(c.id == "caffeine");
    REQUIRE(c.xrefs.chembl == "CHEMBL113");
    REQUIRE(c.xrefs.pubchemCid == 2519);

    // Identity and metadata only: the numeric fields stay at zero because the
    // chem engine is the single source of every descriptor.
    const Molecule m = c.molecule();
    REQUIRE(m.smiles == "CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
    REQUIRE(m.formula.empty());
    REQUIRE(m.molWeight == 0.0);
    REQUIRE(m.logP == 0.0);
    REQUIRE(m.tpsa == 0.0);
    REQUIRE(m.hbd == 0);
    REQUIRE(m.hba == 0);
    REQUIRE(m.rotatableBonds == 0);
    REQUIRE(m.drugClass == "Xanthine stimulant");

    const auto& t = p.targets.front();
    REQUIRE(t.hasBox);
    REQUIRE(t.headline);
    REQUIRE(t.uniprot == "Q01959");
    REQUIRE(t.panels.size() == 1);
    REQUIRE(t.target.box.cx == -2.1);
}

TEST_CASE("A pack that authors descriptor numbers is rejected by name", "[packs]") {
    // Two sources of truth for a molecular weight would silently diverge, with the
    // authored one always losing. The author has to find out, so this is loud.
    bool threw = false;
    try {
        packs::parseString(R"PK({"schemaVersion":1,"id":"a","title":"T","compounds":[
            {"id":"caffeine","name":"Caffeine","smiles":"CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
             "properties":{"molWeight":194.19}}]})PK",
                           "authored.json");
    } catch (const Error& e) {
        threw = true;
        REQUIRE(e.code == Error::Code::Parse);
        REQUIRE(e.message.find("caffeine") != std::string::npos);
        REQUIRE(e.message.find("properties") != std::string::npos);
        REQUIRE(e.message.find("SMILES") != std::string::npos);
    }
    REQUIRE(threw);
}

TEST_CASE("An unknown schemaVersion is a named error, never a silent skip", "[packs]") {
    bool threw = false;
    try {
        packs::parseString(packText("a", "caffeine", "DAT", /*schemaVersion=*/99), "future.json");
    } catch (const Error& e) {
        threw = true;
        REQUIRE(e.code == Error::Code::Unsupported);
        REQUIRE(e.message.find("99") != std::string::npos);
        REQUIRE(e.message.find("future.json") != std::string::npos);
    }
    REQUIRE(threw);

    // A pack with no schemaVersion at all is equally a hard error: an unversioned
    // document cannot be validated against anything.
    REQUIRE_THROWS_AS(packs::parseString(R"PK({"id":"a","title":"T"})PK", "old.json"), Error);
}

TEST_CASE("A binding-site box requires a structure to justify it", "[packs]") {
    // Box without a pdb: an unverifiable binding site, rejected outright.
    REQUIRE_THROWS_AS(
        packs::parseString(R"PK({"schemaVersion":1,"id":"a","title":"T","targets":[
            {"id":"X","name":"X","box":{"cx":0,"cy":0,"cz":0,"sx":10,"sy":10,"sz":10}}]})PK",
                           "nobox.json"),
        Error);

    // No box at all is fine - that is an honest coverage gap, and it is not dockable.
    const auto p = packs::parseString(
        R"PK({"schemaVersion":1,"id":"a","title":"T","targets":[{"id":"X","name":"X"}]})PK",
        "gap.json");
    REQUIRE(p.targets.size() == 1);
    REQUIRE_FALSE(p.targets.front().hasBox);

    // A degenerate box is rejected rather than clamped into something plausible.
    REQUIRE_THROWS_AS(
        packs::parseString(R"PK({"schemaVersion":1,"id":"a","title":"T","targets":[
            {"id":"X","name":"X","pdb":"1ABC","box":{"cx":0,"cy":0,"cz":0,"sx":0,"sy":10,"sz":10}}]})PK",
                          "flat.json"),
        Error);
}

TEST_CASE("Required fields are required", "[packs]") {
    REQUIRE_THROWS_AS(packs::parseString(R"PK({"schemaVersion":1,"title":"T"})PK", "x.json"), Error);
    REQUIRE_THROWS_AS(
        packs::parseString(R"PK({"schemaVersion":1,"id":"a","title":"T",
                              "compounds":[{"id":"c","name":"C"}]})PK", "x.json"),
        Error);  // a compound without SMILES cannot be computed on
    REQUIRE_THROWS_AS(packs::parseString("not json at all", "x.json"), Error);
}

TEST_CASE("A duplicate id inside one pack is rejected", "[packs]") {
    REQUIRE_THROWS_AS(
        packs::parseString(R"PK({"schemaVersion":1,"id":"a","title":"T","compounds":[
            {"id":"c","name":"C","smiles":"C"},{"id":"c","name":"C2","smiles":"CC"}]})PK",
                          "dup.json"),
        Error);
}

TEST_CASE("A later pack overrides an earlier one by pack id", "[packs]") {
    packs::LoadReport report;
    report.packs.push_back(packs::parseString(packText("shared", "caffeine", "DAT"), "builtin.json"));

    // Same pack id, different contents: the user copy wins wholesale, so the user
    // can always name the one document that produced a row.
    report.packs.push_back(
        packs::parseString(packText("shared", "theobromine", "SERT"), "user.json"));

    const auto compounds = report.compounds();
    REQUIRE(compounds.size() == 1);
    REQUIRE(compounds.front().id == "theobromine");

    const auto targets = report.targets();
    REQUIRE(targets.size() == 1);
    REQUIRE(targets.front().target.id == "SERT");
    REQUIRE(report.errors.empty());  // overriding by pack id is not an error
}

TEST_CASE("A duplicate id across different packs is reported, not resolved", "[packs]") {
    packs::LoadReport report;
    report.packs.push_back(packs::parseString(packText("one", "caffeine", "DAT"), "one.json"));
    report.packs.push_back(packs::parseString(packText("two", "caffeine", "DAT"), "two.json"));

    const auto compounds = report.compounds();
    REQUIRE(compounds.size() == 1);          // the later definition is dropped
    REQUIRE_FALSE(report.errors.empty());    // and the user is told why
    REQUIRE(report.errors.front().find("caffeine") != std::string::npos);
    REQUIRE(report.errors.front().find("one") != std::string::npos);
    REQUIRE(report.errors.front().find("two") != std::string::npos);
}

TEST_CASE("findTarget resolves by id, by display name and by prefix", "[packs]") {
    packs::LoadReport report;
    report.packs.push_back(packs::parseString(packText("one", "caffeine", "DAT"), "one.json"));

    REQUIRE(report.findTarget("DAT").has_value());
    REQUIRE(report.findTarget("dat").has_value());
    REQUIRE(report.findTarget("DAT (dopamine transporter)PK").has_value());
    REQUIRE_FALSE(report.findTarget("__absent__").has_value());
}

TEST_CASE("A missing pack directory yields an empty report, not an error", "[packs]") {
    const auto report = packs::loadFrom("__no_such_pack_directory__");
    REQUIRE(report.packs.empty());
    REQUIRE(report.errors.empty());
}
