// tests/test_chem_crippen.cpp - Wildman-Crippen logP / MR.
//
// What these cases defend, and why each one is here:
//  * every heavy atom of every catalog compound gets a class - an unclassified
//    atom silently drops its contribution, which is a wrong number, not a
//    missing one;
//  * Kekule and aromatic spellings of one molecule give one answer - the whole
//    reason crippen() perceives aromaticity itself;
//  * the sum is the table's sum, checked against values computed by hand from
//    assets/packs/descriptors/crippen.json;
//  * MR is positive and grows with the homologous series;
//  * a Heuristic-free honest failure: a missing pack yields ok == false.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <string>

#include "chem/Crippen.h"
#include "chem/Descriptors.h"
#include "chem/Smiles.h"
#include "packs/Pack.h"

using namespace biocad;
using Catch::Matchers::WithinAbs;

namespace {

// The descriptor pack is data on disk. ctest runs from the build directory, so
// point the loader at the in-tree copy. chem::crippen() caches only a SUCCESSFUL
// load, so calling this from every case is enough no matter which test ran first.
void useInTreePack() {
    const std::string dir = std::string(BIOCAD_ASSETS_DIR) + "/packs/descriptors";
#if defined(_WIN32)
    _putenv_s("BIOCAD_DESCRIPTOR_DIR", dir.c_str());
#else
    setenv("BIOCAD_DESCRIPTOR_DIR", dir.c_str(), 1);
#endif
}

chem::Molecule parse(const std::string& smiles) {
    auto m = chem::parseSmiles(smiles);
    REQUIRE(m.has_value());
    return *m;
}

}  // namespace

TEST_CASE("Crippen classifies every heavy atom of every catalog compound", "[chem][crippen]") {
    useInTreePack();
    auto report = packs::loadFrom(std::filesystem::path(BIOCAD_ASSETS_DIR) / "packs");
    REQUIRE(report.errors.empty());
    const auto compounds = report.compounds();
    REQUIRE(compounds.size() >= 60);

    for (const auto& c : compounds) {
        if (c.smiles.empty()) continue;
        const auto r = chem::crippen(parse(c.smiles));
        INFO(c.id << "  " << c.smiles << "  " << r.note);
        CHECK(r.unclassified.empty());
        CHECK(r.ok);
        CHECK(r.atomTypes.size() == static_cast<std::size_t>(chem::heavyAtomCount(parse(c.smiles))));
    }
}

TEST_CASE("Crippen logP is independent of Kekule vs aromatic input", "[chem][crippen]") {
    useInTreePack();
    struct Pair { const char* aromatic; const char* kekule; };
    const Pair pairs[] = {
        {"Cn1cnc2c1c(=O)n(C)c(=O)n2C", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"},  // caffeine
        {"c1ccccc1", "C1=CC=CC=C1"},                                     // benzene
        {"CN1CCCC1c1cccnc1", "CN1CCCC1C1=CC=CN=C1"},                     // nicotine
        {"CC(=O)Oc1ccccc1C(=O)O", "CC(=O)OC1=CC=CC=C1C(=O)O"},           // aspirin
        {"c1ccoc1", "C1=COC=C1"},                                        // furan
        {"c1ccsc1", "C1=CSC=C1"},                                        // thiophene
    };
    for (const auto& p : pairs) {
        const auto a = chem::crippen(parse(p.aromatic));
        const auto k = chem::crippen(parse(p.kekule));
        INFO(p.aromatic << " vs " << p.kekule);
        REQUIRE(a.ok);
        REQUIRE(k.ok);
        CHECK_THAT(k.logP, WithinAbs(a.logP, 1e-12));
        CHECK_THAT(k.molarRefractivity, WithinAbs(a.molarRefractivity, 1e-12));
    }
}

TEST_CASE("Crippen sums the published contributions exactly", "[chem][crippen]") {
    useInTreePack();
    // Benzene: 6 x C18 (0.1581) + 6 x H1 (0.1230) = 1.6866, and
    //          6 x 3.350   + 6 x 1.057            = 26.442.
    const auto benzene = chem::crippen(parse("c1ccccc1"));
    REQUIRE(benzene.ok);
    CHECK_THAT(benzene.logP, WithinAbs(6 * 0.1581 + 6 * 0.1230, 1e-9));
    CHECK_THAT(benzene.molarRefractivity, WithinAbs(6 * 3.350 + 6 * 1.057, 1e-9));

    // Ethanol: C1 (0.1441) + C3 (-0.2035) + O2 (-0.2893) + 5 x H1 + 1 x H2.
    const auto ethanol = chem::crippen(parse("CCO"));
    REQUIRE(ethanol.ok);
    CHECK(ethanol.atomTypes == std::vector<std::string>{"C1", "C3", "O2"});
    CHECK_THAT(ethanol.logP,
               WithinAbs(0.1441 - 0.2035 - 0.2893 + 5 * 0.1230 + (-0.2677), 1e-9));

    // A carboxylic acid OH is H4 (0.2980), not the H2 an alcohol would get:
    // acetic acid = C1 + C5 + O9 + O2 + 3 x H1 + 1 x H4.
    const auto acetic = chem::crippen(parse("CC(=O)O"));
    REQUIRE(acetic.ok);
    CHECK(acetic.atomTypes == std::vector<std::string>{"C1", "C5", "O9", "O2"});
    CHECK_THAT(acetic.logP,
               WithinAbs(0.1441 - 0.2783 - 0.1526 - 0.2893 + 3 * 0.1230 + 0.2980, 1e-9));

    // Phenol's OH is H2, because the pack's H2 rule (O on aromatic carbon) is
    // tested before the acid rule.
    const auto phenol = chem::crippen(parse("Oc1ccccc1"));
    REQUIRE(phenol.ok);
    CHECK_THAT(phenol.logP,
               WithinAbs(-0.2893 + 0.5437 + 5 * 0.1581 + 5 * 0.1230 + (-0.2677), 1e-9));
}

TEST_CASE("Crippen MR is positive and monotone over a homologous series", "[chem][crippen]") {
    useInTreePack();
    double previous = 0.0;
    for (const char* smiles : {"C", "CC", "CCC", "CCCC"}) {
        const auto r = chem::crippen(parse(smiles));
        INFO(smiles);
        REQUIRE(r.ok);
        CHECK(r.molarRefractivity > 0.0);
        CHECK(r.molarRefractivity > previous);
        previous = r.molarRefractivity;
    }
}

TEST_CASE("chem::crippenLogP delegates to the faithful implementation", "[chem][crippen]") {
    useInTreePack();
    for (const char* smiles : {"Cn1cnc2c1c(=O)n(C)c(=O)n2C", "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
                               "NCCc1ccc(O)c(O)c1", "CCO"}) {
        const auto m = parse(smiles);
        CHECK_THAT(chem::crippenLogP(m), WithinAbs(chem::crippen(m).logP, 1e-12));
    }
}

TEST_CASE("Crippen reports a missing parameter pack instead of returning zero", "[chem][crippen]") {
    // Point the loader at a directory that cannot contain the pack. The cache
    // holds successful loads only, so this case cannot poison the others - but it
    // must restore the in-tree directory before it returns.
#if defined(_WIN32)
    _putenv_s("BIOCAD_DESCRIPTOR_DIR", "__nonexistent_biocad_descriptor_dir__");
#else
    setenv("BIOCAD_DESCRIPTOR_DIR", "/nonexistent-biocad-descriptor-dir", 1);
#endif
    // Only meaningful when the working directory does not also hold the pack; the
    // relative fallbacks are what make a dev run work, so this is a conditional
    // assertion rather than a fabricated guarantee.
    const bool relativeFallbackExists = std::filesystem::exists("assets/packs/descriptors/crippen.json");
    const auto r = chem::crippen(parse("CCO"));
    if (!relativeFallbackExists) {
        CHECK_FALSE(r.ok);
        CHECK(r.note.find("crippen.json") != std::string::npos);
        CHECK(r.logP == 0.0);
    }
    useInTreePack();
    CHECK(chem::crippen(parse("CCO")).ok);
}

TEST_CASE("Crippen carries its citation and the method's own error", "[chem][crippen]") {
    const std::string cite = chem::crippenCitation();
    CHECK(cite.find("Wildman-Crippen") != std::string::npos);
    CHECK(cite.find("1999;39:868-873") != std::string::npos);
    CHECK(cite.find("0.67") != std::string::npos);
}
