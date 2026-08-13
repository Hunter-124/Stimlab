// Protein-protein interfaces and the geometric alanine scan, against three committed
// RCSB fixtures: 1BRS (barnase/barstar), 1VFB (anti-lysozyme D1.3) and 1DQJ
// (anti-lysozyme HyHEL-63). Each fixture is the deposited entry with only the chains
// under test and their non-water records kept.
//
// The buried-area expectations are RANGES, not point values, and they are around the
// interface areas commonly reported for these complexes (barnase/barstar buries about
// 1560 A^2 in total, i.e. ~780 A^2 per side; the two antibody/lysozyme interfaces bury
// roughly 700-900 A^2 per side). A point assertion would be wrong on principle: an
// interface area moves by ~10% with the probe radius, the radii set and the hydrogen
// policy, which is exactly why every area here carries its parameter string.
#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <string>

#include "bio/Interface.h"
#include "bio/PdbReader.h"

using namespace biocad;
using namespace biocad::bio;

namespace {

Structure loadFixture(const std::string& id) {
    const std::filesystem::path p = std::filesystem::path(BIOCAD_TEST_FIXTURES) / (id + ".pdb");
    Structure s = readPdbFile(p);
    REQUIRE_FALSE(s.models.empty());
    return s;
}

}  // namespace

TEST_CASE("chain lists parse both spellings", "[biologics]") {
    CHECK(parseChainList("AB") == std::vector<std::string>{"A", "B"});
    CHECK(parseChainList("H,L") == std::vector<std::string>{"H", "L"});
    CHECK(parseChainList("AA,BB") == std::vector<std::string>{"AA", "BB"});
}

TEST_CASE("barnase/barstar buries the published area with its hotspot contacts", "[biologics]") {
    const Structure s = loadFixture("1BRS");
    const InterfaceReport r = interfaceOf(s, "A", "D");

    REQUIRE(r.buriedSurfaceArea.provenance == Provenance::Measured);
    CHECK(r.buriedSurfaceArea.unit == "A^2");
    CHECK_FALSE(r.sasaParameters.empty());
    // BSA is the TOTAL area buried on both sides, so the per-side figure is half.
    const double perSide = r.buriedSurfaceArea.value / 2.0;
    CHECK(perSide > 700.0);
    CHECK(perSide < 900.0);
    CHECK(r.buriedSurfaceArea.value ==
          r.sasaA.value + r.sasaB.value - r.sasaComplex.value);

    int hbonds = 0, saltBridges = 0;
    for (const auto& c : r.contacts) {
        hbonds += c.hydrogenBond ? 1 : 0;
        saltBridges += c.saltBridge ? 1 : 0;
    }
    CHECK(r.contacts.size() > 30);
    CHECK(hbonds >= 10);
    // Barstar Asp39 salt-bridges to barnase Arg83/Arg87 in the active site.
    CHECK(saltBridges >= 1);
    CHECK_FALSE(r.coreResidues.empty());
    CHECK_FALSE(r.supportResidues.empty());
    // No antibody chain was named, so nothing is guessed about paratopes.
    CHECK(r.cdrContacts.empty());
    CHECK(r.paratope.empty());
    CHECK(r.epitope.empty());
}

TEST_CASE("antibody/lysozyme interfaces report CDR contacts, paratope and epitope",
          "[biologics]") {
    struct Case {
        const char* id;
        double lo, hi;
    };
    for (const Case c : {Case{"1VFB", 550.0, 900.0}, Case{"1DQJ", 600.0, 1000.0}}) {
        const Structure s = loadFixture(c.id);
        InterfaceOptions opts;
        opts.antibodyChains = {"A", "B"};
        const InterfaceReport r = interfaceOf(s, "AB", "C", opts);
        const double perSide = r.buriedSurfaceArea.value / 2.0;
        CHECK(perSide > c.lo);
        CHECK(perSide < c.hi);
        CHECK(r.contacts.size() > 20);
        CHECK_FALSE(r.epitope.empty());
        CHECK_FALSE(r.cdrContacts.empty());
        // Framework contacts are REPORTED, not filtered: the paratope is at least the
        // CDR set and any residue outside it is named in a warning.
        CHECK(r.cdrContacts.size() <= r.paratope.size());
        if (r.cdrContacts.size() < r.paratope.size()) CHECK_FALSE(r.warnings.empty());
    }
}

TEST_CASE("a chain that is not in the model is a warning, not a silent empty interface",
          "[biologics]") {
    const Structure s = loadFixture("1BRS");
    const InterfaceReport r = interfaceOf(s, "A", "Z");
    REQUIRE_FALSE(r.warnings.empty());
    CHECK(r.contacts.empty());
}

TEST_CASE("every alanine-scan impact is unit-free and Heuristic", "[biologics]") {
    const Structure s = loadFixture("1BRS");
    const AlanineScanReport scan = alanineScan(s, "A", "D");
    REQUIRE_FALSE(scan.positions.empty());

    for (const auto& p : scan.positions) {
        // The whole point: an empty unit and the Heuristic tier. makeQuantity() throws
        // on a Heuristic that carries a unit, so this cannot regress quietly.
        CHECK(p.impact.unit.empty());
        CHECK(p.impact.provenance == Provenance::Heuristic);
        CHECK(p.impact.value >= 0.0);
        CHECK(p.impact.source.find("NOT a ddG") != std::string::npos);
        // Gly and Ala have nothing beyond C-beta, so they are never scanned.
        CHECK(p.residue.find(":GLY") == std::string::npos);
        CHECK(p.residue.find(":ALA") == std::string::npos);
    }
    // Descending by impact, so the top row is the strongest geometric contributor.
    for (std::size_t i = 1; i < scan.positions.size(); ++i)
        CHECK(scan.positions[i - 1].impact.value >= scan.positions[i].impact.value);

    // No hard-coded benchmark correlation, ever.
    CHECK(scan.benchmarkSpearman.provenance == Provenance::NotComputed);
    CHECK(scan.benchmarkSpearman.value == 0.0);
    CHECK(scan.benchmarkSpearman.source == "a measured benchmark subset");
    CHECK(scan.disclaimer.find("NOT a binding free energy") != std::string::npos);

    // The known barnase/barstar hotspots (barnase Arg59/Lys27/His102, barstar
    // Asp39/Asp35/Trp44) should dominate the top of a geometric ranking. This asserts
    // the RANKING, not an energy.
    std::string top;
    for (std::size_t i = 0; i < scan.positions.size() && i < 8; ++i)
        top += scan.positions[i].residue + " ";
    CHECK(top.find("ARG59") != std::string::npos);
    CHECK(top.find("ASP39") != std::string::npos);
}
