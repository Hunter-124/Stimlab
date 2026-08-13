// Residue-template connectivity: the bond graph a coordinate file does not contain.
//
// The load-bearing tests are the ones about what is NOT produced: no bond across a chain break,
// no bond for a residue with no template, and no distance-inferred bond unless the caller asked
// for one. A connectivity pass that quietly bonds everything within 1.9 A would pass a naive
// "did we get bonds" test and be wrong in exactly the ways that matter.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <string>

#include "bio/Connectivity.h"
#include "bio/PdbReader.h"

using namespace biocad::bio;

namespace {

std::filesystem::path fixture(const char* name) {
    return std::filesystem::path(BIOCAD_TEST_FIXTURES) / name;
}

std::filesystem::path templateFile() {
    return std::filesystem::path(BIOCAD_ASSETS_DIR) / "packs" / "structure" /
           "residue-templates.json";
}

// connect() reads the pack through core::assetDir(). Point that at the in-tree assets so the
// test does not depend on where the test binary was built.
void useInTreeAssets() {
    const std::string dir = BIOCAD_ASSETS_DIR;
#if defined(_WIN32)
    _putenv_s("BIOCAD_ASSETS", dir.c_str());
#else
    setenv("BIOCAD_ASSETS", dir.c_str(), 1);
#endif
}

int countKind(const ConnectivityResult& r, BondKind kind) {
    int n = 0;
    for (const StructureBond& b : r.bonds)
        if (b.kind == kind) ++n;
    return n;
}

double bondLength(const Model& m, const StructureBond& b) {
    const Atom& x = m.chains[static_cast<std::size_t>(b.a.chain)]
                        .residues[static_cast<std::size_t>(b.a.residue)]
                        .atoms[static_cast<std::size_t>(b.a.atom)];
    const Atom& y = m.chains[static_cast<std::size_t>(b.b.chain)]
                        .residues[static_cast<std::size_t>(b.b.residue)]
                        .atoms[static_cast<std::size_t>(b.b.atom)];
    return std::sqrt((x.x - y.x) * (x.x - y.x) + (x.y - y.y) * (x.y - y.y) +
                     (x.z - y.z) * (x.z - y.z));
}

}  // namespace

TEST_CASE("the residue template pack loads and covers the standard residues", "[bio][bonds]") {
    const ResidueTemplatePack pack = loadResidueTemplates(templateFile());
    REQUIRE(pack.errors.empty());
    REQUIRE(pack.ok);
    REQUIRE(pack.schemaVersion == 1);
    for (const char* name : {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
                             "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
                             "TYR", "VAL", "MSE", "A", "C", "G", "U", "DA", "DC", "DG", "DT"}) {
        REQUIRE(pack.find(name) != nullptr);
    }
    // The HIS tautomer spellings differ only in a proton, and this pack has no protons, so they
    // must resolve onto HIS itself rather than being unknown residues.
    REQUIRE(pack.find("HID") == pack.find("HIS"));
    REQUIRE(pack.find("HIE") == pack.find("HIS"));
    REQUIRE(pack.find("HIP") == pack.find("HIS"));
    REQUIRE(pack.isSolvent("HOH"));
    REQUIRE_FALSE(pack.isSolvent("ALA"));
    REQUIRE(pack.find("LIG") == nullptr);
}

TEST_CASE("template bonding reproduces a hand-counted tripeptide", "[bio][bonds]") {
    useInTreeAssets();
    const Structure st = readPdbFile(fixture("4m48_fragment.pdb"));
    const Model* m = st.model(1);
    REQUIRE(m != nullptr);
    const ConnectivityResult r = connect(*m);

    // GLY -1: N-CA, CA-C, C-O = 3. ALA 100: + CA-CB = 4. SER 100A: N-CA, CA-C, C-O, CA-CB and
    // CB-OG for BOTH altLoc copies of OG = 6. Two peptide bonds. CYS 5 in chain B: N-CA, CA-C,
    // C-O, CA-CB, CB-SG = 5. 3 + 4 + 6 + 5 = 18 template bonds and 2 links, 20 in total.
    REQUIRE(countKind(r, BondKind::Template) == 18);
    REQUIRE(countKind(r, BondKind::PeptideLink) == 2);
    REQUIRE(r.bonds.size() == 20);
    REQUIRE(countKind(r, BondKind::DistanceInferred) == 0);
    REQUIRE(r.diagnostics.gaps.empty());
    REQUIRE(r.diagnostics.unknownResidues.empty());
    // The water oxygen and the calcium ion are the only unbonded atoms, and neither is reported
    // as an unknown residue: there is no bond to make in either case.
    REQUIRE(r.diagnostics.unbondedAtoms == 2);
    REQUIRE(r.diagnostics.solventResidues == 1);
    REQUIRE(r.diagnostics.monatomicResidues == 1);

    for (const StructureBond& b : r.bonds) {
        if (b.kind != BondKind::PeptideLink) continue;
        REQUIRE_THAT(bondLength(*m, b), Catch::Matchers::WithinAbs(1.34, 0.05));
    }
}

TEST_CASE("an unknown residue is reported by name and left unbonded", "[bio][bonds]") {
    useInTreeAssets();
    const Structure st = readPdbFile(fixture("connectivity_probe.pdb"));
    const Model* m = st.model(1);
    REQUIRE(m != nullptr);
    const ConnectivityResult r = connect(*m);

    REQUIRE(r.diagnostics.unknownResidues.size() == 1);
    REQUIRE(r.diagnostics.unknownResidues[0].name == "LIG");
    REQUIRE(r.diagnostics.unknownResidues[0].count == 1);
    REQUIRE(r.diagnostics.unknownResidues[0].atoms == 4);
    // Its four atoms contributed no bond of any kind.
    for (const StructureBond& b : r.bonds) {
        const Residue& ra = m->chains[static_cast<std::size_t>(b.a.chain)]
                                .residues[static_cast<std::size_t>(b.a.residue)];
        const Residue& rb = m->chains[static_cast<std::size_t>(b.b.chain)]
                                .residues[static_cast<std::size_t>(b.b.residue)];
        REQUIRE(ra.name != "LIG");
        REQUIRE(rb.name != "LIG");
    }
    bool named = false;
    for (const std::string& w : r.diagnostics.warnings)
        if (w.find("LIG") != std::string::npos) named = true;
    REQUIRE(named);
}

TEST_CASE("a chain break is a reported gap and never a bond", "[bio][bonds]") {
    useInTreeAssets();
    const Structure st = readPdbFile(fixture("connectivity_probe.pdb"));
    const Model* m = st.model(1);
    REQUIRE(m != nullptr);
    const ConnectivityResult r = connect(*m);

    REQUIRE(r.diagnostics.gaps.size() == 1);
    const ChainGap& g = r.diagnostics.gaps[0];
    REQUIRE(g.chainId == "A");
    REQUIRE(g.fromSeqId == 3);
    REQUIRE(g.toSeqId == 11);
    REQUIRE(g.distance > 10.0);
    // Four residues in the chain, one break: exactly two peptide bonds, not three.
    REQUIRE(countKind(r, BondKind::PeptideLink) == 2);
    for (const StructureBond& b : r.bonds) REQUIRE(bondLength(*m, b) < 2.6);
}

TEST_CASE("a disulfide is found from SG-SG geometry with no SSBOND record", "[bio][bonds]") {
    useInTreeAssets();
    const Structure st = readPdbFile(fixture("connectivity_probe.pdb"));
    const Model* m = st.model(1);
    REQUIRE(m != nullptr);
    const ConnectivityResult r = connect(*m);

    REQUIRE(countKind(r, BondKind::Disulfide) == 1);
    REQUIRE(r.diagnostics.disulfides == 1);
    for (const StructureBond& b : r.bonds) {
        if (b.kind != BondKind::Disulfide) continue;
        REQUIRE_THAT(bondLength(*m, b), Catch::Matchers::WithinAbs(2.03, 1e-3));
    }

    ConnectivityOptions off;
    off.findDisulfides = false;
    REQUIRE(countKind(connect(*m, off), BondKind::Disulfide) == 0);
}

TEST_CASE("the distance fallback is off by default and labels what it invents", "[bio][bonds]") {
    useInTreeAssets();
    const Structure st = readPdbFile(fixture("connectivity_probe.pdb"));
    const Model* m = st.model(1);
    REQUIRE(m != nullptr);

    const ConnectivityResult plain = connect(*m);
    REQUIRE(countKind(plain, BondKind::DistanceInferred) == 0);
    REQUIRE(plain.diagnostics.inferredBonds == 0);

    ConnectivityOptions opts;
    opts.inferByDistance = true;
    const ConnectivityResult inferred = connect(*m, opts);
    // The ligand's three short contacts, and nothing else: a templated residue is never
    // re-bonded by distance, so the template count is untouched.
    REQUIRE(countKind(inferred, BondKind::DistanceInferred) == 3);
    REQUIRE(inferred.diagnostics.inferredBonds == 3);
    REQUIRE(countKind(inferred, BondKind::Template) == countKind(plain, BondKind::Template));
    REQUIRE(countKind(inferred, BondKind::PeptideLink) == countKind(plain, BondKind::PeptideLink));
}
