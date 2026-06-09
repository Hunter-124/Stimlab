#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <thread>

#include "chem/Embed3D.h"
#include "chem/Smiles.h"
#include "modules/RealBackend.h"
#include "modules/docking/Backends.h"
#include "modules/docking/EngineLocator.h"
#include "modules/docking/PdbqtWriter.h"
#include "modules/docking/Presets.h"
#include "modules/docking/Provisioning.h"
#include "modules/docking/ReceptorPrep.h"
#include "modules/docking/VinaParser.h"

using namespace stimlab;
using Catch::Matchers::WithinAbs;

namespace {
chem::Molecule graph(const char* smi) {
    auto m = chem::parseSmiles(smi);
    REQUIRE(m.has_value());
    return *m;
}

// Emit one strictly column-aligned PDB ATOM/HETATM record so the fixed-column
// receptor parser sees coordinates and the element field exactly where RCSB writes
// them (x at col 31, element at cols 77-78).
std::string pdbAtom(const char* rec, int serial, const char* name, const char* res, char chain,
                    int resSeq, double x, double y, double z, const char* elem) {
    char b[100];
    std::snprintf(b, sizeof(b),
                  "%-6.6s%5d %-4.4s %-3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s\n",
                  rec, serial, name, res, chain, resSeq, x, y, z, 1.0, 20.0, elem);
    return b;
}

// A tiny synthetic receptor fragment: a PHE (aromatic ring), a SER hydroxyl, a MET
// thioether, one explicit H, plus a water and a sulfate to be stripped.
std::string syntheticPdb() {
    std::string p;
    p += pdbAtom("ATOM", 1, "N", "PHE", 'A', 1, 11.104, 13.207, 10.000, "N");
    p += pdbAtom("ATOM", 2, "CA", "PHE", 'A', 1, 12.560, 13.207, 10.000, "C");
    p += pdbAtom("ATOM", 3, "C", "PHE", 'A', 1, 13.000, 14.500, 10.500, "C");
    p += pdbAtom("ATOM", 4, "O", "PHE", 'A', 1, 12.300, 15.400, 10.900, "O");
    p += pdbAtom("ATOM", 5, "CB", "PHE", 'A', 1, 13.100, 12.000, 10.700, "C");
    p += pdbAtom("ATOM", 6, "CG", "PHE", 'A', 1, 14.600, 12.000, 10.700, "C");   // aromatic
    p += pdbAtom("ATOM", 7, "CD1", "PHE", 'A', 1, 15.300, 13.100, 11.100, "C");  // aromatic
    p += pdbAtom("ATOM", 8, "CE1", "PHE", 'A', 1, 16.700, 13.100, 11.100, "C");  // aromatic
    p += pdbAtom("ATOM", 9, "CZ", "PHE", 'A', 1, 17.400, 12.000, 10.700, "C");   // aromatic
    p += pdbAtom("ATOM", 10, "OG", "SER", 'A', 2, 18.000, 10.000, 9.000, "O");
    p += pdbAtom("ATOM", 11, "SD", "MET", 'A', 3, 19.000, 8.000, 8.000, "S");
    p += pdbAtom("ATOM", 12, "H", "PHE", 'A', 1, 11.000, 12.500, 9.500, "H");  // dropped (nonpolar)
    p += pdbAtom("HETATM", 13, "O", "HOH", 'A', 101, 20.000, 20.000, 20.000, "O");  // water, stripped
    p += pdbAtom("HETATM", 14, "S", "SO4", 'A', 102, 21.000, 21.000, 21.000, "S");  // additive, stripped
    p += "END\n";
    return p;
}
}  // namespace

TEST_CASE("CNS receptor presets are populated with boxes", "[docking][presets]") {
    const auto& presets = docking::cnsPresets();
    REQUIRE(presets.size() >= 25);
    for (const auto& t : presets) {
        INFO(t.name);
        REQUIRE_FALSE(t.id.empty());
        REQUIRE_FALSE(t.pdb.empty());
        REQUIRE(t.box.sx > 1.0);   // non-degenerate search box
        REQUIRE(t.box.sy > 1.0);
        REQUIRE(t.box.sz > 1.0);
    }
    // The historical transporters resolve by id and by display name.
    REQUIRE(docking::findPreset("DAT") != nullptr);
    REQUIRE(docking::findPreset("SERT") != nullptr);
    REQUIRE(docking::presetNames().size() == presets.size());
}

TEST_CASE("Rigid PDBQT writer emits a clean ROOT ligand", "[docking][pdbqt]") {
    const auto m = graph("CC(N)Cc1ccccc1");  // amphetamine
    const auto conf = chem::embed3D(m);
    const auto lig = docking::writeRigidPdbqt(m, conf);

    REQUIRE(lig.atomCount > 0);
    REQUIRE(lig.atomCount == lig.heavyCount + lig.polarH);
    REQUIRE(lig.text.find("ROOT") != std::string::npos);
    REQUIRE(lig.text.find("ENDROOT") != std::string::npos);
    REQUIRE(lig.text.find("TORSDOF 0") != std::string::npos);  // rigid, no Meeko tree
    // No NaN/inf leaked into the coordinate columns.
    REQUIRE(lig.text.find("nan") == std::string::npos);
    REQUIRE(lig.text.find("inf") == std::string::npos);
}

TEST_CASE("Vina output parser ranks poses with conformers", "[docking][parser]") {
    const auto m = graph("CC(N)Cc1ccccc1");
    const auto conf = chem::embed3D(m);
    const auto lig = docking::writeRigidPdbqt(m, conf);

    // Synthesize a 2-MODEL Vina output by wrapping the same ATOM block under two
    // REMARK VINA RESULT scores (weaker one written first to prove sorting).
    const std::string model =
        "MODEL %d\nREMARK VINA RESULT:    %s      0.000      0.000\n" + lig.text + "ENDMDL\n";
    auto fill = [&](int n, const char* score) {
        std::string s = model;
        s.replace(s.find("%d"), 2, std::to_string(n));
        s.replace(s.find("%s"), 2, score);
        return s;
    };
    const std::string out = fill(1, "-6.20") + fill(2, "-7.80");

    const auto poses = docking::parseVinaPdbqt(out, conf);
    REQUIRE(poses.size() == 2);
    REQUIRE(poses[0].affinityKcalPerMol < poses[1].affinityKcalPerMol);   // strongest first
    REQUIRE_THAT(poses[0].affinityKcalPerMol, WithinAbs(-7.80, 0.01));
    REQUIRE(poses[0].rank == 1);
    REQUIRE_FALSE(poses[0].ligand.empty());
    REQUIRE(poses[0].ligand.heavyCount == conf.heavyCount);
}

TEST_CASE("Docking degrades gracefully to a labeled estimate", "[docking][fallback]") {
    // No engine/receptor is provisioned in this environment, so the real backends
    // are unavailable and the descriptor estimate must carry the result.
    const auto m = graph("CC(N)Cc1ccccc1");
    const auto conf = chem::embed3D(m);
    docking::EstimateBackend est;
    ReceptorTarget tgt;
    tgt.id = "DAT";
    const auto e = est.dock(m, conf, tgt);
    REQUIRE_FALSE(e.real);
    REQUIRE(e.poses.size() == 6);
    for (const auto& p : e.poses) {
        REQUIRE(std::isfinite(p.affinityKcalPerMol));
        REQUIRE_FALSE(p.ligand.empty());   // viewer always has geometry to show
    }
    REQUIRE(e.poses.front().affinityKcalPerMol <= e.poses.back().affinityKcalPerMol);
}

TEST_CASE("RealDocking wires backends and keeps both views consistent", "[docking][real]") {
    RealBackend backend;
    Services s = backend.services();
    REQUIRE(s.docking->presets().size() >= 25);
    REQUIRE_FALSE(s.docking->targets().empty());

    const auto mol = s.library->byId("amphetamine");
    REQUIRE(mol.has_value());

    // Dock into a target with NO prepared receptor so the path is the deterministic,
    // hermetic descriptor estimate (the engine returns early before any subprocess),
    // independent of whatever a dev box may have provisioned under runtime/receptors.
    // The live real-engine path is covered by the --selftest-dock acceptance run.
    const std::string target = "__unprepared_test_target__";
    const auto detail = s.docking->dockDetailed(*mol, target);
    REQUIRE_FALSE(detail.poses.empty());
    REQUIRE(std::isfinite(detail.bestAffinity()));
    REQUIRE_FALSE(detail.real);                        // no receptor -> labeled estimate
    REQUIRE(detail.engine == "descriptor-estimate");

    const auto legacy = s.docking->dock(*mol, target);
    REQUIRE_THAT(legacy.bestAffinity, WithinAbs(detail.bestAffinity(), 1e-9));
    REQUIRE(legacy.poses.size() == detail.poses.size());
}

// ---------------------------------------------------------------- WP-F receptor prep
TEST_CASE("Built-in receptor prep strips solvent and types a rigid receptor", "[docking][receptor]") {
    const auto rec = docking::pdbToRigidReceptor(syntheticPdb());

    // 11 protein heavy atoms kept; the explicit H is dropped (united-atom receptor),
    // the water + sulfate HETATM (2 atoms) are stripped.
    REQUIRE(rec.atomCount == 11);
    REQUIRE(rec.keptHetero == 0);
    REQUIRE(rec.droppedHetero == 2);

    // Solvent / additive residues never leak into the receptor.
    REQUIRE(rec.text.find("HOH") == std::string::npos);
    REQUIRE(rec.text.find("SO4") == std::string::npos);

    // Aromatic PHE ring carbons get the AutoDock 'A' type; oxygens get 'OA'.
    REQUIRE(rec.text.find(" A \n") != std::string::npos);
    REQUIRE(rec.text.find("OA") != std::string::npos);

    // A receptor is rigid by definition: no ligand torsion tree, no NaN/inf.
    REQUIRE(rec.text.find("ROOT") == std::string::npos);
    REQUIRE(rec.text.find("TORSDOF") == std::string::npos);
    REQUIRE(rec.text.find("nan") == std::string::npos);
    REQUIRE(rec.text.find("inf") == std::string::npos);
    REQUIRE(rec.text.rfind("REMARK", 0) == 0);  // starts with our provenance header

    // A box center is derived in the structure's own frame. With no large co-crystal
    // ligand here (SO4 has < 8 atoms), it falls back to the receptor centroid, which
    // must be finite and inside the coordinate span of the kept atoms.
    REQUIRE(rec.hasBox);
    REQUIRE(rec.boxSource == "receptor centroid");
    REQUIRE(std::isfinite(rec.cx));
    REQUIRE(std::isfinite(rec.cy));
    REQUIRE(std::isfinite(rec.cz));
    REQUIRE(rec.cx > 10.0);
    REQUIRE(rec.cx < 20.0);
}

TEST_CASE("Box center is taken from the co-crystal ligand when present", "[docking][receptor]") {
    // A 10-atom HETATM ligand far from the protein should win the box center over the
    // protein centroid (its centroid is ~ (50,50,50)).
    std::string pdb;
    pdb += pdbAtom("ATOM", 1, "CA", "ALA", 'A', 1, 0.0, 0.0, 0.0, "C");
    pdb += pdbAtom("ATOM", 2, "CB", "ALA", 'A', 1, 1.5, 0.0, 0.0, "C");
    for (int i = 0; i < 10; ++i)
        pdb += pdbAtom("HETATM", 100 + i, "C", "LIG", 'A', 500, 50.0, 50.0, 50.0 + i * 0.01, "C");
    pdb += "END\n";

    double cx = 0, cy = 0, cz = 0;
    std::string src;
    REQUIRE(docking::receptorBoxFromPdb(pdb, cx, cy, cz, src));
    REQUIRE(src.find("co-crystal ligand") != std::string::npos);
    REQUIRE(src.find("LIG") != std::string::npos);
    REQUIRE_THAT(cx, WithinAbs(50.0, 0.5));
    REQUIRE_THAT(cy, WithinAbs(50.0, 0.5));
}

TEST_CASE("Receptor prep can retain a named cofactor", "[docking][receptor]") {
    // Keeping SO4 by name retains its atom instead of stripping it.
    const auto rec = docking::pdbToRigidReceptor(syntheticPdb(), {"SO4"});
    REQUIRE(rec.atomCount == 12);
    REQUIRE(rec.keptHetero == 1);
    REQUIRE(rec.droppedHetero == 1);  // only the water now
}

TEST_CASE("Prepared-receptor cache lookup is honest", "[docking][receptor]") {
    // A target with no cached PDBQT is reported not-ready (no network on this path).
    docking::ReceptorPrepResult miss = docking::locatePreparedReceptor("__stimlab_absent_id__");
    REQUIRE_FALSE(miss.ready);
    REQUIRE(miss.note.find("no prepared receptor") != std::string::npos);

    // locate-only ensureReceptor never reaches the network and stays in fallback.
    ReceptorTarget t;
    t.id = "__stimlab_absent_id__";
    t.pdb = "0XYZ";
    REQUIRE_FALSE(docking::ensureReceptor(t, /*allowDownload=*/false).ready);

    // When a prepared PDBQT IS on disk, the lookup finds it and returns its path.
    std::error_code ec;
    std::filesystem::create_directories(docking::receptorsDir(), ec);
    const auto path = docking::receptorsDir() / "__stimlab_unit__.pdbqt";
    {
        std::ofstream o(path, std::ios::binary);
        o << docking::pdbToRigidReceptor(syntheticPdb()).text;
    }
    docking::ReceptorPrepResult hit = docking::locatePreparedReceptor("__stimlab_unit__");
    REQUIRE(hit.ready);
    REQUIRE(hit.path == path.string());
    std::filesystem::remove(path, ec);  // clean up the unit artifact
}

TEST_CASE("Vina integrity pin is well-formed", "[docking][provision]") {
    // No env/file/compile pin is configured by default, so this is empty; if a pin
    // is present it must be a 64-hex SHA-256.
    const std::string h = docking::expectedVinaSha256();
    REQUIRE((h.empty() || h.size() == 64));
}

TEST_CASE("Provisioner runs a locate-only probe without blocking or networking",
          "[docking][provision]") {
    docking::Provisioner prov;
    REQUIRE_FALSE(prov.status().empty());

    // Locate-only (no download): completes near-instantly and touches no network.
    prov.start(/*allowDownload=*/false, docking::headlinePresets());
    REQUIRE(prov.everRun());
    REQUIRE(prov.receptorsTotal() == 4);  // DAT/NET/SERT/TAAR1

    for (int i = 0; i < 400 && prov.running(); ++i)
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    REQUIRE_FALSE(prov.running());
    REQUIRE(prov.receptorsReady() <= prov.receptorsTotal());
    REQUIRE_FALSE(prov.status().empty());
}
