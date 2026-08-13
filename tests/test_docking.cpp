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

using namespace biocad;
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

TEST_CASE("Receptor presets load from the data packs with real boxes", "[docking][presets]") {
    const auto& presets = docking::targetPresets();
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

TEST_CASE("Every listed target resolves for on-demand provisioning", "[docking][presets]") {
    // The Docking/Workflows combos list presetNames(); selecting one and provisioning
    // it on demand resolves that display name back to a preset via findPreset - exactly
    // the path AppShell::provisionTarget / receptorReady take. Guard that round-trip for
    // every preset so a pack edit cannot silently break the on-demand UI.
    const auto names = docking::presetNames();
    REQUIRE(names.size() >= 25);
    for (const auto& n : names) {
        INFO(n);
        const auto* byName = docking::findPreset(n);          // combo value -> preset
        REQUIRE(byName != nullptr);
        REQUIRE_FALSE(byName->id.empty());
        REQUIRE_FALSE(byName->pdb.empty());                   // a PDB to fetch from RCSB
        REQUIRE(byName->box.sx > 1.0);                        // a real, non-degenerate box
        REQUIRE(docking::findPreset(byName->id) != nullptr);  // id also resolves (CLI/selftest)
    }
    // The 4 headlines resolve, and a representative NON-headline target (D2) is a fully
    // provisionable preset - the new on-demand path covers the other 25 presets too.
    for (const char* id : {"DAT", "NET", "SERT", "TAAR1"})
        REQUIRE(docking::findPreset(id) != nullptr);
    const auto* d2 = docking::findPreset("D2");
    REQUIRE(d2 != nullptr);
    REQUIRE(d2->id == "D2");
    REQUIRE_FALSE(d2->pdb.empty());
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

TEST_CASE("Flexible PDBQT writer builds a valid torsion tree", "[docking][pdbqt]") {
    const auto m = graph("CC(N)Cc1ccccc1");  // amphetamine
    const auto conf = chem::embed3D(m);
    const auto lig = docking::writeFlexiblePdbqt(m, conf);
    const auto rigid = docking::writeRigidPdbqt(m, conf);

    REQUIRE(lig.atomCount > 0);
    REQUIRE(lig.atomCount == lig.heavyCount + lig.polarH);
    REQUIRE(lig.atomCount == rigid.atomCount);                       // same atoms, different tree
    REQUIRE(static_cast<int>(lig.serialToConf.size()) == lig.atomCount);
    REQUIRE(lig.text.find("ROOT") != std::string::npos);
    REQUIRE(lig.text.find("ENDROOT") != std::string::npos);

    // Amphetamine has exactly two acyclic rotatable bonds (Calpha-CH2, CH2-aryl);
    // the aromatic ring and the terminal methyl/amine bonds are not torsions.
    REQUIRE(lig.torsions == 2);
    REQUIRE(lig.text.find("TORSDOF 2") != std::string::npos);

    auto count = [](const std::string& s, const std::string& sub) {
        size_t n = 0, p = 0;
        while ((p = s.find(sub, p)) != std::string::npos) { ++n; p += sub.size(); }
        return n;
    };
    // Every BRANCH is matched by an ENDBRANCH, one pair per active torsion.
    REQUIRE(count(lig.text, "\nBRANCH ") == 2);
    REQUIRE(count(lig.text, "\nENDBRANCH ") == 2);
    REQUIRE(lig.text.find("nan") == std::string::npos);
    REQUIRE(lig.text.find("inf") == std::string::npos);

    // Round-trip: the engine preserves atom order in its output, so a pose parsed
    // with serialToConf must scatter coordinates back onto the ORIGINAL topology
    // (tree order != conformer order). Wrap the written ATOM block as one MODEL.
    const std::string out =
        "MODEL 1\nREMARK VINA RESULT:    -7.00      0.000      0.000\n" + lig.text + "ENDMDL\n";
    const auto poses = docking::parseVinaPdbqt(out, conf, &lig.serialToConf);
    REQUIRE(poses.size() == 1);
    const auto& pose = poses[0].ligand;
    REQUIRE(pose.heavyCount == conf.heavyCount);          // full topology from the reference
    REQUIRE(pose.z.size() == conf.z.size());
    for (int k = 0; k < lig.atomCount; ++k) {
        const int idx = lig.serialToConf[k];
        REQUIRE(idx >= 0);
        REQUIRE(idx < static_cast<int>(pose.pos.size()));
        REQUIRE_THAT(pose.pos[idx].x, WithinAbs(conf.pos[idx].x, 1e-2));
        REQUIRE_THAT(pose.pos[idx].y, WithinAbs(conf.pos[idx].y, 1e-2));
        REQUIRE_THAT(pose.pos[idx].z, WithinAbs(conf.pos[idx].z, 1e-2));
    }
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
    REQUIRE_FALSE(e.fromEngine());
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
    REQUIRE_FALSE(detail.fromEngine());                        // no receptor -> labeled estimate
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

TEST_CASE("Modified residues (MSE) are kept as protein, not stripped", "[docking][receptor]") {
    // MSE (selenomethionine) is recorded as HETATM but is covalently part of the
    // chain - a large fraction of crystal structures use it. Stripping it would tear
    // a gap in the backbone and disconnect a pocket residue. It must survive cleanup,
    // be retyped as MET, and its Se-delta must dock as the thioether sulfur (S).
    std::string pdb;
    pdb += pdbAtom("ATOM",   1, "N",  "ALA", 'A', 1,  0.000, 0.000, 0.000, "N");
    pdb += pdbAtom("ATOM",   2, "CA", "ALA", 'A', 1,  1.458, 0.000, 0.000, "C");
    pdb += pdbAtom("HETATM", 3, "N",  "MSE", 'A', 2,  2.000, 1.400, 0.000, "N");
    pdb += pdbAtom("HETATM", 4, "CA", "MSE", 'A', 2,  3.400, 1.700, 0.000, "C");
    pdb += pdbAtom("HETATM", 5, "CB", "MSE", 'A', 2,  3.900, 3.100, 0.000, "C");
    pdb += pdbAtom("HETATM", 6, "CG", "MSE", 'A', 2,  5.400, 3.300, 0.000, "C");
    pdb += pdbAtom("HETATM", 7, "SE", "MSE", 'A', 2,  6.100, 5.000, 0.000, "Se");  // selenium
    pdb += pdbAtom("HETATM", 8, "CE", "MSE", 'A', 2,  7.900, 4.800, 0.000, "C");
    pdb += pdbAtom("HETATM", 9, "O",  "HOH", 'A', 99, 40.000, 40.000, 40.000, "O");  // still stripped
    pdb += "END\n";

    const auto rec = docking::pdbToRigidReceptor(pdb);

    // 2 ALA + 6 MSE heavy atoms kept; the water is the only thing dropped.
    REQUIRE(rec.atomCount == 8);
    REQUIRE(rec.keptModified == 6);
    REQUIRE(rec.keptHetero == 0);
    REQUIRE(rec.droppedHetero == 1);

    // The residue is canonicalised to MET (no MSE label leaks through) and the
    // selenium is emitted as a sulfur atom type, never the AutoDock default carbon.
    REQUIRE(rec.text.find("MSE") == std::string::npos);
    REQUIRE(rec.text.find("MET") != std::string::npos);
    REQUIRE(rec.text.find(" SD ") != std::string::npos);  // SE renamed to S-delta
    REQUIRE(rec.text.find("Se") == std::string::npos);    // no stray selenium element
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
    docking::ReceptorPrepResult miss = docking::locatePreparedReceptor("__biocad_absent_id__");
    REQUIRE_FALSE(miss.ready);
    REQUIRE(miss.note.find("no prepared receptor") != std::string::npos);

    // locate-only ensureReceptor never reaches the network and stays in fallback.
    ReceptorTarget t;
    t.id = "__biocad_absent_id__";
    t.pdb = "0XYZ";
    REQUIRE_FALSE(docking::ensureReceptor(t, /*allowDownload=*/false).ready);

    // When a prepared PDBQT IS on disk, the lookup finds it and returns its path.
    std::error_code ec;
    std::filesystem::create_directories(docking::receptorsDir(), ec);
    const auto path = docking::receptorsDir() / "__biocad_unit__.pdbqt";
    {
        std::ofstream o(path, std::ios::binary);
        o << docking::pdbToRigidReceptor(syntheticPdb()).text;
    }
    docking::ReceptorPrepResult hit = docking::locatePreparedReceptor("__biocad_unit__");
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

TEST_CASE("Vina-GPU backend locates, labels, and degrades honestly", "[docking][vinagpu]") {
    // Identity is stable - realEngines() ordering and the UI labels depend on it.
    docking::VinaGpuBackend gpu;
    REQUIRE(gpu.id() == "vina-gpu");
    REQUIRE(gpu.displayName() == "GPU (OpenCL, Vina-GPU)");

    // The engine gets its own subfolder of runtime/engines (it carries companion files).
    REQUIRE(docking::vinaGpuDir() == docking::enginesDir() / "vina-gpu");

    // Locate-only provisioning never touches the network and never throws; whatever it
    // reports as fetched must agree with the backend's own availability (exe + kernel).
    const auto probe = docking::ensureVinaGpu(/*allowDownload=*/false);
    REQUIRE_FALSE(probe.note.empty());
    REQUIRE(probe.fetched == gpu.available());

    // available() is exactly "exe located AND a compiled kernel beside it" - true or
    // false depending on whether this host has provisioned Vina-GPU, but always consistent.
    const auto bin = docking::locateEngine(docking::Engine::VinaGpu);
    const bool kernel = bin && std::filesystem::exists(bin->parent_path() / "Kernel2_Opt.bin");
    REQUIRE(gpu.available() == (bin.has_value() && kernel));

    // With no prepared receptor the dock degrades to real=false (never throws), so the
    // caller can fall back to the next engine / the labeled estimate.
    const auto m = graph("CC(N)Cc1ccccc1");
    const auto conf = chem::embed3D(m);
    ReceptorTarget tgt;
    tgt.id = "DAT";
    tgt.name = "DAT";  // receptorPath intentionally empty -> early, hermetic real=false
    const auto d = gpu.dock(m, conf, tgt);
    REQUIRE_FALSE(d.fromEngine());
    REQUIRE_FALSE(d.log.empty());
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
