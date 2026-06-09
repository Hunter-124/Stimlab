#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <string>

#include "chem/Embed3D.h"
#include "chem/Smiles.h"
#include "modules/RealBackend.h"
#include "modules/docking/Backends.h"
#include "modules/docking/PdbqtWriter.h"
#include "modules/docking/Presets.h"
#include "modules/docking/VinaParser.h"

using namespace stimlab;
using Catch::Matchers::WithinAbs;

namespace {
chem::Molecule graph(const char* smi) {
    auto m = chem::parseSmiles(smi);
    REQUIRE(m.has_value());
    return *m;
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
    const auto detail = s.docking->dockDetailed(*mol, "DAT");
    REQUIRE_FALSE(detail.poses.empty());
    REQUIRE(std::isfinite(detail.bestAffinity()));
    REQUIRE_FALSE(detail.real);  // graceful fallback in this environment

    const auto legacy = s.docking->dock(*mol, "DAT");
    REQUIRE_THAT(legacy.bestAffinity, WithinAbs(detail.bestAffinity(), 1e-9));
    REQUIRE(legacy.poses.size() == detail.poses.size());
}
