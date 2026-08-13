#include <catch2/catch_test_macros.hpp>

#include <cstdint>

#include "chem/Embed3D.h"
#include "chem/Smiles.h"
#include "render/MolViewport.h"

using namespace biocad;

namespace {
chem::Conformer embed(const char* smi) {
    auto m = chem::parseSmiles(smi);
    REQUIRE(m.has_value());
    return chem::embed3D(*m);
}
}  // namespace

// The CPU scene builder is the unit-testable half of the viewer (the GPU half
// needs a device). These assertions pin its instance accounting + CPK coloring.
TEST_CASE("buildMolScene accounts atoms and bonds per style", "[render][scene]") {
    const auto conf = embed("CC(N)Cc1ccccc1");  // amphetamine
    REQUIRE(conf.size() > 0);

    // Ball-and-stick, hydrogens shown: one sphere per atom, one cylinder per bond.
    const auto bs = render::buildMolScene(conf, /*spacefill=*/false, /*showH=*/true);
    REQUIRE(bs.atoms.size() == static_cast<std::size_t>(conf.size()));
    REQUIRE(bs.bonds.size() == conf.bonds.size());
    REQUIRE(bs.radius > 0.0f);

    // Hydrogens hidden: only heavy atoms, and no bond touching a hydrogen.
    const auto noH = render::buildMolScene(conf, false, false);
    REQUIRE(noH.atoms.size() == static_cast<std::size_t>(conf.heavyCount));
    REQUIRE(noH.bonds.size() < bs.bonds.size());

    // Spacefill: van-der-Waals spheres, no cylinders.
    const auto sf = render::buildMolScene(conf, true, true);
    REQUIRE(sf.bonds.empty());
    REQUIRE(sf.atoms.size() == static_cast<std::size_t>(conf.size()));
}

TEST_CASE("buildMolScene tolerates an empty conformer", "[render][scene]") {
    const auto sc = render::buildMolScene(chem::Conformer{}, false, true);
    REQUIRE(sc.atoms.empty());
    REQUIRE(sc.bonds.empty());
    REQUIRE(sc.radius >= 1.0f);  // never zero (camera divides by it)
}

TEST_CASE("cpkColor distinguishes elements and is opaque", "[render][cpk]") {
    REQUIRE(render::cpkColor(7) != render::cpkColor(8));   // N != O
    REQUIRE(render::cpkColor(6) != render::cpkColor(1));   // C != H
    REQUIRE(render::cpkColor(8) != render::cpkColor(6));   // O != C
    // Packed as 0xAABBGGRR with a fully-opaque alpha byte.
    REQUIRE((render::cpkColor(6) & 0xFF000000u) == 0xFF000000u);
    REQUIRE((render::cpkColor(99) & 0xFF000000u) == 0xFF000000u);  // default still opaque
}
