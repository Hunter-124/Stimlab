#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <string>
#include <vector>

#include "chem/Embed3D.h"
#include "chem/Smiles.h"
#include "modules/RealBackend.h"

using namespace stimlab::chem;
using Catch::Matchers::WithinAbs;

namespace {

double dist(const Vec3& a, const Vec3& b) {
    const double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Ideal length of a bond from covalent radii + order (mirrors Embed3D's model).
double idealBond(const Molecule& m, const Bond& b) {
    const double base = covalentRadius(m.atoms[b.a].z) + covalentRadius(m.atoms[b.b].z);
    if (b.aromatic) return base * 0.92;
    if (b.order == 3.0) return base * 0.78;
    if (b.order == 2.0) return base * 0.87;
    return base;
}

Molecule mustParse(const std::string& smi) {
    auto m = parseSmiles(smi);
    REQUIRE(m.has_value());
    return *m;
}

// Smallest distance between any two distinct atoms in the conformer.
double minPairwise(const Conformer& c) {
    double mn = 1e9;
    for (int i = 0; i < c.size(); ++i)
        for (int j = i + 1; j < c.size(); ++j) mn = std::min(mn, dist(c.pos[i], c.pos[j]));
    return mn;
}

}  // namespace

TEST_CASE("Embed3D produces correct bond lengths", "[chem][embed]") {
    // Heavy bonds occupy the first m.bonds.size() entries of the conformer bond
    // list, in graph order, so each maps back to its parsed bond + order.
    for (const char* smi : {"CC(N)Cc1ccccc1",                    // amphetamine
                            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"}) {   // caffeine
        const auto m = mustParse(smi);
        const auto c = embed3D(m);
        REQUIRE(c.heavyCount == m.heavyAtomCount());
        REQUIRE(static_cast<int>(c.bonds.size()) >= static_cast<int>(m.bonds.size()));

        std::vector<double> ccSingle;
        for (size_t k = 0; k < m.bonds.size(); ++k) {
            const Bond& b = m.bonds[k];
            const double d = dist(c.pos[b.a], c.pos[b.b]);
            const double ideal = idealBond(m, b);
            INFO(smi << " bond " << b.a << "-" << b.b << " ideal=" << ideal << " got=" << d);
            REQUIRE_THAT(d, WithinAbs(ideal, 0.20));
            if (!b.aromatic && b.order == 1.0 && m.atoms[b.a].z == 6 && m.atoms[b.b].z == 6)
                ccSingle.push_back(d);
        }
        // The brief's acceptance: single C-C bonds land near 1.5 A.
        for (double d : ccSingle) REQUIRE_THAT(d, WithinAbs(1.52, 0.15));
    }
}

TEST_CASE("Embed3D leaves no atom overlaps", "[chem][embed]") {
    for (const char* smi : {"CC(N)Cc1ccccc1",                     // amphetamine
                            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",       // caffeine
                            "COC(=O)C1C(OC(=O)c2ccccc2)CC3CCC1N3C"}) {  // cocaine
        const auto c = embed3D(mustParse(smi));
        INFO(smi << " min pairwise = " << minPairwise(c));
        REQUIRE(minPairwise(c) > 0.80);  // no coincident / clashing atoms
    }
}

TEST_CASE("Embed3D embeds every library compound without NaN", "[chem][embed][library]") {
    stimlab::RealBackend backend;
    auto svc = backend.services();
    REQUIRE(svc.library->count() >= 28);

    for (const auto& dm : svc.library->all()) {
        const auto pm = parseSmiles(dm.smiles);
        REQUIRE(pm.has_value());
        const auto c = embed3D(*pm);
        INFO(dm.name << "  (" << dm.smiles << ")");
        REQUIRE(c.heavyCount == pm->heavyAtomCount());
        REQUIRE_FALSE(c.empty());
        for (const auto& p : c.pos) {
            REQUIRE(std::isfinite(p.x));
            REQUIRE(std::isfinite(p.y));
            REQUIRE(std::isfinite(p.z));
        }
        REQUIRE(minPairwise(c) > 0.70);
    }
}
