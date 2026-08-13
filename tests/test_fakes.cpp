#include <catch2/catch_test_macros.hpp>

#include <optional>
#include <string>

#include "fakes/FakeBackend.h"

using namespace biocad;

namespace {
Molecule mol(Services& s, const std::string& id) {
    auto m = s.library->byId(id);
    REQUIRE(m.has_value());
    return *m;
}
bool hasEndpointContaining(const AdmetReport& r, const std::string& needle) {
    for (const auto& e : r.endpoints)
        if (e.name.find(needle) != std::string::npos) return true;
    return false;
}
}  // namespace

TEST_CASE("Default library is populated", "[fakes][library]") {
    FakeBackend backend;
    Services s = backend.services();
    REQUIRE(s.valid());
    REQUIRE(s.library->count() >= 14);
    REQUIRE(s.library->byId("amphetamine").has_value());
    REQUIRE(s.library->byId("caffeine").has_value());
}

// Golden ADMET facts backed by validated reference values.
TEST_CASE("ADMET goldens hold", "[fakes][admet][golden]") {
    FakeBackend backend;
    Services s = backend.services();

    const auto amp = s.admet->screen(mol(s, "amphetamine"));
    REQUIRE(amp.overall == Verdict::Warn);
    REQUIRE(hasEndpointContaining(amp, "MAO"));
    REQUIRE(hasEndpointContaining(amp, "CYP2D6"));

    const auto meth = s.admet->screen(mol(s, "methamphetamine"));
    REQUIRE(meth.overall == Verdict::Warn);

    const auto mdma = s.admet->screen(mol(s, "mdma"));
    REQUIRE(mdma.overall == Verdict::Warn);

    const auto dopa = s.admet->screen(mol(s, "dopamine"));
    REQUIRE(dopa.overall == Verdict::Warn);
    REQUIRE(hasEndpointContaining(dopa, "COMT"));

    const auto apap = s.admet->screen(mol(s, "acetaminophen"));
    REQUIRE(hasEndpointContaining(apap, "NAPQI"));

    const auto caf = s.admet->screen(mol(s, "caffeine"));
    REQUIRE(caf.overall == Verdict::Info);
}

TEST_CASE("Stability flags ester hydrolysis as the limiting factor", "[fakes][stability]") {
    FakeBackend backend;
    Services s = backend.services();

    const auto cocaine = s.stability->analyze(mol(s, "cocaine"));
    const auto caffeine = s.stability->analyze(mol(s, "caffeine"));
    // Ester-bearing cocaine should be less stable overall than caffeine.
    REQUIRE(cocaine.overallScore < caffeine.overallScore);
    REQUIRE_FALSE(cocaine.degradants.empty());
    REQUIRE_FALSE(cocaine.shelfLifeEstimate.empty());
}

TEST_CASE("Absorption model produces bounded, CNS-aware metrics", "[fakes][absorption]") {
    FakeBackend backend;
    Services s = backend.services();

    const auto meth = s.absorption->predict(mol(s, "methamphetamine"));
    REQUIRE(meth.bioavailabilityPct >= 0.0);
    REQUIRE(meth.bioavailabilityPct <= 100.0);
    REQUIRE(meth.hiaPct >= 0.0);
    REQUIRE(meth.hiaPct <= 100.0);
    REQUIRE(meth.cnsPenetrant);          // low TPSA, lipophilic -> crosses BBB
    REQUIRE_FALSE(meth.metrics.empty());

    // Dopamine (polar catechol) should not be predicted CNS-penetrant.
    const auto dopa = s.absorption->predict(mol(s, "dopamine"));
    REQUIRE_FALSE(dopa.cnsPenetrant);
}

TEST_CASE("Similarity ranks a near analog at the top", "[fakes][similarity]") {
    FakeBackend backend;
    Services s = backend.services();

    const auto sim = s.similarity->search(mol(s, "methcathinone"));
    REQUIRE_FALSE(sim.hits.empty());
    REQUIRE(sim.nearestScore <= 1.0);
    REQUIRE(sim.nearestScore >= 0.0);
    // hits are sorted descending by tanimoto
    for (size_t i = 1; i < sim.hits.size(); ++i)
        REQUIRE(sim.hits[i - 1].tanimoto >= sim.hits[i].tanimoto);
}
