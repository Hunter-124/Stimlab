#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <optional>
#include <string>

#include "modules/RealBackend.h"

using namespace stimlab;
using Catch::Matchers::WithinAbs;

namespace {
Molecule mol(Services& s, const std::string& id) {
    auto m = s.library->byId(id);
    REQUIRE(m.has_value());
    return *m;
}
bool hasEndpoint(const AdmetReport& r, const std::string& needle) {
    for (const auto& e : r.endpoints)
        if (e.name.find(needle) != std::string::npos) return true;
    return false;
}
}  // namespace

TEST_CASE("Real library computes properties from structure", "[modules][real]") {
    RealBackend backend;
    Services s = backend.services();
    REQUIRE(s.valid());
    REQUIRE(s.library->count() >= 28);

    const auto amp = mol(s, "amphetamine");
    REQUIRE(amp.formula == "C9H13N");                 // computed, not stored
    REQUIRE_THAT(amp.molWeight, WithinAbs(135.21, 0.1));
    REQUIRE_THAT(amp.tpsa, WithinAbs(26.02, 0.5));
    REQUIRE(amp.hbd == 1);
    REQUIRE(amp.hba == 1);
}

// Golden ADMET facts (plan.md / handoff.md section 9), now derived from REAL
// structure perception rather than hardcoded.
TEST_CASE("Real ADMET goldens hold (structure-derived)", "[modules][admet][golden]") {
    RealBackend backend;
    Services s = backend.services();

    const auto amp = s.admet->screen(mol(s, "amphetamine"));
    REQUIRE(amp.overall == Verdict::Warn);
    REQUIRE(hasEndpoint(amp, "MAO"));
    REQUIRE(hasEndpoint(amp, "CYP2D6"));

    REQUIRE(s.admet->screen(mol(s, "methamphetamine")).overall == Verdict::Warn);
    REQUIRE(s.admet->screen(mol(s, "mdma")).overall == Verdict::Warn);

    const auto dopa = s.admet->screen(mol(s, "dopamine"));
    REQUIRE(dopa.overall == Verdict::Warn);
    REQUIRE(hasEndpoint(dopa, "COMT"));

    REQUIRE(hasEndpoint(s.admet->screen(mol(s, "acetaminophen")), "NAPQI"));
    REQUIRE(s.admet->screen(mol(s, "caffeine")).overall == Verdict::Info);
}

TEST_CASE("Real similarity uses fingerprints and ranks analogs", "[modules][similarity]") {
    RealBackend backend;
    Services s = backend.services();

    const auto sim = s.similarity->search(mol(s, "amphetamine"));
    REQUIRE_FALSE(sim.hits.empty());
    REQUIRE(sim.nearestScore > 0.3);
    for (std::size_t i = 1; i < sim.hits.size(); ++i)
        REQUIRE(sim.hits[i - 1].tanimoto >= sim.hits[i].tanimoto);
}

TEST_CASE("Real absorption is CNS-aware", "[modules][absorption]") {
    RealBackend backend;
    Services s = backend.services();
    REQUIRE(s.absorption->predict(mol(s, "methamphetamine")).cnsPenetrant);
    REQUIRE_FALSE(s.absorption->predict(mol(s, "norepinephrine")).cnsPenetrant);  // polar catechol
}

TEST_CASE("Real stability flags ester hydrolysis", "[modules][stability]") {
    RealBackend backend;
    Services s = backend.services();
    const auto cocaine = s.stability->analyze(mol(s, "cocaine"));
    const auto caffeine = s.stability->analyze(mol(s, "caffeine"));
    REQUIRE(cocaine.overallScore < caffeine.overallScore);
}
