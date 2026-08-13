// tests/test_backend.cpp - the REAL backend, under test.
//
// There is deliberately no fake backend in this tree. A second implementation of
// the same science is a second set of answers, and a suite that validates the
// double while the product ships the original is a suite that proves nothing.
// Every assertion below runs the shipping code: RealBackend, the in-house chem
// engine, and the data packs on disk.
//
// This stays hermetic without a double: the library is loaded from
// assets/packs/*.json, every property is computed from SMILES, and the docking
// module's hot path is cache-only (it never downloads and never spawns an engine
// that has not already been provisioned).
#include <catch2/catch_test_macros.hpp>

#include <optional>
#include <string>

#include "modules/RealBackend.h"

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

TEST_CASE("The library is loaded from the data packs", "[backend][library]") {
    RealBackend backend;
    Services s = backend.services();
    REQUIRE(s.valid());

    // Four built-in packs ship 67 compounds; assert a floor rather than the exact
    // count so adding a pack is not a test failure.
    REQUIRE(s.library->count() >= 60);
    REQUIRE(s.library->byId("amphetamine").has_value());
    REQUIRE(s.library->byId("caffeine").has_value());
    REQUIRE(s.library->byId("ibuprofen").has_value());        // analgesics-otc pack
    REQUIRE(s.library->byId("metformin").has_value());        // metabolic pack
    REQUIRE_FALSE(s.library->byId("__absent_compound__").has_value());
}

TEST_CASE("Every property is computed from the structure, not authored",
          "[backend][library][chem]") {
    RealBackend backend;
    Services s = backend.services();

    // Caffeine, C8H10N4O2, MW 194.19 g/mol. The pack supplies only identity and a
    // SMILES; formula and mass come out of chem::Descriptors, so this asserts the
    // engine, not a transcribed literal.
    const auto caf = mol(s, "caffeine");
    REQUIRE(caf.formula == "C8H10N4O2");
    REQUIRE(caf.molWeight > 194.0);
    REQUIRE(caf.molWeight < 194.5);
    REQUIRE(caf.hba > 0);

    // Amphetamine, C9H13N, MW 135.21.
    const auto amp = mol(s, "amphetamine");
    REQUIRE(amp.formula == "C9H13N");
    REQUIRE(amp.molWeight > 135.0);
    REQUIRE(amp.molWeight < 135.5);
    REQUIRE(amp.hbd >= 1);

    // Ibuprofen, C13H18O2, MW 206.28 - a pack that authors no properties at all.
    const auto ibu = mol(s, "ibuprofen");
    REQUIRE(ibu.formula == "C13H18O2");
    REQUIRE(ibu.molWeight > 206.0);
    REQUIRE(ibu.molWeight < 206.5);
}

TEST_CASE("ADMET endpoints follow the perceived structure", "[backend][admet]") {
    RealBackend backend;
    Services s = backend.services();

    const auto amp = s.admet->screen(mol(s, "amphetamine"));
    REQUIRE(amp.overall == Verdict::Warn);
    REQUIRE(hasEndpointContaining(amp, "MAO"));
    REQUIRE(hasEndpointContaining(amp, "CYP2D6"));

    REQUIRE(s.admet->screen(mol(s, "methamphetamine")).overall == Verdict::Warn);
    REQUIRE(s.admet->screen(mol(s, "mdma")).overall == Verdict::Warn);

    // Dopamine is a catechol, so COMT must appear.
    const auto dopa = s.admet->screen(mol(s, "dopamine"));
    REQUIRE(dopa.overall == Verdict::Warn);
    REQUIRE(hasEndpointContaining(dopa, "COMT"));

    // Acetaminophen's bioactivation route is the reason it is in the library.
    REQUIRE(hasEndpointContaining(s.admet->screen(mol(s, "acetaminophen")), "NAPQI"));

    // Caffeine has none of those liabilities.
    REQUIRE(s.admet->screen(mol(s, "caffeine")).overall == Verdict::Info);
}

TEST_CASE("Stability ranks an ester below a xanthine", "[backend][stability]") {
    RealBackend backend;
    Services s = backend.services();

    const auto cocaine = s.stability->analyze(mol(s, "cocaine"));    // two esters
    const auto caffeine = s.stability->analyze(mol(s, "caffeine"));  // no labile group
    REQUIRE(cocaine.overallScore < caffeine.overallScore);
    REQUIRE_FALSE(cocaine.degradants.empty());
    // A shelf life is an extrapolation of measured degradation rates, and analyze()
    // is given only a structure - so it must come back NotComputed naming the missing
    // measurement rather than mapping a flag count onto "~12 months".
    REQUIRE(cocaine.shelfLife.provenance == Provenance::NotComputed);
    REQUIRE(cocaine.shelfLife.unit.empty());
    REQUIRE(cocaine.shelfLife.source.find("three or more temperatures") != std::string::npos);
}

TEST_CASE("Absorption is bounded and tracks polarity", "[backend][absorption]") {
    RealBackend backend;
    Services s = backend.services();

    const auto meth = s.absorption->predict(mol(s, "methamphetamine"));
    REQUIRE(meth.bioavailabilityPct >= 0.0);
    REQUIRE(meth.bioavailabilityPct <= 100.0);
    REQUIRE(meth.hiaPct >= 0.0);
    REQUIRE(meth.hiaPct <= 100.0);
    REQUIRE(meth.cnsPenetrant);  // low TPSA, lipophilic
    REQUIRE_FALSE(meth.metrics.empty());

    // A polar catechol must not read as blood-brain-barrier permeant.
    REQUIRE_FALSE(s.absorption->predict(mol(s, "dopamine")).cnsPenetrant);
}

TEST_CASE("Similarity is bounded and sorted descending", "[backend][similarity]") {
    RealBackend backend;
    Services s = backend.services();

    const auto sim = s.similarity->search(mol(s, "methcathinone"));
    REQUIRE_FALSE(sim.hits.empty());
    REQUIRE(sim.nearestScore >= 0.0);
    REQUIRE(sim.nearestScore <= 1.0);
    for (size_t i = 1; i < sim.hits.size(); ++i)
        REQUIRE(sim.hits[i - 1].tanimoto >= sim.hits[i].tanimoto);
}

TEST_CASE("The docking hot path is cache-only and labels its provenance",
          "[backend][docking][provenance]") {
    RealBackend backend;
    Services s = backend.services();

    // No receptor is prepared for this id and no engine can claim it, so this must
    // come back as the labelled descriptor estimate - never as a silent failure and
    // never as a Model-tier score. If an engine IS provisioned on this machine the
    // target still has no receptor, so the fallback is what is exercised either way.
    const auto d = s.docking->dockDetailed(mol(s, "amphetamine"), "__unprepared_target__");
    REQUIRE_FALSE(d.fromEngine());
    REQUIRE(d.provenance == Provenance::Heuristic);
    REQUIRE_FALSE(d.log.empty());

    // Targets come from the packs, and every one of them resolves.
    const auto names = s.docking->targets();
    REQUIRE_FALSE(names.empty());
    REQUIRE(s.docking->presets().size() == names.size());
}
