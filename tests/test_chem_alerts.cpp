// Tests for the structural-alert screen (chem/Alerts + RealAlerts).
//
// Two things are under test here and the second matters more than the first:
// (1) the alerts match the motifs they claim to, on the canonical literature
//     examples, and do not match compounds that lack them; and
// (2) NOTHING in this feature can produce a toxicity verdict. An alert is a
//     liability flag, so Verdict::Danger must be unreachable - asserted over
//     every compound in the shipped library, not just the interesting ones.
#include <catch2/catch_test_macros.hpp>

#include <string>
#include <vector>

#include "chem/Alerts.h"
#include "chem/Smiles.h"
#include "modules/AlertsModule.h"
#include "modules/RealBackend.h"

using namespace biocad;

namespace {

std::vector<chem::AlertHit> screen(const std::string& smiles) {
    auto m = chem::parseSmiles(smiles);
    REQUIRE(m.has_value());
    return chem::screenAlerts(*m);
}

bool flagged(const std::vector<chem::AlertHit>& hits, const std::string& key) {
    for (const auto& h : hits) {
        if (h.key == key) return true;
    }
    return false;
}

}  // namespace

TEST_CASE("The bioactivation alert pack loads cleanly", "[chem][alerts]") {
    // A rule that failed to parse is dropped and reported. If the pack has errors
    // the screen is quietly weaker than it looks, which is the failure mode this
    // assertion exists to catch.
    for (const auto& e : chem::alertPackErrors()) WARN(e);
    REQUIRE(chem::alertPackErrors().empty());
    REQUIRE(chem::alertRuleCount() >= 12);
}

TEST_CASE("Acetaminophen raises the para-aminophenol quinone-imine flag", "[chem][alerts]") {
    // The canonical case: CYP oxidation of the 4-aminophenol gives NAPQI. This is
    // why the compound is in the library at all.
    const auto hits = screen("CC(=O)Nc1ccc(O)cc1");
    REQUIRE(flagged(hits, "para-aminophenol"));
    for (const auto& h : hits) {
        REQUIRE(!h.atoms.empty());
        REQUIRE(!h.citation.empty());
    }
}

TEST_CASE("Nitroaromatic matches; benzene and caffeine match nothing", "[chem][alerts]") {
    REQUIRE(flagged(screen("c1ccc(cc1)[N+](=O)[O-]"), "nitroaromatic"));
    // Both nitro tautomer spellings must hit: a rule pack that depends on how the
    // author drew the charge separation is a rule pack that silently misses.
    REQUIRE(flagged(screen("O=N(=O)c1ccccc1"), "nitroaromatic"));
    REQUIRE(screen("c1ccccc1").empty());
    REQUIRE(screen("Cn1cnc2c1c(=O)n(C)c(=O)n2C").empty());
}

TEST_CASE("A Michael acceptor matches and its saturated analogue does not",
          "[chem][alerts]") {
    REQUIRE(flagged(screen("C=CC(=O)C"), "michael-acceptor"));
    REQUIRE(screen("CCC(=O)C").empty());
}

TEST_CASE("Acyl-glucuronide precursor needs an acid on a lipophilic scaffold",
          "[chem][alerts]") {
    REQUIRE(flagged(screen("CC(C)Cc1ccc(cc1)C(C)C(=O)O"), "acyl-glucuronide-precursor"));
    // Ethanol has no acid at all; acetic acid has one but no lipophilic scaffold,
    // and the migration liability is a property of the scaffold, not the -COOH.
    REQUIRE(screen("CCO").empty());
    REQUIRE(!flagged(screen("CC(=O)O"), "acyl-glucuronide-precursor"));
}

TEST_CASE("Each authored motif matches its own reference structure", "[chem][alerts]") {
    REQUIRE(flagged(screen("CCOc1ccc(NC(C)=O)cc1"), "para-alkoxyaniline"));   // phenacetin
    REQUIRE(flagged(screen("NCCc1ccc(O)c(O)c1"), "catechol"));                // dopamine
    REQUIRE(flagged(screen("Oc1ccc(O)cc1"), "hydroquinone"));
    REQUIRE(flagged(screen("Nc1ccccc1"), "primary-aniline"));                 // aniline
    REQUIRE(flagged(screen("c1ccsc1"), "thiophene"));
    REQUIRE(flagged(screen("c1ccoc1"), "furan"));
    REQUIRE(flagged(screen("C#CCO"), "terminal-alkyne"));
    REQUIRE(flagged(screen("NNC(=O)c1ccncc1"), "hydrazine"));                 // isoniazid
    REQUIRE(flagged(screen("C1OC1"), "epoxide"));
    REQUIRE(flagged(screen("NC(=S)N"), "thiourea"));
    REQUIRE(flagged(screen("ClCC(=O)NC"), "alpha-haloamide"));
}

TEST_CASE("Every alert carries a non-empty citation", "[chem][alerts]") {
    // An uncitable alert is exactly the fabricated rule this project removes, so
    // the pack loader drops one - and this asserts none survives into a hit.
    const std::vector<std::string> probes = {
        "CC(=O)Nc1ccc(O)cc1", "c1ccc(cc1)[N+](=O)[O-]", "Nc1ccccc1", "c1ccsc1", "c1ccoc1",
        "C#CCO", "NNC(=O)c1ccncc1", "C1OC1", "C=CC(=O)C", "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "NC(=S)N", "ClCC(=O)NC", "NCCc1ccc(O)c(O)c1", "CCOc1ccc(NC(C)=O)cc1"};
    std::size_t total = 0;
    for (const auto& smiles : probes) {
        for (const auto& h : screen(smiles)) {
            REQUIRE(!h.citation.empty());
            REQUIRE(!h.label.empty());
            REQUIRE(!h.mechanism.empty());
            ++total;
        }
    }
    REQUIRE(total >= probes.size());
}

TEST_CASE("No alert yields Verdict::Danger for any library compound", "[chem][alerts]") {
    RealBackend backend;
    Services svc = backend.services();
    REQUIRE(svc.alerts != nullptr);
    REQUIRE(svc.library != nullptr);

    const auto compounds = svc.library->all();
    REQUIRE(compounds.size() > 20);

    std::size_t flags = 0;
    for (const auto& m : compounds) {
        const AlertReport r = svc.alerts->screen(m);
        REQUIRE(!r.summary.empty());
        for (const auto& f : r.flags) {
            // The whole point: a substructure match cannot support a toxicity
            // verdict, so only these two levels are reachable.
            REQUIRE((f.severity == Verdict::Info || f.severity == Verdict::Warn));
            REQUIRE(f.severity != Verdict::Danger);
            REQUIRE(f.severity != Verdict::Good);
            REQUIRE(!f.citation.empty());
            REQUIRE(f.atomCount > 0);
            ++flags;
        }
    }
    // The library contains acetaminophen and several NSAIDs, so a run that raised
    // nothing at all would mean the pack silently failed to load.
    REQUIRE(flags > 0);
}

TEST_CASE("An unmatched screen states the non-claim explicitly", "[chem][alerts]") {
    RealAlerts alerts;
    Molecule benzene;
    benzene.id = "benzene";
    benzene.smiles = "c1ccccc1";
    const AlertReport r = alerts.screen(benzene);
    REQUIRE(r.flags.empty());
    REQUIRE(r.summary.find("not a safety claim") != std::string::npos);

    // An unparseable structure was not screened, and must not read as "clean".
    Molecule broken;
    broken.id = "broken";
    broken.smiles = "C(((";
    const AlertReport b = alerts.screen(broken);
    REQUIRE(b.flags.empty());
    REQUIRE(b.summary.find("could not be parsed") != std::string::npos);
}
