// tests/test_mechanism.cpp - Phase 7: retrieved mechanism, panel coverage, pathway
// context, interaction flags and pharmacogenomic notes, under test.
//
// Every assertion runs the shipping code: the real packs on disk, the real
// RealMechanism, and the real docking module through RealBackend. There is no
// double, and the cases that would need a network assert the OFFLINE contract -
// networkAvailable == false with a stated reason and no invented content - because
// that is the contract the default build actually ships.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "modules/MechanismModule.h"
#include "modules/RealBackend.h"

using namespace biocad;

namespace {

std::filesystem::path packDir() {
    return std::filesystem::path(BIOCAD_ASSETS_DIR) / "packs" / "mechanism";
}

MechanismPacks loaded() {
    MechanismPacks p = loadMechanismPacks(packDir());
    REQUIRE(p.errors.empty());
    return p;
}

std::string lower(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

bool has(const std::string& haystack, const std::string& needle) {
    return haystack.find(needle) != std::string::npos;
}

Molecule caffeine() {
    Molecule m;
    m.id = "caffeine";
    m.name = "Caffeine";
    m.smiles = "Cn1cnc2c1c(=O)n(C)c(=O)n2C";
    return m;
}

}  // namespace

TEST_CASE("Mechanism packs load from the shipped asset tree", "[mechanism]") {
    const MechanismPacks p = loaded();
    REQUIRE(p.panels.size() == 3);
    REQUIRE(p.actionTypes.values.size() >= 30);
    REQUIRE(p.interactions.members.size() >= 40);
    REQUIRE_FALSE(p.interactions.boundaryNote.empty());
    REQUIRE_FALSE(p.interactions.classRules.empty());
    REQUIRE_FALSE(p.pharmacogenomics.notes.empty());

    // At least 40 hand-curated supplement / food / lifestyle entries, which is the
    // half of the interaction pack the FDA tables do not cover.
    const auto& ms = p.interactions.members;
    const long curated = std::count_if(ms.begin(), ms.end(), [](const StackMember& m) {
        return m.kind == "supplement" || m.kind == "food" || m.kind == "lifestyle";
    });
    REQUIRE(curated >= 40);
    // Every member carries a citation, enforced by the loader.
    for (const StackMember& m : ms) REQUIRE_FALSE(m.citation.empty());
}

TEST_CASE("Panel packs declare their own coverage honestly", "[mechanism]") {
    const MechanismPacks p = loaded();

    const MechanismPanelPack* ss44 = p.panel("safetyscreen44");
    REQUIRE(ss44 != nullptr);
    REQUIRE(ss44->targets.size() == 44);
    REQUIRE(ss44->declaredSize == 44);

    // The extended panel declares 87 and enumerates fewer. The declared size is what
    // coverage is measured against, so the rows BioCAD cannot even name stay counted.
    const MechanismPanelPack* ss87 = p.panel("safetyscreen87");
    REQUIRE(ss87 != nullptr);
    REQUIRE(ss87->declaredSize == 87);
    REQUIRE(static_cast<int>(ss87->targets.size()) < ss87->declaredSize);

    const MechanismPanelPack* cipa = p.panel("cipa-currents");
    REQUIRE(cipa != nullptr);
    REQUIRE(cipa->targets.size() == 5);
}

TEST_CASE("An unsupported pack schemaVersion is an error, never a silent skip", "[mechanism]") {
    MechanismPacks p;
    parseMechanismDocument(R"({"schemaVersion":99,"kind":"panel","panelId":"x"})", "<test>", p);
    REQUIRE(p.panels.empty());
    REQUIRE(p.errors.size() == 1);
    REQUIRE(has(p.errors[0], "schemaVersion 99"));
}

TEST_CASE("A pack shipping deprecated phenotype nomenclature is rejected", "[mechanism]") {
    MechanismPacks p;
    parseMechanismDocument(
        R"({"schemaVersion":1,"kind":"pharmacogenomics","id":"t","boundaryStatement":"b",
            "phenotypes":[{"code":"EM","term":"extensive metabolizer"},
                          {"code":"NM","term":"normal metabolizer"}]})",
        "<test>", p);
    REQUIRE(p.pharmacogenomics.phenotypes.size() == 1);
    REQUIRE(p.pharmacogenomics.phenotypes[0].code == "NM");
    REQUIRE(p.errors.size() == 1);
    REQUIRE(has(p.errors[0], "extensive metabolizer"));
}

TEST_CASE("The CPIC vocabulary is UM/RM/NM/IM/PM with no deprecated term", "[mechanism]") {
    const MechanismPacks p = loaded();
    REQUIRE(p.pharmacogenomics.phenotypes.size() == 5);
    std::string codes;
    for (const PgxPhenotype& ph : p.pharmacogenomics.phenotypes) {
        codes += ph.code;
        REQUIRE_FALSE(has(lower(ph.term), "extensive metabolizer"));
    }
    REQUIRE(codes == "UMRMNMIMPM");
    // The CYP2D6 activity-score bands are carried as data so no panel restates them.
    bool pm = false, um = false;
    for (const PgxActivityBand& b : p.pharmacogenomics.bands) {
        if (b.gene == "CYP2D6" && b.phenotype == "PM") { pm = true; REQUIRE(b.band == "AS = 0"); }
        if (b.gene == "CYP2D6" && b.phenotype == "UM") { um = true; REQUIRE(b.band == "AS > 2.25"); }
    }
    REQUIRE(pm);
    REQUIRE(um);
}

TEST_CASE("No forbidden data source appears in a mechanism pack", "[mechanism]") {
    std::string text;
    for (const auto& e : std::filesystem::directory_iterator(packDir())) {
        std::ifstream in(e.path());
        std::ostringstream ss;
        ss << in.rdbuf();
        text += ss.str();
    }
    const std::string low = lower(text);
    // Licence-blocked sources: not queried, not bundled, not named as a source.
    REQUIRE_FALSE(has(low, "kegg.jp/rest"));
    REQUIRE_FALSE(has(low, "rest.kegg.jp"));
    REQUIRE_FALSE(has(low, "string-db"));
    REQUIRE_FALSE(has(low, "drugbank"));
    REQUIRE_FALSE(has(low, "pdbbind"));
    REQUIRE_FALSE(has(low, "biolip"));
    REQUIRE_FALSE(has(low, "guidetopharmacology"));

    // The distinction that matters: a kegg.jp PATHWAY PAGE deep link is allowed,
    // because linking is not redistribution. It is a page URL, not an API URL.
    REQUIRE(keggPathwayUrl("hsa04726") == "https://www.kegg.jp/pathway/hsa04726");
    REQUIRE_FALSE(has(keggPathwayUrl("hsa04726"), "/rest"));
    REQUIRE(has(keggDeepLinkNote(), "linking is not redistribution"));
}

TEST_CASE("Services::valid() holds with the mechanism module wired", "[mechanism]") {
    RealBackend backend;
    Services s = backend.services();
    REQUIRE(s.mechanism != nullptr);
    REQUIRE(s.valid());
}

TEST_CASE("A panel screen reports unscreened = panelSize - screened", "[mechanism]") {
    RealBackend backend;
    Services s = backend.services();
    const PanelScreenReport r = s.mechanism->screenPanel(caffeine(), "safetyscreen44");

    REQUIRE(r.panelSize == 44);
    REQUIRE(r.results.size() == 44);
    REQUIRE(r.unscreened == r.panelSize - r.screened);
    // The coverage statement names BOTH numbers, so a reader cannot see the hits
    // without seeing the size of the unknown.
    REQUIRE(has(r.coverageStatement, std::to_string(r.unscreened)));
    REQUIRE(has(r.coverageStatement, "44"));
    REQUIRE(has(r.coverageStatement, "NOT screened"));
    REQUIRE(has(r.coverageStatement, "no composite safety score"));
    // Each row carries its own preparation and box, or a stated reason for neither.
    for (const PanelTargetResult& row : r.results) {
        if (row.screened) {
            REQUIRE_FALSE(row.receptorPreparation.empty());
            REQUIRE_FALSE(row.boxDefinition.empty());
            REQUIRE(row.affinity.provenance == Provenance::Model);
        } else {
            REQUIRE_FALSE(row.skipReason.empty());
            REQUIRE(row.affinity.provenance == Provenance::NotComputed);
        }
    }
}

TEST_CASE("A missing panel roster is reported, not fabricated", "[mechanism]") {
    RealMechanism m(nullptr, {}, loaded());
    const PanelScreenReport r = m.screenPanel(caffeine(), "no-such-panel");
    REQUIRE(r.results.empty());
    REQUIRE(r.screened == 0);
    REQUIRE(r.hergSafetyMargin.provenance == Provenance::NotComputed);
    REQUIRE_FALSE(r.warnings.empty());
}

TEST_CASE("The hERG margin exists only from user-supplied measurements", "[mechanism]") {
    RealMechanism m(nullptr, {}, loaded());

    SECTION("no measured IC50 -> NotComputed naming the prerequisite") {
        const PanelScreenReport r = m.screenPanel(caffeine(), "cipa-currents");
        REQUIRE(r.hergSafetyMargin.provenance == Provenance::NotComputed);
        REQUIRE(has(r.hergSafetyMargin.source, "measured hERG IC50"));
        REQUIRE(has(r.hergSafetyMargin.source, "does not predict"));
    }
    SECTION("a measured IC50 with no free Cmax is still NotComputed") {
        HergInput in;
        in.measuredIc50Molar = 1.0e-5;
        m.setHergInput(in);
        const PanelScreenReport r = m.screenPanel(caffeine(), "cipa-currents");
        REQUIRE(r.hergSafetyMargin.provenance == Provenance::NotComputed);
        REQUIRE(has(r.hergSafetyMargin.source, "free Cmax"));
    }
    SECTION("both measurements -> a dimensionless Measured ratio, flagged below 30") {
        HergInput in;
        in.measuredIc50Molar = 1.0e-5;   // 10 uM
        in.freeCmaxMolar = 1.0e-6;       // 1 uM free
        in.citation = "a patch-clamp measurement entered by the user";
        m.setHergInput(in);
        const PanelScreenReport r = m.screenPanel(caffeine(), "cipa-currents");
        REQUIRE(r.hergSafetyMargin.provenance == Provenance::Measured);
        REQUIRE(r.hergSafetyMargin.unit.empty());   // a ratio has no unit
        REQUIRE(r.hergSafetyMargin.value > 9.999);
        REQUIRE(r.hergSafetyMargin.value < 10.001);
        REQUIRE(r.hergMarginFlagBelow == 30.0);
        const bool flagged = std::any_of(r.warnings.begin(), r.warnings.end(),
                                         [](const std::string& w) {
                                             return w.find("below the flag threshold") !=
                                                    std::string::npos;
                                         });
        REQUIRE(flagged);
    }
}

TEST_CASE("The stack checker flags a mechanism with a citation, never a severity", "[mechanism]") {
    RealMechanism m(nullptr, {}, loaded());
    const StackReport r = m.checkStack({"caffeine", "fluvoxamine", "unobtainium-xr"});

    // Caffeine is a CYP1A2 substrate; fluvoxamine is a CYP1A2 inhibitor.
    const auto flag = std::find_if(r.flags.begin(), r.flags.end(), [](const InteractionFlag& f) {
        return f.leftId == "fluvoxamine" && f.rightId == "caffeine" &&
               f.mechanism.find("CYP1A2 inhibition") != std::string::npos;
    });
    REQUIRE(flag != r.flags.end());
    REQUIRE(has(flag->citation, "FDA"));
    REQUIRE(has(flag->citation, "public domain"));
    REQUIRE(has(flag->direction, "perpetrator"));
    REQUIRE_FALSE(flag->boundaryNote.empty());
    REQUIRE(has(flag->boundaryNote, "not advice"));
    REQUIRE(has(flag->boundaryNote, "no severity score"));

    // An unknown member is listed, never silently ignored.
    REQUIRE(r.unknownMembers.size() == 1);
    REQUIRE(r.unknownMembers[0] == "unobtainium-xr");
    REQUIRE(std::find(r.members.begin(), r.members.end(), "unobtainium-xr") != r.members.end());
}

TEST_CASE("Hand-curated supplement and mechanism-class flags carry their own citations",
          "[mechanism]") {
    RealMechanism m(nullptr, {}, loaded());
    const StackReport r =
        m.checkStack({"grapefruit-juice", "simvastatin", "syrian-rue", "5-htp"});
    REQUIRE(r.unknownMembers.empty());

    const bool grapefruit =
        std::any_of(r.flags.begin(), r.flags.end(), [](const InteractionFlag& f) {
            return f.leftId == "grapefruit-juice" &&
                   f.mechanism.find("CYP3A4 inhibition") != std::string::npos &&
                   f.citation.find("furanocoumarin") != std::string::npos;
        });
    REQUIRE(grapefruit);

    const bool maoi = std::any_of(r.flags.begin(), r.flags.end(), [](const InteractionFlag& f) {
        return f.leftId == "syrian-rue" &&
               f.mechanism.find("MAO inhibition") != std::string::npos;
    });
    REQUIRE(maoi);

    for (const InteractionFlag& f : r.flags) {
        REQUIRE_FALSE(f.citation.empty());
        REQUIRE_FALSE(f.boundaryNote.empty());
    }
}

TEST_CASE("Pharmacogenomic notes are conditional and carry their activity-score band",
          "[mechanism]") {
    RealMechanism m(nullptr, {}, loaded());
    const PharmacogenomicReport r = m.pharmacogenomics("codeine");
    // The CPIC content is a bundled CC0 pack, so no live query is ever claimed.
    REQUIRE_FALSE(r.networkAvailable);
    REQUIRE_FALSE(r.notes.empty());
    REQUIRE_FALSE(r.boundaryStatement.empty());
    REQUIRE(has(r.boundaryStatement, "no dose"));
    for (const PharmacogenomicNote& n : r.notes) {
        REQUIRE_FALSE(n.citation.empty());
        REQUIRE_FALSE(n.activityScoreBand.empty());
        REQUIRE_FALSE(has(lower(n.implication), "extensive metabolizer"));
    }
    const PharmacogenomicReport none = m.pharmacogenomics("unobtainium-xr");
    REQUIRE(none.notes.empty());
    REQUIRE_FALSE(none.warnings.empty());
}

TEST_CASE("Retrieval reports state whether anything was retrieved, and why not", "[mechanism]") {
    RealBackend backend;
    Services s = backend.services();
    const MechanismReport mech = s.mechanism->mechanisms("caffeine");
    const PathwayContext path = s.mechanism->pathways("P29274");

    REQUIRE(mech.networkAvailable == mechanismNetworkAvailable());
    REQUIRE(path.networkAvailable == mechanismNetworkAvailable());
    REQUIRE_FALSE(mech.coverageNote.empty());
    REQUIRE_FALSE(mech.warnings.empty());
    // No pathway result, retrieved or not, may carry an impact score - and the report
    // says so every time rather than leaving its absence to be noticed.
    REQUIRE(std::any_of(path.warnings.begin(), path.warnings.end(), [](const std::string& w) {
        return w.find("no pathway impact score") != std::string::npos;
    }));

    if (!mechanismNetworkAvailable()) {
        // The offline contract: nothing retrieved, nothing invented, reason stated.
        REQUIRE(mech.entries.empty());
        REQUIRE_FALSE(mech.retrievalAttempted);
        REQUIRE(path.pathways.empty());
        REQUIRE(std::any_of(mech.warnings.begin(), mech.warnings.end(), [](const std::string& w) {
            return w.find("no network transport") != std::string::npos;
        }));
        REQUIRE(std::any_of(path.warnings.begin(), path.warnings.end(), [](const std::string& w) {
            return w.find("no network transport") != std::string::npos;
        }));
        REQUIRE(has(mech.coverageNote, "says nothing at all about the compound"));
    } else {
        // Online, every retrieved record is Measured and carries its source.
        for (const MechanismEntry& e : mech.entries) {
            REQUIRE(e.provenance == Provenance::Measured);
            REQUIRE(has(e.source, "ChEMBL"));
            // BioCAD's own modality axis is never inferred from an action type.
            REQUIRE(e.modality == InhibitionModality::Unknown);
        }
    }
}

TEST_CASE("A compound with no cross-reference yields a query statement, not a verdict",
          "[mechanism]") {
    RealMechanism m(nullptr, {}, loaded());
    const MechanismReport r = m.mechanisms("unobtainium-xr");
    REQUIRE(r.entries.empty());
    REQUIRE_FALSE(r.retrievalAttempted);
    REQUIRE(has(r.coverageNote, "says nothing at all about the compound"));
}
