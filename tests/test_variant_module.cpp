// tests/test_variant_module.cpp - IVariantModule end to end: what it computes, and
// more importantly what it refuses.
//
// The refusals are the contract. Below the homolog minimum every conservation-derived
// quantity must be NotComputed and must NAME the shortfall, because the failure mode
// of this whole area is a plausible-looking 0.03 computed from five sequences. The
// BLOSUM62 delta is the one number that survives, because it is a table lookup that
// does not depend on the homolog set at all.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <string>
#include <vector>

#include "bio/Conservation.h"
#include "bio/PdbReader.h"
#include "bio/Rotamer.h"
#include "contracts/Services.h"
#include "modules/VariantModule.h"

using namespace biocad;
using Catch::Approx;

namespace {

template <class T>
concept HasEnergy = requires(T t) { t.energy; };
template <class T>
concept HasScore = requires(T t) { t.score; };
template <class T>
concept HasDeltaDeltaG = requires(T t) { t.deltaDeltaG; };

// A homolog set large enough to be usable: 20 sequences over a 12-residue query,
// each differing from the query at a different position, so the clustering does not
// collapse them and position 1 stays perfectly conserved.
std::vector<std::string> usableHomologs(const std::string& query) {
    std::vector<std::string> out;
    for (std::size_t i = 0; i < query.size(); ++i) {
        for (char sub : {'G', 'K'}) {
            std::string h = query;
            h[i] = sub;
            out.push_back(h);
        }
    }
    return out;
}

}  // namespace

TEST_CASE("Below the homolog minimum every dependent score is NotComputed",
          "[variants][module]") {
    const RealVariants v;
    const std::string query = "MKTAYIAKQRQ";
    const ConservationProfile p = v.conservation(query, {query, query, query});
    REQUIRE_FALSE(p.usable);
    REQUIRE(p.homologs.sequenceCount == 3);
    REQUIRE(p.minimumHomologsRequired == bio::kMinimumHomologs);

    const VariantScore s = v.score(p, 1, 'P');
    CHECK(s.siftScore.provenance == Provenance::NotComputed);
    CHECK(s.proveanScore.provenance == Provenance::NotComputed);
    CHECK(s.columnEntropy.provenance == Provenance::NotComputed);
    // The reason names the number that was missing, not just "unavailable".
    CHECK(s.siftScore.source.find("15 homologs") != std::string::npos);
    CHECK(s.siftScore.source.find("3 supplied") != std::string::npos);
    CHECK(s.interpretation.find("not a pathogenicity call") != std::string::npos);
    // The published thresholds travel with the DTO whether or not it was scored.
    CHECK(s.siftDeleteriousBelow == 0.05);
    CHECK(s.proveanDeleteriousBelow == -2.282);
}

TEST_CASE("A usable profile scores, and says what the score is not",
          "[variants][module]") {
    const RealVariants v;
    const std::string query = "MKTAYIAKQRQ";
    const ConservationProfile p = v.conservation(query, usableHomologs(query));
    REQUIRE(p.usable);
    REQUIRE(p.homologs.sequenceCount == 22);
    REQUIRE(p.columns.size() == query.size());

    const VariantScore s = v.score(p, 1, 'P');
    CHECK(s.wildType == 'M');
    CHECK(s.mutant == 'P');
    // BLOSUM62: score(M,P) = -2, score(M,M) = 5, so the delta is -7 exactly.
    CHECK(s.blosum62Delta.provenance == Provenance::Measured);
    CHECK(s.blosum62Delta.value == Approx(-7.0));
    CHECK(s.blosum62Delta.unit == "half-bits");
    // Conservation scores rank; they carry no unit, which makeQuantity enforces.
    CHECK(s.siftScore.provenance == Provenance::Heuristic);
    CHECK(s.siftScore.unit.empty());
    CHECK(s.proveanScore.provenance == Provenance::Heuristic);
    CHECK(s.proveanScore.unit.empty());
    CHECK(s.columnEntropy.unit == "bits");
    // The alignment the score came from travels with it.
    CHECK(s.homologs.sequenceCount == 22);
    CHECK(s.homologs.medianIdentityPct > 80.0);
    CHECK(s.interpretation.find("not a pathogenicity call") != std::string::npos);
    CHECK(s.siftScore.source.find("median identity") != std::string::npos);

    // The wild type is the consensus of a mostly-conserved column, so its own SIFT
    // score is exactly 1.
    CHECK(bio::siftScore(p.columns[0], 'M') == 1.0);
}

TEST_CASE("A position outside the profile is refused by name", "[variants][module]") {
    const RealVariants v;
    const std::string query = "MKTAYIAKQRQ";
    const ConservationProfile p = v.conservation(query, usableHomologs(query));
    const VariantScore s = v.score(p, 999, 'A');
    CHECK(s.siftScore.provenance == Provenance::NotComputed);
    CHECK(s.siftScore.source.find("999") != std::string::npos);
    CHECK(s.blosum62Delta.provenance == Provenance::NotComputed);
}

TEST_CASE("A rebuilt side chain is Model, names its library and carries no energy",
          "[variants][module]") {
    const RealVariants v;
    const bio::Structure s =
        bio::readPdbFile(std::filesystem::path(BIOCAD_TEST_FIXTURES) / "1VFB.pdb");
    REQUIRE(s.model(1) != nullptr);
    const std::string chain = s.model(1)->chains.front().id;
    const int resNum = s.model(1)->chains.front().residues[10].authSeqId;

    const RotamerRebuild r = v.rebuild(s, chain, resNum, 'L');
    CHECK(r.provenance == Provenance::Model);
    CHECK(r.mutant == 'L');
    CHECK(r.position == resNum);
    CHECK_FALSE(r.rotamerLibrarySource.empty());
    CHECK(r.rotamerLibrarySource.find("PDB-derived") != std::string::npos);
    CHECK(r.chiAngles.size() == 2u);          // leucine has chi1 and chi2
    CHECK(r.rotamerProbability > 0.0);
    CHECK(r.rotamerProbability <= 1.0);
    CHECK(r.clashCount >= 0);
    CHECK_FALSE(r.assumptions.empty());
    bool saysNoEnergy = false, saysBackboneFixed = false, saysDee = false;
    for (const auto& a : r.assumptions) {
        if (a.find("no force field") != std::string::npos) saysNoEnergy = true;
        if (a.find("backbone is held fixed") != std::string::npos) saysBackboneFixed = true;
        if (a.find("dead-end elimination") != std::string::npos) saysDee = true;
    }
    CHECK(saysNoEnergy);
    CHECK(saysBackboneFixed);
    CHECK(saysDee);
    // RotamerRebuild has no energy field at all - the rebuild is a construction, not
    // a calculation. Asserted at compile time so that adding one breaks the build
    // rather than passing review.
    static_assert(!HasEnergy<RotamerRebuild>);
    static_assert(!HasScore<RotamerRebuild>);
    static_assert(!HasDeltaDeltaG<RotamerRebuild>);
    // ...while the type that IS allowed to carry a ddG does, and is not produced here.
    static_assert(HasDeltaDeltaG<StabilityPrediction>);
}

TEST_CASE("Proline and an unknown residue are refused with a reason",
          "[variants][module]") {
    const RealVariants v;
    const bio::Structure s =
        bio::readPdbFile(std::filesystem::path(BIOCAD_TEST_FIXTURES) / "1VFB.pdb");
    const std::string chain = s.model(1)->chains.front().id;
    const int resNum = s.model(1)->chains.front().residues[10].authSeqId;

    const RotamerRebuild pro = v.rebuild(s, chain, resNum, 'P');
    CHECK(pro.chiAngles.empty());
    REQUIRE_FALSE(pro.warnings.empty());
    CHECK(pro.warnings[0].find("Proline") != std::string::npos);

    const RotamerRebuild bad = v.rebuild(s, chain, resNum, 'Z');
    CHECK(bad.chiAngles.empty());
    REQUIRE_FALSE(bad.warnings.empty());

    const RotamerRebuild missing = v.rebuild(s, "ZZ", 1, 'L');
    REQUIRE_FALSE(missing.warnings.empty());
    CHECK(missing.warnings[0].find("not in this structure") != std::string::npos);

    // Glycine is not an error: it has no side chain, and the result says so.
    const RotamerRebuild gly = v.rebuild(s, chain, resNum, 'G');
    CHECK(gly.warnings.empty());
    CHECK(gly.chiAngles.empty());
    REQUIRE_FALSE(gly.assumptions.empty());
    CHECK(gly.assumptions[0].find("no side chain") != std::string::npos);
}

TEST_CASE("The rebuild repacks a real neighbourhood", "[variants][module]") {
    const RealVariants v;
    const bio::Structure s =
        bio::readPdbFile(std::filesystem::path(BIOCAD_TEST_FIXTURES) / "1VFB.pdb");
    const std::string chain = s.model(1)->chains.front().id;
    const int resNum = s.model(1)->chains.front().residues[30].authSeqId;

    const RotamerRebuild r = v.rebuild(s, chain, resNum, 'W');
    CHECK(r.provenance == Provenance::Model);
    CHECK_FALSE(r.repackedNeighbours.empty());
    CHECK(r.repackedNeighbours.size() <= static_cast<std::size_t>(kMaxRepackedNeighbours));
    for (const auto& n : r.repackedNeighbours) CHECK(n.find(':') != std::string::npos);
}

TEST_CASE("Services::valid() requires the variant module", "[variants][module]") {
    RealVariants variants;
    Services s;
    // valid() only tests for null, so distinct non-null addresses exercise the exact
    // startup predicate without dragging in the Windows-only modules.
    char sentinel = 0;
    auto* p = &sentinel;
    s.library = reinterpret_cast<ILibrary*>(p);
    s.stability = reinterpret_cast<IStabilityModule*>(p);
    s.admet = reinterpret_cast<IAdmetModule*>(p);
    s.absorption = reinterpret_cast<IAbsorptionModule*>(p);
    s.similarity = reinterpret_cast<ISimilarityModule*>(p);
    s.legal = reinterpret_cast<ILegalModule*>(p);
    s.docking = reinterpret_cast<IDockingModule*>(p);
    s.runs = reinterpret_cast<IRunStore*>(p);
    s.pharmacodynamics = reinterpret_cast<IPharmacodynamicsModule*>(p);
    s.sequence = reinterpret_cast<ISequenceModule*>(p);
    s.structure = reinterpret_cast<IStructureModule*>(p);
    s.metabolismFacts = reinterpret_cast<IMetabolismFactsModule*>(p);
    s.alerts = reinterpret_cast<IAlertsModule*>(p);
    s.ionization = reinterpret_cast<IIonizationModule*>(p);
    s.nucleicAcid = reinterpret_cast<INucleicAcidModule*>(p);
    s.assay = reinterpret_cast<IAssayModule*>(p);
    s.biologics = reinterpret_cast<IBiologicsModule*>(p);
    s.populationPk = reinterpret_cast<IPopulationPkModule*>(p);
    s.simulation = reinterpret_cast<ISimulationModule*>(p);
    s.enrichment = reinterpret_cast<IEnrichmentModule*>(p);
    s.mechanism = reinterpret_cast<IMechanismModule*>(p);
    CHECK_FALSE(s.valid());
    s.variants = &variants;
    CHECK(s.valid());
}
