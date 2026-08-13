// Codon usage, CAI, and the constraint-based codon optimizer.
//
// The optimizer makes exactly two hard promises - the translation is preserved
// character for character, and no forbidden pattern occurs anywhere in the result,
// including across codon boundaries and on the complementary strand - so those
// are tested over many random inputs rather than one example. The CAI convention
// is pinned by the sequence built from each table's most frequent codons, whose
// CAI must be exactly 1.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <filesystem>
#include <random>
#include <string>

#include "bio/Codon.h"
#include "bio/NucSeq.h"
#include "core/Error.h"

using namespace biocad;
using namespace biocad::bio;

namespace {

std::filesystem::path nucleicPack(const char* name) {
    return std::filesystem::path(BIOCAD_ASSETS_DIR) / "packs" / "nucleic" / name;
}

const CodonUsageSet& usage() {
    static const CodonUsageSet set = loadCodonUsage(nucleicPack("codon-usage.json"));
    return set;
}

const GeneticCodeTables& codes() {
    static const GeneticCodeTables tables = loadGeneticCodes(nucleicPack("genetic-codes.json"));
    return tables;
}

}  // namespace

TEST_CASE("the codon usage pack carries complete, cited tables", "[bio][codon]") {
    REQUIRE(usage().find("ecoli-k12") != nullptr);
    REQUIRE(usage().find("hsapiens") != nullptr);
    for (const char* id : {"ecoli-k12", "hsapiens"}) {
        const CodonUsageTable& table = *usage().find(id);
        REQUIRE(table.entries().size() == 64);
        REQUIRE_FALSE(table.source().empty());
        // Every family must contain exactly one codon with w == 1, or the
        // relative-adaptiveness normalisation is broken.
        for (char aa : std::string("ACDEFGHIKLMNPQRSTVWY")) {
            int ones = 0, members = 0;
            for (const auto& e : table.entries()) {
                if (e.aminoAcid.empty() || e.aminoAcid[0] != aa) continue;
                ++members;
                if (e.relativeAdaptiveness == 1.0) ++ones;
            }
            REQUIRE(members >= 1);
            REQUIRE(ones == 1);
        }
    }
}

TEST_CASE("CAI of a most-frequent-codon-only sequence is exactly 1", "[bio][codon]") {
    for (const char* id : {"ecoli-k12", "hsapiens"}) {
        const CodonUsageTable& table = *usage().find(id);
        std::string best;
        for (char aa : std::string("ACDEFGHIKLNPQRSTVY")) best += table.mostFrequentCodon(aa);
        const CodonMetrics metrics = codonMetrics(best, table);
        REQUIRE(metrics.cai.value == 1.0);
        REQUIRE(metrics.cai.provenance == Provenance::Measured);
        REQUIRE(metrics.cai.unit.empty());
        REQUIRE(metrics.usageTableName == table.name());
        // Met and Trp are excluded by the stated convention, so adding them cannot
        // move the score.
        const CodonMetrics withMW = codonMetrics(best + table.mostFrequentCodon('M') +
                                                     table.mostFrequentCodon('W'),
                                                 table);
        REQUIRE(withMW.cai.value == 1.0);
        REQUIRE_FALSE(withMW.warnings.empty());
        // A rare-codon-only sequence must score well below 1.
        std::string worst;
        for (char aa : std::string("ACDEFGHIKLNPQRSTVY")) {
            const CodonUsageEntry* pick = nullptr;
            for (const auto& e : table.entries()) {
                if (e.aminoAcid.empty() || e.aminoAcid[0] != aa) continue;
                if (e.relativeAdaptiveness <= 0.0) continue;
                if (!pick || e.relativeAdaptiveness < pick->relativeAdaptiveness) pick = &e;
            }
            REQUIRE(pick != nullptr);
            worst += pick->codon;
        }
        REQUIRE(codonMetrics(worst, table).cai.value < 0.5);
    }
}

TEST_CASE("GC and GC3 are measured, and a short input is notComputed", "[bio][codon]") {
    const CodonUsageTable& table = *usage().find("ecoli-k12");
    const CodonMetrics m = codonMetrics("ATGGGCCCCATG", table);
    REQUIRE_THAT(m.gcPercent.value, Catch::Matchers::WithinAbs(gcPercent("ATGGGCCCCATG"), 1e-12));
    REQUIRE_THAT(m.gc3Percent.value, Catch::Matchers::WithinAbs(gc3Percent("ATGGGCCCCATG"), 1e-12));
    const CodonMetrics tooShort = codonMetrics("AT", table);
    REQUIRE(tooShort.cai.provenance == Provenance::NotComputed);
}

TEST_CASE("findForbidden sees sites across codon boundaries and on both strands",
          "[bio][codon]") {
    // GAATTC placed so it straddles two codons.
    REQUIRE(findForbidden("ATGGAATTCAAA", {"GAATTC"}, 0).size() == 1);
    // A non-palindromic pattern must be caught on the complementary strand too.
    REQUIRE(findForbidden("ATGCCCTTTGGG", {"CCCAAA"}, 0).size() == 1);
    // Degenerate patterns match by IUPAC consistency.
    REQUIRE(findForbidden("ATGGACGTCATG", {"GRCGYC"}, 0).size() >= 1);
    // Homopolymer limits are expressed as patterns: a run of six A is forbidden at
    // maxHomopolymer 5 and permitted at 6.
    const std::string sixA = "CTG" + std::string(6, 'A') + "CTG";
    REQUIRE(findForbidden(sixA, {}, 5).size() >= 1);
    REQUIRE(findForbidden(sixA, {}, 6).empty());
    REQUIRE(findForbidden("ATGGGGCCCATG", {"GRCGYC"}, 0).empty());
    REQUIRE_THROWS_AS(findForbidden("ACGT", {"XYZ"}, 0), Error);
}

TEST_CASE("the optimizer preserves translation and avoids every forbidden site",
          "[bio][codon]") {
    const CodonUsageTable& table = *usage().find("ecoli-k12");
    const GeneticCode& code = *codes().find(11);
    OptimizeOptions options;
    options.geneticCodeId = 11;
    options.forbiddenPatterns = {"GAATTC", "AAGCTT", "GGATCC", "CTGCAG", "GTCGAC",
                                 "CCCGGG", "GCGGCCGC", "GGTACC", "TCTAGA", "GRCGYC"};
    options.maxHomopolymer = 5;

    const std::string alphabet = "ACDEFGHIKLMNPQRSTVWY";
    std::mt19937 rng(4242);
    std::uniform_int_distribution<int> lenDist(30, 160);
    std::uniform_int_distribution<int> aaDist(0, static_cast<int>(alphabet.size()) - 1);
    for (int i = 0; i < 200; ++i) {
        std::string protein;
        const int n = lenDist(rng);
        for (int k = 0; k < n; ++k) {
            protein.push_back(alphabet[static_cast<std::size_t>(aaDist(rng))]);
        }
        const CodonOptimizationResult r = optimizeCodons(protein, table, options);
        REQUIRE(r.remainingViolations.empty());
        REQUIRE(r.translationPreserved);
        REQUIRE(r.optimized.size() == protein.size() * 3);
        REQUIRE(translate(r.optimized, code) == protein);
        REQUIRE(findForbidden(r.optimized, options.forbiddenPatterns,
                              options.maxHomopolymer).empty());
        REQUIRE(r.after.cai.value > 0.0);
    }
}

TEST_CASE("the optimizer reports infeasibility instead of breaking a constraint",
          "[bio][codon]") {
    const CodonUsageTable& table = *usage().find("ecoli-k12");
    OptimizeOptions impossible;
    impossible.geneticCodeId = 1;
    impossible.forbiddenPatterns = {"TGG"};   // the only Trp codon
    impossible.maxHomopolymer = 0;
    const CodonOptimizationResult r = optimizeCodons("MWM", table, impossible);
    REQUIRE(r.optimized.empty());
    REQUIRE_FALSE(r.translationPreserved);
    REQUIRE_FALSE(r.remainingViolations.empty());

    // An unknown residue is reported per position, not silently dropped.
    OptimizeOptions plain;
    plain.geneticCodeId = 1;
    plain.maxHomopolymer = 0;
    const CodonOptimizationResult bad = optimizeCodons("MZM", table, plain);
    REQUIRE(bad.optimized.empty());
    REQUIRE_FALSE(bad.remainingViolations.empty());
    REQUIRE(optimizeCodons("", table, plain).remainingViolations.size() == 1);

    // A GC window that cannot be reached is reported as a violation; nothing is
    // returned that breaks the hard constraints.
    OptimizeOptions narrow;
    narrow.geneticCodeId = 11;
    narrow.maxHomopolymer = 0;
    narrow.enforceGcWindow = true;
    narrow.minGcPercent = 10.0;
    narrow.maxGcPercent = 20.0;
    const CodonOptimizationResult impossibleGc = optimizeCodons("MKKAAWLLEEGGGG", table, narrow);
    REQUIRE_FALSE(impossibleGc.remainingViolations.empty());
}

TEST_CASE("a reachable GC window is satisfied and reported", "[bio][codon]") {
    const CodonUsageTable& table = *usage().find("ecoli-k12");
    const GeneticCode& code = *codes().find(11);
    const std::string protein = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVK";
    OptimizeOptions options;
    options.geneticCodeId = 11;
    options.maxHomopolymer = 6;
    options.enforceGcWindow = true;
    options.minGcPercent = 45.0;
    options.maxGcPercent = 55.0;
    const CodonOptimizationResult r = optimizeCodons(protein, table, options);
    REQUIRE(r.remainingViolations.empty());
    REQUIRE(translate(r.optimized, code) == protein);
    REQUIRE(r.after.gcPercent.value >= 45.0);
    REQUIRE(r.after.gcPercent.value <= 55.0);
    // The result documents what it did, and never claims an expression outcome.
    REQUIRE_FALSE(r.assumptions.empty());
    bool saysConstraintSatisfaction = false;
    for (const auto& a : r.assumptions) {
        if (a.find("constraint satisfaction") != std::string::npos) saysConstraintSatisfaction = true;
    }
    REQUIRE(saysConstraintSatisfaction);
}
