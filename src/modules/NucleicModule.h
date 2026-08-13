// modules/NucleicModule.h - the DNA/RNA workbench adapter.
//
// WHAT THIS FILE IS FOR. bio/NucSeq.*, bio/NucIo.*, bio/Restriction.*,
// bio/OligoThermo.* and bio/Codon.* know sequence arithmetic and nearest-neighbour
// thermodynamics but know nothing about a workflow. This adapter is the single
// place where those pieces are composed into the two things a user actually asks
// for and neither of them provides on its own: a PCR primer pair, and a CRISPR
// guide with an honest off-target scope. Everything else here delegates - there is
// exactly one parser, one digest, one Tm.
//
// WHY THE GUIDE SEARCH LOOKS PARANOID. An off-target count is a count *inside a
// reference*, and the reference is whatever the caller handed over. Two bases of
// context change the meaning of "0 off-targets" from a useful screen to a
// falsehood: 0 in 2.7 kb of plasmid says nothing about a 3.1 Gb genome. So
// GuideSearchResult carries basesSearched, referenceName and a prose
// scopeStatement, and genomeWideClaimPossible is a deliberately conservative gate
// that is false by default and false for every input smaller than any real genome.
// No surface in BioCAD may render an off-target count without that scope, and the
// agent tool refuses genome-wide specificity questions outright.
//
// BIOSECURITY BOUNDARY, permanent and repeated here because this is the file a
// reader lands on: no synthesis-vendor integration, no order-sheet export, no
// pathogen-driven batch gene design, no therapeutic or germline CRISPR framing,
// and no "synthesise this" affordance. Export is FASTA and GenBank only. Codon
// optimization is constraint satisfaction, never an expression prediction.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "contracts/IModules.h"

namespace biocad {

// Primer-design limits. These are hard rejections, not weights: a pair that
// violates one is not returned at any rank, because ranking a primer that will
// dimerise with itself is worse than returning nothing.
struct PrimerDesignLimits {
    int    minLength = 18;
    int    maxLength = 30;
    // How far the primer 5' ends may be pushed outside the requested interval.
    // The product always CONTAINS the requested interval; it is never trimmed.
    int    maxFlank = 18;
    double minGcPercent = 40.0;
    double maxGcPercent = 60.0;
    // Each primer's Tm must sit inside this half-window around the target Tm.
    double tmWindowC = 4.0;
    // Pairs above this |Tm(F) - Tm(R)| are rejected: one primer would anneal
    // while the other does not, which no amount of ranking fixes.
    double maxTmDifferenceC = 2.0;
    // A self-dimer or hairpin at or below this dG37 (kcal/mol) rejects the primer.
    // -6.0 is the conventional working line for a 20-mer primer: it is roughly a
    // 3-4 base-pair duplex, and anything more stable competes with the template.
    double minSelfDeltaG37 = -6.0;
    // The two primers together. A cross-dimer costs both primers at once, so it
    // gets a slightly looser threshold applied to a strictly worse consequence.
    double minCrossDeltaG37 = -8.0;
    // Longest tolerated single-base run. Five identical bases slips.
    int    maxHomopolymer = 4;
    // At most this many G/C in the last five bases: a GC-saturated 3' end
    // primes mismatched templates.
    int    maxTerminalGc = 3;
    std::size_t maxPairs = 5;
};

// Guide-search limits and the scope gate.
struct GuideSearchLimits {
    int    protospacerLength = 20;
    double minGcPercent = 25.0;
    double maxGcPercent = 80.0;
    int    maxMismatches = 2;      // the pigeonhole seeding below assumes <= 2
    std::size_t maxGuides = 200;
    // A reference smaller than this cannot be a genome: the smallest free-living
    // bacterial genomes are ~0.6 Mb and the smallest endosymbiont ~0.11 Mb, so
    // 1 Mb is already generous as a floor while still excluding every plasmid,
    // amplicon and contig anyone will paste in.
    std::int64_t genomeWideMinBases = 1000000;
};

class RealNucleicAcid final : public INucleicAcidModule {
public:
    RealNucleicAcid() = default;
    explicit RealNucleicAcid(PrimerDesignLimits primer, GuideSearchLimits guide)
        : primer_(primer), guide_(guide) {}

    std::optional<NucRecord> parse(const std::string& text) const override;
    std::string              toFasta(const NucRecord& r) const override;
    std::string              toGenBank(const NucRecord& r) const override;
    std::string              reverseComplement(const std::string& seq) const override;

    TranslationResult translate(const NucRecord& r, int geneticCodeId,
                                int minOrfAminoAcids) const override;

    RestrictionDigest digest(const NucRecord& r,
                             const std::vector<std::string>& enzymes) const override;

    OligoThermo oligo(const std::string& seq, double naMolar, double mgMolar, double oligoMolar,
                      double dntpMolar) const override;

    std::vector<SecondaryStructure> selfStructures(const std::string& seq,
                                                   double naMolar) const override;

    std::vector<PrimerPair> designPrimers(const NucRecord& r, int begin, int end,
                                          double targetTmC) const override;

    CodonMetrics codonMetrics(const std::string& cds,
                              const std::string& usageTable) const override;

    CodonOptimizationResult optimizeCodons(
        const std::string& cds, const std::string& usageTable,
        const std::vector<std::string>& forbiddenSites) const override;

    GuideSearchResult findGuides(const NucRecord& target, const NucRecord& reference,
                                const std::string& pam) const override;

    [[nodiscard]] const PrimerDesignLimits& primerLimits() const { return primer_; }
    [[nodiscard]] const GuideSearchLimits&  guideLimits() const { return guide_; }

private:
    PrimerDesignLimits primer_{};
    GuideSearchLimits  guide_{};
};

// The fixed statement every surface renders verbatim, so the framing cannot drift
// between the panel, the agent tools and the docs.
const char* nucleicScopeNote();

}  // namespace biocad
