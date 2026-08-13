// tests/test_nucleic_module.cpp - the DNA/RNA workbench adapter.
//
// The two things asserted hardest here are the two things this module owns rather
// than delegates: primer design, where a violated threshold must make a pair
// ABSENT rather than low-ranked, and the guide search, where every count must
// arrive with the scope it was counted in. A guide test that checks a count and
// not its scope would pass while the product lied.
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "contracts/Services.h"
#include "modules/NucleicModule.h"

using namespace biocad;

namespace {

// A deterministic 900 bp template with no long homopolymer, so primer design has
// somewhere to work.
std::string testSequence(std::size_t n) {
    static const char* unit =
        "GATCTAGCCATGGTACGATCGGTACCATTGCAAGCTTGACGTCTGCAGGAATTCCGTAAGCT"
        "GGATCCTTAGCAACGTTGCCATCGTAGCTAGGCATCGATTGCCAAGGTCATGCGATCAGCTT";
    std::string s;
    while (s.size() < n) s += unit;
    s.resize(n);
    return s;
}

NucRecord record(std::string id, std::string seq, bool circular = false) {
    NucRecord r;
    r.id = std::move(id);
    r.sequence = std::move(seq);
    r.circular = circular;
    return r;
}

// The real 2686 bp pUC19 (L09137.2). A synthetic repeated sequence is useless for a
// guide test: every protospacer in it has dozens of exact off-targets, so the test
// would be measuring the repeat rather than the search.
std::string pUC19Text() {
    std::ifstream in(std::filesystem::path(BIOCAD_TEST_FIXTURES) / "pUC19.gb");
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

}  // namespace

TEST_CASE("Services::valid() requires the nucleic-acid module", "[nucleic]") {
    // Only null-ness is inspected here, never a dereference: the point is that a
    // member left out of valid() is a crash waiting for the first panel to use it.
    const RealNucleicAcid mod;
    Services s;
    auto* p = reinterpret_cast<char*>(1);
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
    s.assay = reinterpret_cast<IAssayModule*>(p);
    REQUIRE_FALSE(s.valid());
    s.nucleicAcid = const_cast<RealNucleicAcid*>(&mod);
    REQUIRE(s.valid());
}

TEST_CASE("FASTA and GenBank inputs round-trip through the same record", "[nucleic]") {
    const RealNucleicAcid mod;
    const std::string seq = testSequence(240);
    const NucRecord src = record("rt1", seq);

    const auto fromFasta = mod.parse(mod.toFasta(src));
    const auto fromGenBank = mod.parse(mod.toGenBank(src));
    REQUIRE(fromFasta.has_value());
    REQUIRE(fromGenBank.has_value());
    REQUIRE(fromFasta->sequence == seq);
    REQUIRE(fromGenBank->sequence == seq);
    // Structural identity: re-exporting either parse yields the same FASTA, which
    // is what makes GenBank -> FASTA -> GenBank safe to do to a user's file.
    REQUIRE(mod.toFasta(*fromFasta) == mod.toFasta(*fromGenBank));
    // Unparseable text is std::nullopt, not an exception escaping into a UI frame.
    REQUIRE(mod.parse("not a sequence at all") == std::nullopt);
}

TEST_CASE("Designed primers satisfy every stated threshold", "[nucleic]") {
    const RealNucleicAcid mod;
    const NucRecord r = record("amp", testSequence(900));
    const std::vector<PrimerPair> pairs = mod.designPrimers(r, 100, 600, 60.0);
    REQUIRE_FALSE(pairs.empty());
    const PrimerDesignLimits lim = mod.primerLimits();
    for (const auto& p : pairs) {
        REQUIRE(p.tmDifference < lim.maxTmDifferenceC);
        REQUIRE(p.forwardOligo.gcPercent >= lim.minGcPercent);
        REQUIRE(p.forwardOligo.gcPercent <= lim.maxGcPercent);
        REQUIRE(p.reverseOligo.gcPercent >= lim.minGcPercent);
        REQUIRE(p.reverseOligo.gcPercent <= lim.maxGcPercent);
        // The product CONTAINS the requested interval: the flank is the freedom,
        // the interval is the contract.
        REQUIRE(p.productBegin <= 100);
        REQUIRE(p.productEnd >= 600);
        // Nothing at or below the rejection threshold survived into the list.
        for (const auto& l : p.liabilities) REQUIRE(l.deltaG37.value > lim.minCrossDeltaG37);
        REQUIRE_FALSE(p.warnings.empty());  // the thresholds travel with the pair
        // A 3' GC clamp on both primers.
        REQUIRE((p.forwardOligo.sequence.back() == 'G' || p.forwardOligo.sequence.back() == 'C'));
        REQUIRE((p.reverseOligo.sequence.back() == 'G' || p.reverseOligo.sequence.back() == 'C'));
    }
}

TEST_CASE("A guide search reports its scope, not just its counts", "[nucleic]") {
    const RealNucleicAcid mod;
    const auto plasmid = mod.parse(pUC19Text());
    REQUIRE(plasmid.has_value());
    REQUIRE(plasmid->sequence.size() == 2686);
    REQUIRE(plasmid->circular);

    const GuideSearchResult g = mod.findGuides(*plasmid, *plasmid, "NGG");
    REQUIRE(g.basesSearched == 2686);
    REQUIRE(g.referenceName == "pUC19");
    REQUIRE_FALSE(g.genomeWideClaimPossible);
    REQUIRE(g.scopeStatement.find("pUC19") != std::string::npos);
    REQUIRE(g.scopeStatement.find("2686") != std::string::npos);
    REQUIRE(g.scopeStatement.find("NOT a genome-wide specificity claim") != std::string::npos);
    REQUIRE_FALSE(g.guides.empty());
    // The enumeration cap is itself a scope limit and must be stated when hit.
    if (g.guides.size() >= mod.guideLimits().maxGuides)
        REQUIRE(g.scopeStatement.find("enumeration stopped") != std::string::npos);
}

TEST_CASE("Off-targets are counted by mismatch class, excluding the on-target", "[nucleic]") {
    const RealNucleicAcid mod;
    const auto plasmid = mod.parse(pUC19Text());
    REQUIRE(plasmid.has_value());
    const std::string target = plasmid->sequence.substr(0, 400);

    const GuideSearchResult base = mod.findGuides(record("t", target), record("t", target), "NGG");
    REQUIRE_FALSE(base.guides.empty());

    // Take a guide with no off-target in the bare target, then plant an exact copy
    // and a 2-mismatch copy in a longer reference.
    const GuideCandidate* clean = nullptr;
    for (const auto& c : base.guides)
        if (c.exactOffTargets == 0 && c.oneMismatchOffTargets == 0 &&
            c.twoMismatchOffTargets == 0) {
            clean = &c;
            break;
        }
    REQUIRE(clean != nullptr);

    std::string twoMismatch = clean->protospacer;
    twoMismatch[3] = (twoMismatch[3] == 'A') ? 'C' : 'A';
    twoMismatch[12] = (twoMismatch[12] == 'G') ? 'T' : 'G';
    REQUIRE(twoMismatch != clean->protospacer);

    const NucRecord ref = record("ref-with-plants",
                                 target + "AAAAAAAAAA" + clean->protospacer + clean->pam +
                                     "AAAAAAAAAA" + twoMismatch + clean->pam + "AAAAAAAAAA");
    const GuideSearchResult g = mod.findGuides(record("t", target), ref, "NGG");

    const GuideCandidate* found = nullptr;
    for (const auto& c : g.guides)
        if (c.protospacer == clean->protospacer && c.strand == clean->strand) found = &c;
    REQUIRE(found != nullptr);
    REQUIRE(found->exactOffTargets == 1);        // the planted duplicate, not itself
    REQUIRE(found->oneMismatchOffTargets == 0);
    REQUIRE(found->twoMismatchOffTargets == 1);  // counted in its own bucket
    REQUIRE(g.basesSearched == static_cast<std::int64_t>(ref.sequence.size()));
    // Still not a genome: a reference that merely contains the target does not
    // license a genome-wide claim.
    REQUIRE_FALSE(g.genomeWideClaimPossible);
}

TEST_CASE("A non-NGG PAM is honoured under IUPAC", "[nucleic]") {
    const RealNucleicAcid mod;
    const auto plasmid = mod.parse(pUC19Text());
    REQUIRE(plasmid.has_value());

    const GuideSearchResult ngg = mod.findGuides(*plasmid, *plasmid, "NGG");
    const GuideSearchResult tttv = mod.findGuides(*plasmid, *plasmid, "TTTV");
    REQUIRE_FALSE(ngg.guides.empty());
    REQUIRE_FALSE(tttv.guides.empty());
    for (const auto& c : ngg.guides) {
        REQUIRE(c.pam.size() == 3);
        REQUIRE(c.pam[1] == 'G');
        REQUIRE(c.pam[2] == 'G');
        REQUIRE(c.protospacer.size() == 20);
    }
    for (const auto& c : tttv.guides) {
        REQUIRE(c.pam.size() == 4);
        REQUIRE(c.pam.substr(0, 3) == "TTT");
        REQUIRE(c.pam[3] != 'T');  // V = A, C or G
    }
}

TEST_CASE("Oligo thermodynamics carry their conditions and their parameter set",
          "[nucleic]") {
    const RealNucleicAcid mod;
    const OligoThermo t = mod.oligo("GTAAAACGACGGCCAGT", 0.05, 0.0, 2.5e-7, 0.0);
    REQUIRE(t.tm.unit == "degC");
    REQUIRE(t.tm.provenance == Provenance::Predicted);
    REQUIRE_FALSE(t.tm.source.empty());     // the parameter set, always
    REQUIRE(t.naMolar == 0.05);             // echoed back, so the Tm is reproducible
    REQUIRE(t.oligoMolar == 2.5e-7);
    REQUIRE_FALSE(t.assumptions.empty());
    // A zero concentration means "the caller did not say", and the documented
    // default stands in rather than a Tm computed at 0 M salt.
    const OligoThermo d = mod.oligo("GTAAAACGACGGCCAGT", 0.0, 0.0, 0.0, 0.0);
    REQUIRE(d.naMolar > 0.0);
    REQUIRE(d.oligoMolar > 0.0);
}

TEST_CASE("Codon metrics refuse an unknown usage table by name", "[nucleic]") {
    const RealNucleicAcid mod;
    const CodonMetrics m = mod.codonMetrics("ATGGCTAGCTAA", "no-such-table");
    REQUIRE(m.cai.provenance == Provenance::NotComputed);
    REQUIRE(m.cai.source.find("no-such-table") != std::string::npos);
    REQUIRE_FALSE(m.warnings.empty());
}
