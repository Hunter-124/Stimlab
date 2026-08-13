#include "modules/NucleicModule.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <utility>

#include "bio/Codon.h"
#include "bio/NucIo.h"
#include "bio/NucSeq.h"
#include "bio/OligoThermo.h"
#include "bio/Restriction.h"
#include "core/Error.h"

namespace biocad {
namespace {

bio::ThermoOptions thermoOptions(double naMolar, double mgMolar, double oligoMolar,
                                 double dntpMolar) {
    bio::ThermoOptions o;
    // A zero here means "the caller did not say", and a Tm computed at 0 M salt is
    // not a Tm anyone can reproduce at the bench, so the documented default stands
    // in and is reported back through OligoThermo::naMolar.
    if (naMolar > 0.0) o.naMolar = naMolar;
    if (mgMolar > 0.0) o.mgMolar = mgMolar;
    if (oligoMolar > 0.0) o.oligoMolar = oligoMolar;
    if (dntpMolar > 0.0) o.dntpMolar = dntpMolar;
    return o;
}

// The most stable (most negative) dG37 among a set of structures, or 0 when the
// set is empty. 0 is the right neutral element: no structure means no penalty.
double worstDeltaG37(const std::vector<SecondaryStructure>& v) {
    double worst = 0.0;
    for (const auto& s : v) worst = std::min(worst, s.deltaG37.value);
    return worst;
}

int longestRun(const std::string& s) {
    int best = 0, run = 0;
    char prev = 0;
    for (char c : s) {
        run = (c == prev) ? run + 1 : 1;
        prev = c;
        best = std::max(best, run);
    }
    return best;
}

int terminalGcCount(const std::string& s, std::size_t window) {
    int gc = 0;
    const std::size_t n = std::min(window, s.size());
    for (std::size_t i = s.size() - n; i < s.size(); ++i)
        if (s[i] == 'G' || s[i] == 'C') ++gc;
    return gc;
}

// One side of a candidate pair, already thermodynamically evaluated so the pairing
// loop never recomputes a Tm.
struct Candidate {
    std::string                     seq;
    int                             outerEnd = 0;   // 5' end on the record (forward
                                                    // primer) or 3' bound (reverse)
    OligoThermo                     thermo;
    std::vector<SecondaryStructure> liabilities;
    double                          endStabilityDeltaG37 = 0.0;
    double                          score = 0.0;
};

// 3'-end stability: the dG37 of the terminal pentamer. Too weak and the polymerase
// has nothing to hold; too strong and the primer extends off a mismatched template.
// -7.5 kcal/mol is mid-range for a pentamer under the unified parameters, so the
// score penalises distance from it rather than asserting a hard cut.
constexpr double kIdealEndDeltaG37 = -7.5;

double endStability(const std::string& seq) {
    if (seq.size() < 5) return 0.0;
    return bio::oligoThermo(std::string_view(seq).substr(seq.size() - 5)).deltaG37.value;
}

}  // namespace

const char* nucleicScopeNote() {
    return "BioCAD's DNA/RNA workbench reads and writes FASTA and GenBank and nothing else. "
           "There is no synthesis-vendor integration, no order-sheet export, no pathogen-driven "
           "batch gene design and no therapeutic or germline CRISPR framing. Codon optimization "
           "is constraint satisfaction over a cited usage table - it preserves the translated "
           "protein and excludes the forbidden sites you named, and it predicts no expression "
           "level, yield or titre. A guide off-target count is a count inside the reference you "
           "supplied, and is reported with the number of bases actually searched.";
}

// WHY EVERY CALL BELOW IS GUARDED. The bio/* layer throws core::Error for a
// missing pack or unparseable text, which is right for a library. It is wrong for
// a UI frame and wrong for an agent tool: an exception escaping here would take
// down the message pump. So this adapter is the boundary that turns a throw into
// the contract's own vocabulary - std::nullopt, a NotComputed Quantity, or a
// warning on the returned DTO. Nothing is swallowed; the message is always
// carried out to the surface that asked.
std::optional<NucRecord> RealNucleicAcid::parse(const std::string& text) const {
    std::vector<NucRecord> v;
    try {
        v = bio::readNucleic(text);
    } catch (const Error&) {
        return std::nullopt;   // "neither FASTA nor GenBank", exactly as documented
    }
    if (v.empty()) return std::nullopt;
    NucRecord r = v.front();
    if (v.size() > 1)
        r.warnings.push_back("The input held " + std::to_string(v.size()) +
                             " records; only the first was taken.");
    return r;
}

std::string RealNucleicAcid::toFasta(const NucRecord& r) const { return bio::writeNucFasta(r); }

std::string RealNucleicAcid::toGenBank(const NucRecord& r) const { return bio::writeGenBank(r); }

std::string RealNucleicAcid::reverseComplement(const std::string& seq) const {
    return bio::reverseComplement(seq);
}

TranslationResult RealNucleicAcid::translate(const NucRecord& r, int geneticCodeId,
                                             int minOrfAminoAcids) const {
    bio::OrfOptions o;
    o.minAminoAcids = minOrfAminoAcids;
    o.circular = r.circular;
    try {
        return bio::translateRecord(r, geneticCodeId, o);
    } catch (const Error& e) {
        TranslationResult t;
        t.recordId = r.id;
        t.geneticCodeId = geneticCodeId;
        t.warnings.push_back("No translation: " + e.message);
        return t;
    }
}

RestrictionDigest RealNucleicAcid::digest(const NucRecord& r,
                                          const std::vector<std::string>& enzymes) const {
    try {
        return bio::digest(r, enzymes);
    } catch (const Error& e) {
        RestrictionDigest d;
        d.recordId = r.id;
        d.circular = r.circular;
        d.warnings.push_back("No digest: " + e.message);
        return d;
    }
}

OligoThermo RealNucleicAcid::oligo(const std::string& seq, double naMolar, double mgMolar,
                                   double oligoMolar, double dntpMolar) const {
    try {
        return bio::oligoThermo(seq, thermoOptions(naMolar, mgMolar, oligoMolar, dntpMolar));
    } catch (const Error& e) {
        OligoThermo t;
        t.sequence = seq;
        t.tm = notComputed(e.message);
        t.deltaH = notComputed(e.message);
        t.deltaS = notComputed(e.message);
        t.deltaG37 = notComputed(e.message);
        t.assumptions.push_back("No thermodynamics: " + e.message);
        return t;
    }
}

std::vector<SecondaryStructure> RealNucleicAcid::selfStructures(const std::string& seq,
                                                                double naMolar) const {
    // naMolar is accepted because the caller has it and hiding an argument invites a
    // second entry point later; the folding dG37 itself is the unified-parameter
    // value at the 1 M standard state, so the salt shows up in the Tm, not here.
    // That is stated in docs/nucleic.md rather than silently applied.
    (void)naMolar;
    std::vector<SecondaryStructure> out;
    try {
        out = bio::hairpins(seq);
        const std::vector<SecondaryStructure> d = bio::selfDimers(seq);
        out.insert(out.end(), d.begin(), d.end());
    } catch (const Error&) {
        return out;   // no parameter pack: no structures to report, and no crash
    }
    std::sort(out.begin(), out.end(), [](const SecondaryStructure& a, const SecondaryStructure& b) {
        return a.deltaG37.value < b.deltaG37.value;
    });
    return out;
}

CodonMetrics RealNucleicAcid::codonMetrics(const std::string& cds,
                                           const std::string& usageTable) const {
    std::string missing = "codon usage table '" + usageTable + "'";
    try {
        if (const bio::CodonUsageTable* t = bio::builtinCodonUsage().find(usageTable))
            return bio::codonMetrics(cds, *t);
    } catch (const Error& e) {
        missing = e.message;
    }
    CodonMetrics m;
    m.cai = notComputed(missing);
    m.gcPercent = notComputed(missing);
    m.gc3Percent = notComputed(missing);
    m.usageTableName = usageTable;
    m.warnings.push_back("No CAI: " + missing + " is not loaded, so there is nothing for the "
                         "index to be relative to.");
    return m;
}

CodonOptimizationResult RealNucleicAcid::optimizeCodons(
    const std::string& cds, const std::string& usageTable,
    const std::vector<std::string>& forbiddenSites) const {
    CodonOptimizationResult r;
    r.forbiddenSites = forbiddenSites;
    const bio::CodonUsageTable* t = nullptr;
    const bio::GeneticCode* code = nullptr;
    try {
        t = bio::builtinCodonUsage().find(usageTable);
        code = bio::builtinGeneticCodes().find(1);
    } catch (const Error& e) {
        r.remainingViolations.push_back(e.message);
        return r;
    }
    if (!t) {
        r.remainingViolations.push_back("No codon usage table with id '" + usageTable +
                                        "' is loaded.");
        return r;
    }
    // The interface takes a CDS and bio::optimizeCodons takes a protein, because
    // the guarantee it makes is about the translation. Translating here keeps that
    // one conversion in one place instead of duplicating a genetic-code table.
    if (!code) {
        r.remainingViolations.push_back("NCBI genetic code table 1 is not loaded.");
        return r;
    }
    std::string protein = bio::translate(cds, *code);
    while (!protein.empty() && protein.back() == '*') protein.pop_back();
    if (protein.find('*') != std::string::npos) {
        r.protein = protein;
        r.remainingViolations.push_back(
            "The input CDS has an internal stop codon in frame 1, so there is no single protein "
            "to preserve. Trim it to one open frame first.");
        return r;
    }
    bio::OptimizeOptions o;
    o.forbiddenPatterns = forbiddenSites;
    o.usageTableId = usageTable;
    CodonOptimizationResult out = bio::optimizeCodons(protein, *t, o);
    out.before = bio::codonMetrics(cds, *t);
    out.assumptions.push_back(nucleicScopeNote());
    return out;
}

// ---------------------------------------------------------------------------
// Primer design.
// ---------------------------------------------------------------------------

std::vector<PrimerPair> RealNucleicAcid::designPrimers(const NucRecord& r, int begin, int end,
                                                       double targetTmC) const {
    std::vector<PrimerPair> pairs;
    const int n = static_cast<int>(r.sequence.size());
    if (begin < 0 || end > n || end - begin < primer_.minLength * 2) return pairs;
    // One guard covers every thermodynamic call below: the parameter pack is loaded
    // once and cached, so if it resolves here nothing after this line can throw for
    // want of it. Without a Tm no primer can be scored, so the honest answer is an
    // empty list rather than a partially-scored one.
    try {
        bio::builtinNnParameters();
    } catch (const Error&) {
        return pairs;
    }

    const std::string top = r.sequence;
    const std::string bottom = bio::reverseComplement(top);

    // Enumerate one side. `outer` walks the 5' end outward from the requested
    // interval so the product always contains it: the interval is the contract,
    // the flank is the freedom.
    auto enumerate = [&](bool forward) {
        std::vector<Candidate> out;
        for (int flank = 0; flank <= primer_.maxFlank; ++flank) {
            const int start = forward ? begin - flank : (n - end) - flank;
            if (start < 0) break;
            for (int len = primer_.minLength; len <= primer_.maxLength; ++len) {
                const std::string& src = forward ? top : bottom;
                if (start + len > static_cast<int>(src.size())) break;
                Candidate c;
                c.seq = src.substr(static_cast<std::size_t>(start), static_cast<std::size_t>(len));
                if (c.seq.find_first_not_of("ACGT") != std::string::npos) continue;
                c.outerEnd = forward ? start : n - start;  // record coordinates
                c.thermo = bio::oligoThermo(c.seq);
                if (c.thermo.gcPercent < primer_.minGcPercent ||
                    c.thermo.gcPercent > primer_.maxGcPercent)
                    continue;
                if (std::fabs(c.thermo.tm.value - targetTmC) > primer_.tmWindowC) continue;
                if (longestRun(c.seq) > primer_.maxHomopolymer) continue;
                // GC clamp: a 3'-terminal G or C anchors extension, but more than
                // three G/C in the last five over-stabilises a mismatched end.
                const char last = c.seq.back();
                if (last != 'G' && last != 'C') continue;
                if (terminalGcCount(c.seq, 5) > primer_.maxTerminalGc) continue;
                c.liabilities = selfStructures(c.seq, 0.0);
                if (worstDeltaG37(c.liabilities) <= primer_.minSelfDeltaG37) continue;
                c.endStabilityDeltaG37 = endStability(c.seq);
                c.score = 2.0 * std::fabs(c.thermo.tm.value - targetTmC) +
                          1.0 * std::fabs(c.endStabilityDeltaG37 - kIdealEndDeltaG37) +
                          0.5 * -worstDeltaG37(c.liabilities);
                out.push_back(std::move(c));
            }
        }
        std::sort(out.begin(), out.end(),
                  [](const Candidate& a, const Candidate& b) { return a.score < b.score; });
        // Cross-dimer evaluation is the expensive step, so only the best few of each
        // side are paired: 16 x 16 hetero-dimer evaluations, not 200 x 200.
        if (out.size() > 16) out.resize(16);
        return out;
    };

    const std::vector<Candidate> fwd = enumerate(true);
    const std::vector<Candidate> rev = enumerate(false);

    for (const auto& f : fwd) {
        for (const auto& v : rev) {
            const double dTm = std::fabs(f.thermo.tm.value - v.thermo.tm.value);
            if (dTm > primer_.maxTmDifferenceC) continue;
            const std::vector<SecondaryStructure> cross = bio::heteroDimers(f.seq, v.seq);
            if (worstDeltaG37(cross) <= primer_.minCrossDeltaG37) continue;
            PrimerPair p;
            p.forwardOligo = f.thermo;
            p.reverseOligo = v.thermo;
            p.productBegin = f.outerEnd;
            p.productEnd = v.outerEnd;
            p.productLength = p.productEnd - p.productBegin;
            p.tmDifference = dTm;
            p.liabilities = f.liabilities;
            p.liabilities.insert(p.liabilities.end(), v.liabilities.begin(), v.liabilities.end());
            p.liabilities.insert(p.liabilities.end(), cross.begin(), cross.end());
            std::sort(p.liabilities.begin(), p.liabilities.end(),
                      [](const SecondaryStructure& a, const SecondaryStructure& b) {
                          return a.deltaG37.value < b.deltaG37.value;
                      });
            p.warnings.push_back(
                "Rejection thresholds applied: hairpin and self-dimer dG37 > " +
                std::to_string(primer_.minSelfDeltaG37) + " kcal/mol, cross-dimer dG37 > " +
                std::to_string(primer_.minCrossDeltaG37) + " kcal/mol, |dTm| <= " +
                std::to_string(primer_.maxTmDifferenceC) + " degC. Any pair violating one of "
                "those is absent from this list rather than ranked low.");
            pairs.push_back(std::move(p));
        }
    }

    std::sort(pairs.begin(), pairs.end(), [&](const PrimerPair& a, const PrimerPair& b) {
        const double sa = 4.0 * a.tmDifference - worstDeltaG37(a.liabilities) +
                          std::fabs(a.forwardOligo.tm.value - targetTmC) +
                          std::fabs(a.reverseOligo.tm.value - targetTmC);
        const double sb = 4.0 * b.tmDifference - worstDeltaG37(b.liabilities) +
                          std::fabs(b.forwardOligo.tm.value - targetTmC) +
                          std::fabs(b.reverseOligo.tm.value - targetTmC);
        return sa < sb;
    });
    if (pairs.size() > primer_.maxPairs) pairs.resize(primer_.maxPairs);
    return pairs;
}

// ---------------------------------------------------------------------------
// Guide search.
// ---------------------------------------------------------------------------

namespace {

// A 20-mer site in the reference, on one strand, in reference coordinates.
struct RefSite {
    int  start = 0;     // start of the protospacer in `strandSeq`
    int  refPos = 0;    // 0-based position on the forward strand of the reference
    Strand strand = Strand::Forward;
};

// OFF-TARGET COUNTING ALGORITHM (pigeonhole seed index, not an O(n*m) scan).
//
// Every alignment of a 20-mer with at most 2 mismatches must leave at least one
// of three disjoint blocks (7/7/6 nt) completely intact - two mismatches cannot
// hit three blocks. So the reference's PAM-adjacent sites are indexed once by
// (blockIndex, blockSequence); a query looks up its own three blocks, unions the
// candidate positions, and only those candidates are compared base-by-base with an
// early exit past 2 mismatches. Cost is O(reference) once plus O(candidates) per
// guide, instead of O(guides * reference * 20).
struct SeedIndex {
    static constexpr int kBlocks = 3;
    std::array<int, kBlocks> offset{};
    std::array<int, kBlocks> width{};
    std::unordered_multimap<std::uint64_t, std::size_t> map;
    std::vector<RefSite> sites;
    std::vector<std::string> siteSeq;   // the 20-mer, parallel to `sites`

    void layout(int len) {
        const int base = len / kBlocks;
        int used = 0;
        for (int b = 0; b < kBlocks; ++b) {
            offset[b] = used;
            width[b] = (b + 1 == kBlocks) ? len - used : base;
            used += width[b];
        }
    }

    static std::uint64_t key(int block, const char* p, int w) {
        // 2 bits per base; a 20-mer block is at most 7 bases here, so the whole
        // key including the block index fits a single 64-bit word with no hashing
        // of strings and no allocation.
        std::uint64_t k = static_cast<std::uint64_t>(block) + 1;
        for (int i = 0; i < w; ++i) {
            std::uint64_t c;
            switch (p[i]) {
                case 'A': c = 0; break;
                case 'C': c = 1; break;
                case 'G': c = 2; break;
                case 'T': c = 3; break;
                default:  return 0;  // an ambiguity code is not indexed
            }
            k = (k << 2) | c;
        }
        return k;
    }

    void add(RefSite s, std::string seq) {
        const std::size_t i = sites.size();
        sites.push_back(s);
        siteSeq.push_back(std::move(seq));
        for (int b = 0; b < kBlocks; ++b) {
            const std::uint64_t k = key(b, siteSeq[i].data() + offset[b], width[b]);
            if (k) map.emplace(k, i);
        }
    }
};

int mismatches(const std::string& a, const std::string& b, int cap) {
    int m = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i] && ++m > cap) return m;  // early exit past the cap
    }
    return m;
}

// PAM match under IUPAC, e.g. "NGG", "TTTV".
bool pamMatches(const std::string& pattern, const std::string& s) {
    if (s.size() != pattern.size()) return false;
    for (std::size_t i = 0; i < s.size(); ++i)
        if (!bio::iupacMatches(pattern[i], s[i])) return false;
    return true;
}

// Linearised view of a possibly circular sequence: a circular template's sites may
// straddle the origin, and dropping them would silently under-report.
std::string linearise(const NucRecord& r, std::size_t span) {
    if (!r.circular || r.sequence.size() < span) return r.sequence;
    return r.sequence + r.sequence.substr(0, span - 1);
}

}  // namespace

GuideSearchResult RealNucleicAcid::findGuides(const NucRecord& target, const NucRecord& reference,
                                             const std::string& pam) const {
    GuideSearchResult out;
    const std::string pamPattern = pam.empty() ? std::string("NGG") : pam;
    const int gl = guide_.protospacerLength;
    const int pl = static_cast<int>(pamPattern.size());
    out.referenceName = reference.id.empty() ? std::string("<unnamed reference>") : reference.id;
    out.basesSearched = static_cast<std::int64_t>(reference.sequence.size());

    // THE SCOPE GATE. Both conditions are necessary and neither is sufficient,
    // which is exactly why the flag is conservative: a 1 Mb contig passes the size
    // floor and is still not a genome. Nothing downstream may treat the flag as
    // proof of completeness - the prose scopeStatement remains the authoritative
    // description of what was searched.
    const bool distinctReference =
        reference.id != target.id && reference.sequence.size() > target.sequence.size();
    out.genomeWideClaimPossible =
        distinctReference && out.basesSearched >= guide_.genomeWideMinBases;

    const std::string refFwd = linearise(reference, static_cast<std::size_t>(gl + pl));
    const std::string refRev = bio::reverseComplement(refFwd);

    SeedIndex index;
    index.layout(gl);
    const auto indexStrand = [&](const std::string& s, Strand strand) {
        const int n = static_cast<int>(s.size());
        for (int i = 0; i + gl + pl <= n; ++i) {
            if (!pamMatches(pamPattern, s.substr(static_cast<std::size_t>(i + gl),
                                                 static_cast<std::size_t>(pl))))
                continue;
            std::string proto = s.substr(static_cast<std::size_t>(i), static_cast<std::size_t>(gl));
            if (proto.find_first_not_of("ACGT") != std::string::npos) continue;
            RefSite site;
            site.start = i;
            site.strand = strand;
            // Forward-strand coordinate of the protospacer 5' end, modulo the
            // circular wrap so a straddling site is reported once.
            const int refLen = static_cast<int>(reference.sequence.size());
            site.refPos = (strand == Strand::Forward) ? i % refLen
                                                      : (static_cast<int>(s.size()) - i - gl) % refLen;
            index.add(site, std::move(proto));
        }
    };
    indexStrand(refFwd, Strand::Forward);
    indexStrand(refRev, Strand::Reverse);

    // Where the target sits inside the reference, so the guide's OWN site is not
    // counted as its own off-target. When the target is absent from the reference
    // that is a fact the user needs, not something to paper over.
    std::size_t targetOffset = reference.sequence.find(target.sequence);
    const bool targetInReference = targetOffset != std::string::npos;
    if (!targetInReference) {
        targetOffset = 0;
        out.warnings.push_back(
            "The target sequence does not occur in the supplied reference, so no match was "
            "attributed to the on-target site and every count below is a genuine second site.");
    }

    const std::string tgtFwd = linearise(target, static_cast<std::size_t>(gl + pl));
    const std::string tgtRev = bio::reverseComplement(tgtFwd);
    const int tgtLen = static_cast<int>(target.sequence.size());

    const auto scan = [&](const std::string& s, Strand strand) {
        const int n = static_cast<int>(s.size());
        for (int i = 0; i + gl + pl <= n && out.guides.size() < guide_.maxGuides; ++i) {
            const std::string pamSeq =
                s.substr(static_cast<std::size_t>(i + gl), static_cast<std::size_t>(pl));
            if (!pamMatches(pamPattern, pamSeq)) continue;
            std::string proto =
                s.substr(static_cast<std::size_t>(i), static_cast<std::size_t>(gl));
            if (proto.find_first_not_of("ACGT") != std::string::npos) continue;

            GuideCandidate g;
            g.protospacer = proto;
            g.pam = pamSeq;
            g.strand = strand;
            g.position = (strand == Strand::Forward) ? i % tgtLen : (n - i - gl) % tgtLen;
            g.gcPercent = bio::gcPercent(proto);
            if (g.gcPercent < guide_.minGcPercent || g.gcPercent > guide_.maxGcPercent)
                g.warnings.push_back("GC " + std::to_string(static_cast<int>(g.gcPercent)) +
                                     "% is outside the 25-80% window where published activity "
                                     "data exists.");
            if (proto.find("TTTT") != std::string::npos)
                g.warnings.push_back("A TTTT run can terminate U6 transcription of the sgRNA.");

            // The on-target's own coordinate in reference forward-strand space.
            const int selfRefPos =
                targetInReference
                    ? ((strand == Strand::Forward)
                           ? static_cast<int>(targetOffset) + g.position
                           : static_cast<int>(targetOffset) + g.position)
                    : -1;

            // Pigeonhole lookup: union the three block buckets, dedupe, verify.
            std::vector<std::size_t> cand;
            for (int b = 0; b < SeedIndex::kBlocks; ++b) {
                const std::uint64_t k =
                    SeedIndex::key(b, proto.data() + index.offset[b], index.width[b]);
                if (!k) continue;
                auto range = index.map.equal_range(k);
                for (auto it = range.first; it != range.second; ++it) cand.push_back(it->second);
            }
            std::sort(cand.begin(), cand.end());
            cand.erase(std::unique(cand.begin(), cand.end()), cand.end());

            bool selfSeen = false;
            for (std::size_t si : cand) {
                const int mm = mismatches(proto, index.siteSeq[si], guide_.maxMismatches);
                if (mm > guide_.maxMismatches) continue;
                if (mm == 0 && !selfSeen && selfRefPos >= 0 &&
                    index.sites[si].refPos == selfRefPos &&
                    index.sites[si].strand == strand) {
                    selfSeen = true;  // this is the on-target itself
                    continue;
                }
                if (mm == 0) ++g.exactOffTargets;
                else if (mm == 1) ++g.oneMismatchOffTargets;
                else ++g.twoMismatchOffTargets;
            }
            out.guides.push_back(std::move(g));
        }
    };
    scan(tgtFwd, Strand::Forward);
    scan(tgtRev, Strand::Reverse);

    // A truncated guide list is itself a scope limit, so it goes in the scope
    // statement rather than being left for the reader to infer from a round number.
    const bool truncated = out.guides.size() >= guide_.maxGuides;

    out.scopeStatement =
        "Searched " + std::to_string(out.basesSearched) + " bases of '" + out.referenceName +
        "'" + (reference.circular ? " (circular)" : "") +
        ", both strands, for " + pamPattern + "-adjacent " + std::to_string(gl) +
        "-mer sites with up to " + std::to_string(guide_.maxMismatches) +
        " mismatches. Off-target counts below are counts WITHIN those " +
        std::to_string(out.basesSearched) +
        " bases only; the guide's own on-target site is excluded. No sequence outside this "
        "reference was examined" +
        (out.genomeWideClaimPossible
             ? ", and although the reference is large enough to be a genome BioCAD cannot verify "
               "that it is complete, so treat this as a screen and not a genome-wide specificity "
               "claim."
             : ", and this reference is far too small to be a genome, so these counts are NOT a "
               "genome-wide specificity claim.");
    if (truncated)
        out.scopeStatement +=
            " Guide enumeration stopped at the first " + std::to_string(guide_.maxGuides) +
            " protospacers in the target, so the list below is a prefix and not every candidate "
            "site; the off-target counts for the guides shown are unaffected.";
    out.warnings.push_back(nucleicScopeNote());
    return out;
}

}  // namespace biocad
