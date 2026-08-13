#pragma once

// Nucleic-acid sequence primitives: the IUPAC alphabet, reverse complement,
// NCBI genetic codes, translation, six-frame reading and ORF finding.
//
// Two decisions drive this file.
//
// First, the alphabet is the full IUPAC ambiguity set, not ACGT. A vector map, a
// primer binding site and a degenerate restriction site all contain ambiguity
// codes, and silently coercing them to N loses information that the restriction
// mapper needs. Reverse complement is therefore defined for every symbol and is
// an involution over arbitrary IUPAC strings - which is a test, not a comment.
//
// Second, the genetic code is data, not code. There are 33 NCBI translation
// tables and they differ in both amino-acid assignment and initiation codons;
// hard-coding table 1 and pretending mitochondria do not exist is the kind of
// thing that silently produces a wrong protein. The tables come from
// assets/packs/nucleic/genetic-codes.json.

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "data/Nucleic.h"

namespace biocad::bio {

// --------------------------------------------------------------- IUPAC alphabet

// True for the 4 unambiguous bases, the 11 ambiguity codes and U (RNA uracil).
bool isIupac(char c);

// Complement of a single IUPAC symbol, preserving case-insensitively-uppercased
// output. R<->Y, K<->M, B<->V, D<->H; S, W and N are self-complementary. Returns
// 'N' for a symbol that is not IUPAC, because a gap in an alignment-derived
// sequence must not silently become a base.
char complementBase(char c);

// Reverse complement. U is complemented to A and A to T: an RNA input therefore
// comes back as DNA, which is the honest answer for a function that cannot know
// which strand chemistry the caller wants. Use reverseComplementRna for RNA.
std::string reverseComplement(std::string_view seq);
std::string reverseComplementRna(std::string_view seq);

// Uppercases, strips whitespace and digits (so a pasted GenBank ORIGIN block or a
// numbered listing can be fed in directly), and reports anything it dropped.
std::string sanitizeSequence(std::string_view raw, std::vector<std::string>* warnings = nullptr);

// GC as a percentage of the *unambiguous* bases; S counts as GC and W as AT
// because those codes are unambiguous about GC content. Other ambiguity codes are
// excluded from both numerator and denominator rather than guessed.
double gcPercent(std::string_view seq);
double gc3Percent(std::string_view codingSeq);   // third codon positions only

// True when `base` is consistent with the IUPAC `pattern` symbol. The asymmetry
// matters: pattern N matches base A, but pattern A does not match base N, since a
// site is only present if the sequence is known to contain it.
bool iupacMatches(char pattern, char base);

// Every 0-based start index at which `pattern` occurs. On a circular template a
// match may wrap the origin, so indices run over the whole length and a match is
// tested against the doubled sequence - which is what makes a site spanning the
// origin found exactly once rather than twice or not at all.
std::vector<int> findPattern(std::string_view seq, std::string_view pattern, bool circular);

// ------------------------------------------------------------- genetic codes

class GeneticCode {
public:
    GeneticCode() = default;
    GeneticCode(int id, std::string name, std::string aminoAcids, std::vector<std::string> starts,
                std::string notes);

    int id() const { return id_; }
    const std::string& name() const { return name_; }
    const std::string& notes() const { return notes_; }
    const std::vector<std::string>& startCodons() const { return starts_; }

    // 'X' when the codon contains an ambiguity code that spans more than one
    // amino acid, the single amino acid when every expansion agrees (so GCN is A
    // and CTN is L), '*' for a stop, and 'X' for a malformed codon.
    char translateCodon(std::string_view codon) const;
    bool isStop(std::string_view codon) const;
    bool isStart(std::string_view codon) const;   // documented initiators only

private:
    int                      id_ = 0;
    std::string              name_;
    std::string              aminoAcids_;   // 64 chars, NCBI codon order
    std::vector<std::string> starts_;
    std::string              notes_;
};

class GeneticCodeTables {
public:
    const GeneticCode* find(int id) const;
    std::vector<int>   ids() const;
    void               add(GeneticCode code) { codes_.push_back(std::move(code)); }
    bool               empty() const { return codes_.empty(); }

private:
    std::vector<GeneticCode> codes_;
};

GeneticCodeTables parseGeneticCodes(const nlohmann::json& j);
GeneticCodeTables loadGeneticCodes(const std::filesystem::path& file);

// Cached, resolved through core::assetDir("packs/nucleic"). Throws Error::io when
// the pack is missing: a translation with a silently substituted table is worse
// than a hard failure.
const GeneticCodeTables& builtinGeneticCodes();
const GeneticCode&       builtinGeneticCode(int id);

// --------------------------------------------------------------- translation

// Translates in frame from the first base. A trailing partial codon is ignored
// (and reported through `warnings` when supplied). Stops render as '*'.
std::string translate(std::string_view dna, const GeneticCode& code,
                      std::vector<std::string>* warnings = nullptr);

struct OrfOptions {
    int  minAminoAcids = 30;      // excluding the stop
    bool requireStop = true;      // false keeps ORFs running off a linear end
    bool circular = false;        // an ORF may wrap the origin
    bool alternativeStarts = false;  // any documented initiator, not just ATG
};

// Six-frame translation: frames[0..2] are the forward frames 0,1,2 and
// frames[3..5] the reverse-strand frames 0,1,2 read from the reverse
// complement's 5' end.
std::vector<std::string> sixFrameTranslation(std::string_view seq, const GeneticCode& code);

// ORFs, longest first. Coordinates are always top-strand 0-based half-open
// [begin, end) with `end` past the stop codon; for a reverse-strand ORF that
// means begin < end still holds and the protein reads from `end` backwards. On a
// circular template `end` may exceed the sequence length, which is how a wrapping
// ORF is expressed without duplicating the sequence.
std::vector<OpenReadingFrame> findOrfs(std::string_view seq, const GeneticCode& code,
                                       const OrfOptions& options);

// Six frames plus ORFs for a record, honouring record.circular.
TranslationResult translateRecord(const NucRecord& record, int geneticCodeId,
                                  const OrfOptions& options);

// The spliced, strand-oriented sequence of a feature: parts are concatenated in
// the order stored, and a reverse-strand feature is reverse-complemented after
// concatenation. Parts that wrap the origin of a circular record are honoured.
std::string featureSequence(const NucRecord& record, const NucFeature& feature);

}  // namespace biocad::bio
