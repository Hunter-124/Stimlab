#include "bio/NucSeq.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <mutex>

#include "core/Assets.h"
#include "core/Error.h"

namespace biocad::bio {

namespace {

// The IUPAC pairing, written once. Anything not listed complements to N.
struct Pair { char from, to; };
constexpr Pair kComplement[] = {
    {'A', 'T'}, {'T', 'A'}, {'U', 'A'}, {'G', 'C'}, {'C', 'G'},
    {'R', 'Y'}, {'Y', 'R'},            // A/G  <-> C/T
    {'K', 'M'}, {'M', 'K'},            // G/T  <-> A/C
    {'B', 'V'}, {'V', 'B'},            // C/G/T <-> A/C/G
    {'D', 'H'}, {'H', 'D'},            // A/G/T <-> A/C/T
    {'S', 'S'}, {'W', 'W'}, {'N', 'N'},
};

// What each IUPAC symbol may actually be. Used both by iupacMatches() and by
// codon translation, so an ambiguous codon whose expansions all agree translates
// to that amino acid instead of X.
const char* expansionOf(char c) {
    switch (c) {
        case 'A': return "A";
        case 'C': return "C";
        case 'G': return "G";
        case 'T': return "T";
        case 'U': return "T";
        case 'R': return "AG";
        case 'Y': return "CT";
        case 'S': return "CG";
        case 'W': return "AT";
        case 'K': return "GT";
        case 'M': return "AC";
        case 'B': return "CGT";
        case 'D': return "AGT";
        case 'H': return "ACT";
        case 'V': return "ACG";
        case 'N': return "ACGT";
        default:  return "";
    }
}

char upper(char c) { return static_cast<char>(std::toupper(static_cast<unsigned char>(c))); }

int baseIndex(char c) {
    switch (c) {
        case 'T': case 'U': return 0;   // NCBI codon order is T, C, A, G
        case 'C': return 1;
        case 'A': return 2;
        case 'G': return 3;
        default:  return -1;
    }
}

}  // namespace

bool isIupac(char c) { return *expansionOf(upper(c)) != '\0'; }

char complementBase(char c) {
    const char u = upper(c);
    for (const auto& p : kComplement) {
        if (p.from == u) return p.to;
    }
    return 'N';
}

std::string reverseComplement(std::string_view seq) {
    std::string out;
    out.reserve(seq.size());
    for (std::size_t i = seq.size(); i-- > 0;) out.push_back(complementBase(seq[i]));
    return out;
}

std::string reverseComplementRna(std::string_view seq) {
    std::string out = reverseComplement(seq);
    for (char& c : out) {
        if (c == 'T') c = 'U';
    }
    return out;
}

std::string sanitizeSequence(std::string_view raw, std::vector<std::string>* warnings) {
    std::string out;
    out.reserve(raw.size());
    int dropped = 0;
    for (char c : raw) {
        const unsigned char u = static_cast<unsigned char>(c);
        if (std::isspace(u) || std::isdigit(u)) continue;
        const char up = upper(c);
        if (isIupac(up)) {
            out.push_back(up);
        } else {
            ++dropped;
        }
    }
    if (dropped > 0 && warnings) {
        warnings->push_back("dropped " + std::to_string(dropped) +
                            " character(s) that are not IUPAC nucleotide symbols");
    }
    return out;
}

double gcPercent(std::string_view seq) {
    // Only symbols that are unambiguous about GC content participate. Counting N
    // as half a GC, which some tools do, invents information.
    int gc = 0, at = 0;
    for (char c : seq) {
        switch (upper(c)) {
            case 'G': case 'C': case 'S': ++gc; break;
            case 'A': case 'T': case 'U': case 'W': ++at; break;
            default: break;
        }
    }
    const int total = gc + at;
    return total == 0 ? 0.0 : 100.0 * gc / total;
}

double gc3Percent(std::string_view codingSeq) {
    std::string third;
    third.reserve(codingSeq.size() / 3 + 1);
    for (std::size_t i = 2; i < codingSeq.size(); i += 3) third.push_back(codingSeq[i]);
    return gcPercent(third);
}

bool iupacMatches(char pattern, char base) {
    const char* pat = expansionOf(upper(pattern));
    const char* obs = expansionOf(upper(base));
    if (*pat == '\0' || *obs == '\0') return false;
    // The base must be *entirely* consistent with the pattern: pattern N accepts
    // base A, but pattern A must not accept base N, or a degenerate template
    // would report sites that may not be there.
    for (const char* o = obs; *o; ++o) {
        bool found = false;
        for (const char* p = pat; *p; ++p) {
            if (*p == *o) { found = true; break; }
        }
        if (!found) return false;
    }
    return true;
}

std::vector<int> findPattern(std::string_view seq, std::string_view pattern, bool circular) {
    std::vector<int> hits;
    const int n = static_cast<int>(seq.size());
    const int m = static_cast<int>(pattern.size());
    if (m == 0 || n == 0) return hits;
    // On a linear template the last possible start is n - m. On a circular one
    // every position is a possible start and the comparison wraps; the modulo is
    // the whole trick behind a site spanning the origin being reported once.
    const int last = circular ? n - 1 : n - m;
    if (!circular && last < 0) return hits;
    for (int start = 0; start <= last; ++start) {
        bool ok = true;
        for (int k = 0; k < m; ++k) {
            const int idx = circular ? (start + k) % n : start + k;
            if (!iupacMatches(pattern[static_cast<std::size_t>(k)],
                              seq[static_cast<std::size_t>(idx)])) {
                ok = false;
                break;
            }
        }
        if (ok) hits.push_back(start);
    }
    return hits;
}

// ------------------------------------------------------------- genetic codes

GeneticCode::GeneticCode(int id, std::string name, std::string aminoAcids,
                         std::vector<std::string> starts, std::string notes)
    : id_(id), name_(std::move(name)), aminoAcids_(std::move(aminoAcids)),
      starts_(std::move(starts)), notes_(std::move(notes)) {
    if (aminoAcids_.size() != 64) {
        throw Error::parse("genetic code " + std::to_string(id_) + ": aminoAcids must be 64 "
                           "characters, got " + std::to_string(aminoAcids_.size()));
    }
}

char GeneticCode::translateCodon(std::string_view codon) const {
    if (codon.size() != 3) return 'X';
    const char* e0 = expansionOf(upper(codon[0]));
    const char* e1 = expansionOf(upper(codon[1]));
    const char* e2 = expansionOf(upper(codon[2]));
    if (*e0 == '\0' || *e1 == '\0' || *e2 == '\0') return 'X';
    // Expand the ambiguity. GCN is unambiguously Ala and CTN unambiguously Leu,
    // and translating those to X would lose real information; a codon whose
    // expansions disagree is genuinely unknown.
    char seen = 0;
    for (const char* a = e0; *a; ++a) {
        for (const char* b = e1; *b; ++b) {
            for (const char* c = e2; *c; ++c) {
                const int idx = baseIndex(*a) * 16 + baseIndex(*b) * 4 + baseIndex(*c);
                const char aa = aminoAcids_[static_cast<std::size_t>(idx)];
                if (seen == 0) seen = aa;
                else if (seen != aa) return 'X';
            }
        }
    }
    return seen == 0 ? 'X' : seen;
}

bool GeneticCode::isStop(std::string_view codon) const { return translateCodon(codon) == '*'; }

bool GeneticCode::isStart(std::string_view codon) const {
    if (codon.size() != 3) return false;
    std::string up;
    for (char c : codon) up.push_back(upper(c) == 'U' ? 'T' : upper(c));
    return std::find(starts_.begin(), starts_.end(), up) != starts_.end();
}

const GeneticCode* GeneticCodeTables::find(int id) const {
    for (const auto& c : codes_) {
        if (c.id() == id) return &c;
    }
    return nullptr;
}

std::vector<int> GeneticCodeTables::ids() const {
    std::vector<int> out;
    out.reserve(codes_.size());
    for (const auto& c : codes_) out.push_back(c.id());
    return out;
}

GeneticCodeTables parseGeneticCodes(const nlohmann::json& j) {
    const int version = j.value("schemaVersion", 0);
    if (version != 1) {
        throw Error::parse("genetic-codes pack: unsupported schemaVersion " +
                           std::to_string(version));
    }
    // The codon order is part of the pack rather than assumed, and it is checked
    // rather than trusted: a table string interpreted in the wrong order is a
    // silently wrong protein.
    const auto order = j.value("codonOrder", std::vector<std::string>{});
    if (order.size() != 64) throw Error::parse("genetic-codes pack: codonOrder must list 64 codons");
    for (int i = 0; i < 64; ++i) {
        const auto& c = order[static_cast<std::size_t>(i)];
        if (c.size() != 3 || baseIndex(c[0]) * 16 + baseIndex(c[1]) * 4 + baseIndex(c[2]) != i) {
            throw Error::parse("genetic-codes pack: codonOrder is not the NCBI order at index " +
                               std::to_string(i));
        }
    }
    GeneticCodeTables out;
    const auto codes = j.find("codes");
    if (codes == j.end() || !codes->is_array() || codes->empty()) {
        throw Error::parse("genetic-codes pack: no codes");
    }
    for (const auto& c : *codes) {
        out.add(GeneticCode(c.value("id", 0), c.value("name", std::string{}),
                            c.value("aminoAcids", std::string{}),
                            c.value("starts", std::vector<std::string>{}),
                            c.value("notes", std::string{})));
    }
    return out;
}

GeneticCodeTables loadGeneticCodes(const std::filesystem::path& file) {
    std::ifstream in(file, std::ios::binary);
    if (!in) throw Error::io("cannot open genetic-codes pack: " + file.string());
    nlohmann::json j;
    try {
        in >> j;
    } catch (const nlohmann::json::exception& e) {
        throw Error::parse(file.string() + ": " + e.what());
    }
    try {
        return parseGeneticCodes(j);
    } catch (const Error& e) {
        throw Error::parse(file.string() + ": " + e.message);
    }
}

const GeneticCodeTables& builtinGeneticCodes() {
    static GeneticCodeTables tables;
    static std::once_flag once;
    static std::string failure;
    std::call_once(once, [] {
        const auto dir = core::assetDir("packs/nucleic");
        if (dir.empty()) {
            failure = "asset root not found";
            return;
        }
        try {
            tables = loadGeneticCodes(dir / "genetic-codes.json");
        } catch (const Error& e) {
            failure = e.message;
        }
    });
    if (tables.empty()) throw Error::io("genetic codes unavailable: " + failure);
    return tables;
}

const GeneticCode& builtinGeneticCode(int id) {
    const auto* code = builtinGeneticCodes().find(id);
    if (!code) {
        throw Error::notFound("no NCBI translation table " + std::to_string(id) + " in the pack");
    }
    return *code;
}

// --------------------------------------------------------------- translation

std::string translate(std::string_view dna, const GeneticCode& code,
                      std::vector<std::string>* warnings) {
    std::string protein;
    protein.reserve(dna.size() / 3);
    std::size_t i = 0;
    for (; i + 3 <= dna.size(); i += 3) protein.push_back(code.translateCodon(dna.substr(i, 3)));
    if (i != dna.size() && warnings) {
        warnings->push_back("trailing " + std::to_string(dna.size() - i) +
                            " base(s) are not a whole codon and were not translated");
    }
    return protein;
}

std::vector<std::string> sixFrameTranslation(std::string_view seq, const GeneticCode& code) {
    const std::string rc = reverseComplement(seq);
    std::vector<std::string> frames;
    frames.reserve(6);
    for (int f = 0; f < 3; ++f) {
        frames.push_back(static_cast<std::size_t>(f) < seq.size()
                             ? translate(seq.substr(static_cast<std::size_t>(f)), code)
                             : std::string{});
    }
    for (int f = 0; f < 3; ++f) {
        frames.push_back(static_cast<std::size_t>(f) < rc.size()
                             ? translate(std::string_view(rc).substr(static_cast<std::size_t>(f)),
                                         code)
                             : std::string{});
    }
    return frames;
}

namespace {

// ORF scan on one strand of `work`, which is the strand's 5'->3' sequence. For a
// circular template `work` is the sequence doubled, and `n` is the real length:
// an ORF may then start before n and end after it (a wrap) but may never be
// longer than one turn, which is what stops a doubled template reporting an
// ORF twice.
void scanStrand(const std::string& work, int n, bool circular, const GeneticCode& code,
                const OrfOptions& options, Strand strand, int totalLength,
                std::vector<OpenReadingFrame>& out) {
    const int len = static_cast<int>(work.size());
    for (int frame = 0; frame < 3; ++frame) {
        for (int start = frame; start + 3 <= len; start += 3) {
            if (circular && start >= n) break;   // every wrap start is already covered
            const std::string_view startCodon(work.data() + start, 3);
            const bool isStart = options.alternativeStarts ? code.isStart(startCodon)
                                                           : startCodon == "ATG";
            if (!isStart) continue;
            std::string protein;
            int pos = start;
            bool stopped = false;
            for (; pos + 3 <= len; pos += 3) {
                if (circular && pos - start >= n) break;   // one turn only
                const std::string_view codon(work.data() + pos, 3);
                if (code.isStop(codon)) {
                    stopped = true;
                    pos += 3;
                    break;
                }
                protein.push_back(code.translateCodon(codon));
            }
            if (options.requireStop && !stopped) continue;
            if (static_cast<int>(protein.size()) < options.minAminoAcids) continue;
            OpenReadingFrame orf;
            orf.strand = strand;
            orf.frame = frame;
            orf.protein = protein;
            orf.stopped = stopped;
            if (strand == Strand::Forward) {
                orf.begin = start;
                orf.end = pos;
            } else {
                // Reverse-strand coordinates are converted back to top-strand
                // indices so every coordinate in the DTO means the same thing.
                orf.begin = totalLength - pos;
                orf.end = totalLength - start;
                if (orf.begin < 0) {              // a wrap on the reverse strand
                    orf.begin += totalLength;
                    orf.end += totalLength;
                }
            }
            out.push_back(std::move(orf));
        }
    }
}

}  // namespace

std::vector<OpenReadingFrame> findOrfs(std::string_view seq, const GeneticCode& code,
                                       const OrfOptions& options) {
    std::vector<OpenReadingFrame> out;
    const int n = static_cast<int>(seq.size());
    if (n < 3) return out;
    std::string fwd(seq);
    std::string rev = reverseComplement(seq);
    if (options.circular) {
        // Doubling is the whole circular implementation: a wrapping ORF is a
        // plain linear ORF in the doubled string, capped at one turn above.
        fwd += fwd;
        rev += rev;
    }
    scanStrand(fwd, n, options.circular, code, options, Strand::Forward, n, out);
    scanStrand(rev, n, options.circular, code, options, Strand::Reverse, n, out);
    std::sort(out.begin(), out.end(), [](const OpenReadingFrame& a, const OpenReadingFrame& b) {
        if (a.protein.size() != b.protein.size()) return a.protein.size() > b.protein.size();
        if (a.begin != b.begin) return a.begin < b.begin;
        return a.strand < b.strand;
    });
    return out;
}

TranslationResult translateRecord(const NucRecord& record, int geneticCodeId,
                                  const OrfOptions& options) {
    TranslationResult out;
    out.recordId = record.id;
    out.geneticCodeId = geneticCodeId;
    const GeneticCode& code = builtinGeneticCode(geneticCodeId);
    out.frames = sixFrameTranslation(record.sequence, code);
    OrfOptions opts = options;
    opts.circular = record.circular;
    out.orfs = findOrfs(record.sequence, code, opts);
    if (record.circular) {
        out.warnings.push_back("circular template: an ORF may wrap the origin, so its end "
                               "coordinate can exceed the sequence length");
    }
    if (!opts.alternativeStarts && (geneticCodeId == 2 || geneticCodeId == 3 ||
                                    geneticCodeId == 4 || geneticCodeId == 5)) {
        out.warnings.push_back("translation table " + std::to_string(geneticCodeId) +
                               " documents initiation codons other than ATG; only ATG starts "
                               "were scanned because alternativeStarts was not requested");
    }
    return out;
}

std::string featureSequence(const NucRecord& record, const NucFeature& feature) {
    const int n = static_cast<int>(record.sequence.size());
    std::string joined;
    for (const auto& part : feature.parts) {
        for (int i = part.first; i < part.second; ++i) {
            int idx = i;
            if (record.circular && n > 0) {
                idx = ((i % n) + n) % n;
            } else if (idx < 0 || idx >= n) {
                continue;   // a partial location beyond a linear end contributes nothing
            }
            joined.push_back(record.sequence[static_cast<std::size_t>(idx)]);
        }
    }
    // GenBank complement(join(a,b)) is the reverse complement of the joined
    // forward sequence, so the complement is applied after concatenation, never
    // per part.
    return feature.strand == Strand::Reverse ? reverseComplement(joined) : joined;
}

}  // namespace biocad::bio
