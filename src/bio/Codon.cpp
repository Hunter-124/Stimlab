#include "bio/Codon.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <mutex>
#include <queue>
#include <unordered_set>

#include "core/Assets.h"
#include "core/Error.h"

namespace biocad::bio {

namespace {

constexpr char kBases[4] = {'A', 'C', 'G', 'T'};

int baseSlot(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return -1;
    }
}

// Homopolymer limits are expressed as ordinary forbidden patterns so the automaton
// below is the single place any constraint is enforced. One pattern per base of
// length maxHomopolymer + 1.
std::vector<std::string> expandPatterns(const std::vector<std::string>& patterns,
                                        int maxHomopolymer) {
    std::vector<std::string> out;
    const auto push = [&out](std::string p) {
        if (p.empty()) return;
        for (char& c : p) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        for (char c : p) {
            if (!isIupac(c)) {
                throw Error::invalidArgument("forbidden pattern '" + p + "' contains '" +
                                             std::string(1, c) + "', which is not an IUPAC symbol");
            }
        }
        if (std::find(out.begin(), out.end(), p) == out.end()) out.push_back(std::move(p));
    };
    for (const auto& p : patterns) {
        push(p);
        // A recognition site is a site on either strand, so the reverse complement
        // is forbidden too. For a palindromic site this is a no-op.
        push(reverseComplement(p));
    }
    if (maxHomopolymer > 0) {
        for (char b : kBases) push(std::string(static_cast<std::size_t>(maxHomopolymer) + 1, b));
    }
    return out;
}

// Aho-Corasick over the four DNA bases, with IUPAC pattern symbols expanded at
// edge-construction time: a pattern character matches every base consistent with
// it, so a degenerate site costs states but no special cases downstream.
class PatternAutomaton {
public:
    explicit PatternAutomaton(const std::vector<std::string>& patterns) {
        nodes_.push_back(Node{});
        for (std::size_t p = 0; p < patterns.size(); ++p) {
            // A pattern with an ambiguity code becomes several trie paths; the
            // number of paths is the product of its symbol expansions, which is
            // small for real recognition sites.
            insertExpansions(patterns[p], 0, std::string{}, p);
        }
        build();
    }

    int next(int state, char base) const {
        const int slot = baseSlot(base);
        if (slot < 0) return 0;   // an ambiguity code in the input cannot match a site
        return nodes_[static_cast<std::size_t>(state)].go[slot];
    }
    bool terminal(int state) const { return nodes_[static_cast<std::size_t>(state)].terminal; }
    int  patternAt(int state) const { return nodes_[static_cast<std::size_t>(state)].pattern; }
    int  depthAt(int state) const { return nodes_[static_cast<std::size_t>(state)].depth; }
    std::size_t size() const { return nodes_.size(); }

private:
    struct Node {
        int  go[4] = {-1, -1, -1, -1};
        int  fail = 0;
        int  depth = 0;
        bool terminal = false;
        int  pattern = -1;
    };

    void insertExpansions(const std::string& pattern, std::size_t index, std::string acc,
                          std::size_t patternIndex) {
        if (index == pattern.size()) {
            int node = 0;
            for (char c : acc) {
                const int slot = baseSlot(c);
                if (nodes_[static_cast<std::size_t>(node)].go[slot] < 0) {
                    nodes_.push_back(Node{});
                    nodes_.back().depth = nodes_[static_cast<std::size_t>(node)].depth + 1;
                    nodes_[static_cast<std::size_t>(node)].go[slot] =
                        static_cast<int>(nodes_.size()) - 1;
                }
                node = nodes_[static_cast<std::size_t>(node)].go[slot];
            }
            nodes_[static_cast<std::size_t>(node)].terminal = true;
            nodes_[static_cast<std::size_t>(node)].pattern = static_cast<int>(patternIndex);
            return;
        }
        for (char b : kBases) {
            if (iupacMatches(pattern[index], b)) {
                insertExpansions(pattern, index + 1, acc + b, patternIndex);
            }
        }
    }

    void build() {
        std::queue<int> q;
        for (int s = 0; s < 4; ++s) {
            if (nodes_[0].go[s] < 0) {
                nodes_[0].go[s] = 0;
            } else {
                nodes_[static_cast<std::size_t>(nodes_[0].go[s])].fail = 0;
                q.push(nodes_[0].go[s]);
            }
        }
        while (!q.empty()) {
            const int v = q.front();
            q.pop();
            // Terminality propagates along failure links: a state that ends a
            // shorter pattern is terminal even when reached through a longer path.
            if (nodes_[static_cast<std::size_t>(nodes_[static_cast<std::size_t>(v)].fail)].terminal) {
                nodes_[static_cast<std::size_t>(v)].terminal = true;
                if (nodes_[static_cast<std::size_t>(v)].pattern < 0) {
                    nodes_[static_cast<std::size_t>(v)].pattern =
                        nodes_[static_cast<std::size_t>(nodes_[static_cast<std::size_t>(v)].fail)]
                            .pattern;
                }
            }
            for (int s = 0; s < 4; ++s) {
                const int u = nodes_[static_cast<std::size_t>(v)].go[s];
                if (u < 0) {
                    nodes_[static_cast<std::size_t>(v)].go[s] =
                        nodes_[static_cast<std::size_t>(nodes_[static_cast<std::size_t>(v)].fail)]
                            .go[s];
                } else {
                    nodes_[static_cast<std::size_t>(u)].fail =
                        nodes_[static_cast<std::size_t>(nodes_[static_cast<std::size_t>(v)].fail)]
                            .go[s];
                    q.push(u);
                }
            }
        }
    }

    std::vector<Node> nodes_;
};

}  // namespace

// --------------------------------------------------------------- usage tables

CodonUsageTable::CodonUsageTable(std::string id, std::string name, std::string source,
                                 int geneticCodeId, std::vector<CodonUsageEntry> entries)
    : id_(std::move(id)), name_(std::move(name)), source_(std::move(source)),
      geneticCodeId_(geneticCodeId), entries_(std::move(entries)) {
    if (entries_.size() != 64) {
        throw Error::parse("codon usage table " + id_ + ": expected 64 codons, got " +
                           std::to_string(entries_.size()));
    }
}

const CodonUsageEntry* CodonUsageTable::find(std::string_view codon) const {
    for (const auto& e : entries_) {
        if (e.codon == codon) return &e;
    }
    return nullptr;
}

std::string CodonUsageTable::mostFrequentCodon(char aminoAcid) const {
    const CodonUsageEntry* best = nullptr;
    for (const auto& e : entries_) {
        if (e.aminoAcid.empty() || e.aminoAcid[0] != aminoAcid) continue;
        if (!best || e.frequencyPerThousand > best->frequencyPerThousand) best = &e;
    }
    return best ? best->codon : std::string{};
}

const CodonUsageTable* CodonUsageSet::find(std::string_view id) const {
    for (const auto& t : tables_) {
        if (t.id() == id) return &t;
    }
    return nullptr;
}

std::vector<std::string> CodonUsageSet::ids() const {
    std::vector<std::string> out;
    out.reserve(tables_.size());
    for (const auto& t : tables_) out.push_back(t.id());
    return out;
}

CodonUsageSet parseCodonUsage(const nlohmann::json& j) {
    const int version = j.value("schemaVersion", 0);
    if (version != 1) {
        throw Error::parse("codon-usage pack: unsupported schemaVersion " +
                           std::to_string(version));
    }
    const auto tables = j.find("tables");
    if (tables == j.end() || !tables->is_array() || tables->empty()) {
        throw Error::parse("codon-usage pack: no tables");
    }
    CodonUsageSet out;
    for (const auto& t : *tables) {
        std::vector<CodonUsageEntry> entries;
        for (const auto& e : t.value("entries", nlohmann::json::array())) {
            CodonUsageEntry entry;
            entry.codon = e.value("codon", std::string{});
            entry.aminoAcid = e.value("aminoAcid", std::string{});
            entry.relativeAdaptiveness = e.value("relativeAdaptiveness", 0.0);
            entry.frequencyPerThousand = e.value("frequencyPerThousand", 0.0);
            if (entry.codon.size() != 3) {
                throw Error::parse("codon-usage pack: bad codon '" + entry.codon + "'");
            }
            if (entry.relativeAdaptiveness < 0.0 || entry.relativeAdaptiveness > 1.0) {
                throw Error::parse("codon-usage pack: " + entry.codon +
                                   " has relativeAdaptiveness outside [0, 1]");
            }
            entries.push_back(std::move(entry));
        }
        out.add(CodonUsageTable(t.value("id", std::string{}), t.value("name", std::string{}),
                                t.value("source", std::string{}), t.value("geneticCodeId", 1),
                                std::move(entries)));
    }
    return out;
}

CodonUsageSet loadCodonUsage(const std::filesystem::path& file) {
    std::ifstream in(file, std::ios::binary);
    if (!in) throw Error::io("cannot open codon-usage pack: " + file.string());
    nlohmann::json j;
    try {
        in >> j;
    } catch (const nlohmann::json::exception& e) {
        throw Error::parse(file.string() + ": " + e.what());
    }
    try {
        return parseCodonUsage(j);
    } catch (const Error& e) {
        throw Error::parse(file.string() + ": " + e.message);
    }
}

const CodonUsageSet& builtinCodonUsage() {
    static CodonUsageSet set;
    static std::once_flag once;
    static std::string failure;
    std::call_once(once, [] {
        const auto dir = core::assetDir("packs/nucleic");
        if (dir.empty()) {
            failure = "asset root not found";
            return;
        }
        try {
            set = loadCodonUsage(dir / "codon-usage.json");
        } catch (const Error& e) {
            failure = e.message;
        }
    });
    if (set.empty()) throw Error::io("codon usage tables unavailable: " + failure);
    return set;
}

const CodonUsageTable& builtinCodonUsageTable(std::string_view id) {
    const auto* t = builtinCodonUsage().find(id);
    if (!t) throw Error::notFound("no codon usage table '" + std::string(id) + "' in the pack");
    return *t;
}

// ---------------------------------------------------------------- CAI and GC

CodonMetrics codonMetrics(std::string_view cds, const CodonUsageTable& table) {
    CodonMetrics out;
    out.usageTableName = table.name();
    const std::string source = "arithmetic on " + table.source();
    if (cds.size() < 3) {
        out.cai = notComputed("at least one complete codon");
        out.gcPercent = notComputed("at least one base");
        out.gc3Percent = notComputed("at least one complete codon");
        return out;
    }
    if (cds.size() % 3 != 0) {
        out.warnings.push_back("the sequence is not a multiple of three; the trailing " +
                               std::to_string(cds.size() % 3) +
                               " base(s) were excluded from CAI and GC3");
    }
    double sumLog = 0.0;
    int counted = 0, skippedSingle = 0, skippedZero = 0, skippedUnknown = 0;
    for (std::size_t i = 0; i + 3 <= cds.size(); i += 3) {
        const std::string codon(cds.substr(i, 3));
        const CodonUsageEntry* e = table.find(codon);
        if (!e || e->aminoAcid.empty()) {
            ++skippedUnknown;
            continue;
        }
        const char aa = e->aminoAcid[0];
        if (aa == '*') continue;                 // stops are outside the convention
        if (aa == 'M' || aa == 'W') {            // single-codon families
            ++skippedSingle;
            continue;
        }
        if (e->relativeAdaptiveness <= 0.0) {
            ++skippedZero;
            continue;
        }
        sumLog += std::log(e->relativeAdaptiveness);
        ++counted;
    }
    if (counted == 0) {
        out.cai = notComputed("at least one codon in a multi-codon family with a non-zero w");
    } else {
        out.cai = makeQuantity(std::exp(sumLog / counted), "", 0.0, Provenance::Measured,
                               source + "; CAI is the geometric mean of w over " +
                                   std::to_string(counted) +
                                   " codons, excluding Met, Trp and stops");
    }
    out.gcPercent = makeQuantity(gcPercent(cds), "%", 0.0, Provenance::Measured,
                                 "base composition of the given sequence");
    out.gc3Percent = makeQuantity(gc3Percent(cds), "%", 0.0, Provenance::Measured,
                                  "third-codon-position composition of the given sequence");
    if (skippedSingle > 0) {
        out.warnings.push_back(std::to_string(skippedSingle) +
                               " Met/Trp codon(s) excluded from CAI: a single-codon family has "
                               "w = 1 by construction");
    }
    if (skippedZero > 0) {
        out.warnings.push_back(std::to_string(skippedZero) +
                               " codon(s) with w = 0 in this table excluded from CAI; including "
                               "them would drive the geometric mean to zero");
    }
    if (skippedUnknown > 0) {
        out.warnings.push_back(std::to_string(skippedUnknown) +
                               " codon(s) are not in the table (ambiguity codes?) and were "
                               "excluded");
    }
    return out;
}

std::vector<std::string> findForbidden(std::string_view seq,
                                       const std::vector<std::string>& patterns,
                                       int maxHomopolymer) {
    const auto expanded = expandPatterns(patterns, maxHomopolymer);
    std::vector<std::string> hits;
    if (expanded.empty()) return hits;
    const PatternAutomaton automaton(expanded);
    int state = 0;
    for (std::size_t i = 0; i < seq.size(); ++i) {
        state = automaton.next(state, seq[i]);
        if (automaton.terminal(state)) {
            const int p = automaton.patternAt(state);
            const int depth = automaton.depthAt(state);
            hits.push_back((p >= 0 ? expanded[static_cast<std::size_t>(p)] : std::string("?")) +
                           "@" + std::to_string(static_cast<int>(i) - depth + 1));
        }
    }
    return hits;
}

// --------------------------------------------------------------- optimization

namespace {

struct Choice { std::string codon; double logW = 0.0; int gc = 0; };

// The synonymous codons of one amino acid under the requested table, with the
// weight the DP maximises. A family whose codons all have w = 0 in the usage table
// would make every path -inf, so a floor is applied and reported.
std::vector<Choice> choicesFor(char aa, const CodonUsageTable& table, const GeneticCode& code,
                               bool& usedFloor) {
    std::vector<Choice> out;
    for (const auto& e : table.entries()) {
        if (code.translateCodon(e.codon) != aa) continue;
        Choice c;
        c.codon = e.codon;
        double w = e.relativeAdaptiveness;
        if (w <= 0.0) {
            w = 1e-6;   // representable, still ranked last
            usedFloor = true;
        }
        c.logW = std::log(w);
        c.gc = 0;
        for (char b : e.codon) {
            if (b == 'G' || b == 'C') ++c.gc;
        }
        out.push_back(std::move(c));
    }
    std::sort(out.begin(), out.end(),
              [](const Choice& a, const Choice& b) { return a.logW > b.logW; });
    return out;
}

struct DpResult {
    bool        feasible = false;
    std::string sequence;
    double      score = 0.0;    // sum of log w plus the GC bias term
    int         firstInfeasiblePosition = -1;
};

// One pass of the dynamic program. `gcBias` is the Lagrange multiplier on GC
// content: adding gcBias per G/C base tilts the choice without ever relaxing the
// pattern constraint, which stays hard.
DpResult runDp(const std::string& protein, const std::vector<std::vector<Choice>>& choices,
               const PatternAutomaton& automaton, double gcBias) {
    const std::size_t n = protein.size();
    const std::size_t states = automaton.size();
    constexpr double kNegInf = -std::numeric_limits<double>::infinity();
    std::vector<double> best(states, kNegInf), nextBest(states, kNegInf);
    // Back-pointers: for each position and state, which previous state and which
    // codon got here. Stored flat; n is a protein length, so this is small.
    std::vector<std::vector<int>> fromState(n), fromChoice(n);
    best[0] = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        fromState[i].assign(states, -1);
        fromChoice[i].assign(states, -1);
        std::fill(nextBest.begin(), nextBest.end(), kNegInf);
        for (std::size_t s = 0; s < states; ++s) {
            if (best[s] == kNegInf) continue;
            for (std::size_t c = 0; c < choices[i].size(); ++c) {
                const Choice& choice = choices[i][c];
                int state = static_cast<int>(s);
                bool blocked = false;
                for (char b : choice.codon) {
                    state = automaton.next(state, b);
                    if (automaton.terminal(state)) {
                        blocked = true;
                        break;
                    }
                }
                if (blocked) continue;
                const double score =
                    best[s] + choice.logW + gcBias * static_cast<double>(choice.gc);
                if (score > nextBest[static_cast<std::size_t>(state)]) {
                    nextBest[static_cast<std::size_t>(state)] = score;
                    fromState[i][static_cast<std::size_t>(state)] = static_cast<int>(s);
                    fromChoice[i][static_cast<std::size_t>(state)] = static_cast<int>(c);
                }
            }
        }
        best.swap(nextBest);
        if (std::none_of(best.begin(), best.end(), [](double v) { return v > kNegInf; })) {
            DpResult out;
            out.firstInfeasiblePosition = static_cast<int>(i);
            return out;
        }
    }
    DpResult out;
    std::size_t bestState = 0;
    double bestScore = kNegInf;
    for (std::size_t s = 0; s < states; ++s) {
        if (best[s] > bestScore) {
            bestScore = best[s];
            bestState = s;
        }
    }
    out.feasible = bestScore > kNegInf;
    if (!out.feasible) return out;
    out.score = bestScore;
    std::string reversed;
    std::size_t state = bestState;
    for (std::size_t i = n; i-- > 0;) {
        const int c = fromChoice[i][state];
        const int prev = fromState[i][state];
        if (c < 0 || prev < 0) {
            out.feasible = false;
            return out;
        }
        const std::string& codon = choices[i][static_cast<std::size_t>(c)].codon;
        reversed.append(codon.rbegin(), codon.rend());
        state = static_cast<std::size_t>(prev);
    }
    out.sequence.assign(reversed.rbegin(), reversed.rend());
    return out;
}

}  // namespace

CodonOptimizationResult optimizeCodons(std::string_view protein, const CodonUsageTable& table,
                                       const OptimizeOptions& options) {
    CodonOptimizationResult out;
    out.protein = std::string(protein);
    const GeneticCode& code = builtinGeneticCode(options.geneticCodeId);
    const auto expanded = expandPatterns(options.forbiddenPatterns, options.maxHomopolymer);
    out.forbiddenSites = expanded;
    out.assumptions.push_back("codon usage table: " + table.name() + " (" + table.source() + ")");
    out.assumptions.push_back("NCBI translation table " + std::to_string(options.geneticCodeId));
    out.assumptions.push_back("this is constraint satisfaction, not a prediction: nothing here "
                              "estimates expression level, yield or titre");
    out.assumptions.push_back("the reverse complement of every forbidden pattern is forbidden "
                              "too, because a site on either strand is a site");
    if (options.maxHomopolymer > 0) {
        out.assumptions.push_back("homopolymer runs longer than " +
                                  std::to_string(options.maxHomopolymer) +
                                  " were forbidden as four explicit patterns");
    }

    if (protein.empty()) {
        out.remainingViolations.push_back("no protein sequence was supplied");
        return out;
    }
    std::vector<std::vector<Choice>> choices(protein.size());
    bool usedFloor = false;
    for (std::size_t i = 0; i < protein.size(); ++i) {
        choices[i] = choicesFor(protein[i], table, code, usedFloor);
        if (choices[i].empty()) {
            out.remainingViolations.push_back(
                std::string("residue ") + std::to_string(i + 1) + " ('" + protein[i] +
                "') has no codon in NCBI translation table " +
                std::to_string(options.geneticCodeId));
        }
    }
    if (usedFloor) {
        out.assumptions.push_back("at least one synonymous codon has w = 0 in this table; it was "
                                  "given a 1e-6 floor so it stays representable and last-ranked "
                                  "instead of making every solution infeasible");
    }
    if (!out.remainingViolations.empty()) return out;

    const PatternAutomaton automaton(expanded.empty() ? std::vector<std::string>{"$"} : expanded);
    DpResult best = runDp(out.protein, choices, automaton, 0.0);
    if (!best.feasible) {
        out.remainingViolations.push_back(
            "no synonymous assignment avoids every forbidden pattern" +
            (best.firstInfeasiblePosition >= 0
                 ? "; the search dead-ends at residue " +
                       std::to_string(best.firstInfeasiblePosition + 1)
                 : std::string{}));
        return out;
    }

    if (options.enforceGcWindow) {
        // Lagrangian relaxation: the GC term is separable over codons, so a single
        // multiplier bisection walks GC into the window while the DP keeps the
        // pattern constraint hard. If the window cannot be reached, that is
        // reported - the sequence returned is still constraint-clean.
        const auto gcOf = [](const std::string& s) { return gcPercent(s); };
        double gc = gcOf(best.sequence);
        if (gc < options.minGcPercent || gc > options.maxGcPercent) {
            const double target = gc < options.minGcPercent ? options.minGcPercent
                                                            : options.maxGcPercent;
            double lo = -5.0, hi = 5.0;
            DpResult candidate = best;
            for (int iter = 0; iter < 40; ++iter) {
                const double mid = 0.5 * (lo + hi);
                DpResult trial = runDp(out.protein, choices, automaton, mid);
                if (!trial.feasible) break;
                const double trialGc = gcOf(trial.sequence);
                if (std::abs(trialGc - target) < std::abs(gcOf(candidate.sequence) - target)) {
                    candidate = trial;
                }
                if (trialGc < target) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            best = candidate;
            gc = gcOf(best.sequence);
            out.assumptions.push_back("GC window enforced by a Lagrangian bias on G/C bases; the "
                                      "translation and pattern constraints remained hard");
            if (gc < options.minGcPercent || gc > options.maxGcPercent) {
                out.remainingViolations.push_back(
                    "GC content " + std::to_string(gc) + "% is outside the requested window " +
                    std::to_string(options.minGcPercent) + "-" +
                    std::to_string(options.maxGcPercent) +
                    "%; no synonymous assignment reaches it without breaking a hard constraint");
            }
        }
    }

    // Verify the two hard guarantees against the produced string rather than
    // trusting the DP. This is cheap and it is the whole promise of the type.
    const std::string back = translate(best.sequence, code);
    out.translationPreserved = back == out.protein;
    if (!out.translationPreserved) {
        out.remainingViolations.push_back("internal error: the assignment does not back-translate "
                                          "to the input protein");
        out.optimized.clear();
        return out;
    }
    const auto hits = findForbidden(best.sequence, options.forbiddenPatterns,
                                    options.maxHomopolymer);
    if (!hits.empty()) {
        for (const auto& h : hits) out.remainingViolations.push_back("forbidden pattern " + h);
        out.optimized.clear();
        out.translationPreserved = false;
        return out;
    }
    out.optimized = best.sequence;

    // Before/after metrics use the most-frequent-codon back-translation as the
    // "before", because there is no input DNA: the caller gave a protein.
    std::string naive;
    for (char aa : out.protein) naive += table.mostFrequentCodon(aa);
    out.before = codonMetrics(naive, table);
    out.before.warnings.push_back("'before' is the most-frequent-codon back-translation of the "
                                  "input protein, since the input was a protein and there is no "
                                  "original DNA to measure");
    out.after = codonMetrics(out.optimized, table);
    return out;
}

}  // namespace biocad::bio
