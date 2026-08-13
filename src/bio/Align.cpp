#include "bio/Align.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>

#include "core/Error.h"

namespace biocad::bio {
namespace {

// Well below any reachable score, but far from INT_MIN so that subtracting a gap
// penalty from it cannot overflow.
constexpr int kNeg = std::numeric_limits<int>::min() / 4;

// Traceback states. M = residue aligned to residue, Ix = gap in b (consumes a),
// Iy = gap in a (consumes b).
enum State : unsigned char { kM = 0, kIx = 1, kIy = 2, kStop = 3 };

struct Dp {
    std::size_t rows = 0, cols = 0;
    std::vector<int> m, ix, iy;
    std::vector<unsigned char> fm, fix, fiy;   // predecessor state per cell

    Dp(std::size_t n, std::size_t k)
        : rows(n + 1), cols(k + 1),
          m(rows * cols, kNeg), ix(rows * cols, kNeg), iy(rows * cols, kNeg),
          fm(rows * cols, kStop), fix(rows * cols, kStop), fiy(rows * cols, kStop) {}

    [[nodiscard]] std::size_t at(std::size_t i, std::size_t j) const { return i * cols + j; }
};

// Fills the three Gotoh matrices. `local` switches on the Smith-Waterman zero
// floor in M (a local alignment may restart anywhere) and suppresses the free
// leading-gap initialisation.
void fill(Dp& dp, std::string_view a, std::string_view b, const SubstitutionMatrix& mat,
          GapCost gaps, bool local) {
    const int openExtend = gaps.open + gaps.extend;   // NCBI: gap of k costs open + k*extend
    dp.m[dp.at(0, 0)] = 0;
    if (!local) {
        for (std::size_t i = 1; i < dp.rows; ++i) {
            dp.ix[dp.at(i, 0)] = -(gaps.open + static_cast<int>(i) * gaps.extend);
            dp.fix[dp.at(i, 0)] = (i == 1) ? kM : kIx;
        }
        for (std::size_t j = 1; j < dp.cols; ++j) {
            dp.iy[dp.at(0, j)] = -(gaps.open + static_cast<int>(j) * gaps.extend);
            dp.fiy[dp.at(0, j)] = (j == 1) ? kM : kIy;
        }
    }

    for (std::size_t i = 1; i < dp.rows; ++i) {
        for (std::size_t j = 1; j < dp.cols; ++j) {
            const std::size_t cur = dp.at(i, j);
            const std::size_t diag = dp.at(i - 1, j - 1);
            const int s = mat.score(a[i - 1], b[j - 1]);

            int best = dp.m[diag];
            unsigned char from = kM;
            if (dp.ix[diag] > best) { best = dp.ix[diag]; from = kIx; }
            if (dp.iy[diag] > best) { best = dp.iy[diag]; from = kIy; }
            best = (best <= kNeg) ? kNeg : best + s;
            if (local && best < 0) { best = 0; from = kStop; }
            dp.m[cur] = best;
            dp.fm[cur] = from;

            const std::size_t up = dp.at(i - 1, j);
            const int openIx = (dp.m[up] <= kNeg) ? kNeg : dp.m[up] - openExtend;
            const int extIx = (dp.ix[up] <= kNeg) ? kNeg : dp.ix[up] - gaps.extend;
            if (extIx > openIx) { dp.ix[cur] = extIx; dp.fix[cur] = kIx; }
            else                { dp.ix[cur] = openIx; dp.fix[cur] = kM; }

            const std::size_t left = dp.at(i, j - 1);
            const int openIy = (dp.m[left] <= kNeg) ? kNeg : dp.m[left] - openExtend;
            const int extIy = (dp.iy[left] <= kNeg) ? kNeg : dp.iy[left] - gaps.extend;
            if (extIy > openIy) { dp.iy[cur] = extIy; dp.fiy[cur] = kIy; }
            else                { dp.iy[cur] = openIy; dp.fiy[cur] = kM; }
        }
    }
}

// Walks back from (i, j) in `state`, prepending columns. Stops at the origin
// (global) or at the first zero-valued M cell (local). Returns the start cell.
void traceback(const Dp& dp, std::string_view a, std::string_view b, std::size_t i,
               std::size_t j, unsigned char state, bool local, AlignedRows& rows,
               std::size_t& iStart, std::size_t& jStart) {
    std::string ra, rb;
    while (true) {
        if (i == 0 && j == 0) break;
        if (state == kM) {
            if (i == 0 || j == 0) break;
            if (local && dp.m[dp.at(i, j)] == 0) break;
            ra.push_back(a[i - 1]);
            rb.push_back(b[j - 1]);
            const unsigned char next = dp.fm[dp.at(i, j)];
            --i; --j;
            if (next == kStop) break;
            state = next;
        } else if (state == kIx) {
            if (i == 0) break;
            ra.push_back(a[i - 1]);
            rb.push_back('-');
            const unsigned char next = dp.fix[dp.at(i, j)];
            --i;
            if (next == kStop) break;
            state = next;
        } else {
            if (j == 0) break;
            ra.push_back('-');
            rb.push_back(b[j - 1]);
            const unsigned char next = dp.fiy[dp.at(i, j)];
            --j;
            if (next == kStop) break;
            state = next;
        }
    }
    std::reverse(ra.begin(), ra.end());
    std::reverse(rb.begin(), rb.end());
    rows.a = std::move(ra);
    rows.b = std::move(rb);
    iStart = i;
    jStart = j;
}

AlignmentStats statsOf(AlignedRows& rows, const SubstitutionMatrix& mat) {
    AlignmentStats st;
    rows.midline.assign(rows.a.size(), ' ');
    bool inGapA = false, inGapB = false;
    for (std::size_t k = 0; k < rows.a.size(); ++k) {
        const char ca = rows.a[k];
        const char cb = rows.b[k];
        if (ca == '-' || cb == '-') {
            ++st.gapColumns;
            if (ca == '-' && !inGapA) { ++st.gapOpens; }
            if (cb == '-' && !inGapB) { ++st.gapOpens; }
            inGapA = (ca == '-');
            inGapB = (cb == '-');
            continue;
        }
        inGapA = inGapB = false;
        ++st.alignedColumns;
        if (ca == cb) {
            ++st.identical;
            ++st.positive;
            rows.midline[k] = '|';
        } else if (mat.score(ca, cb) > 0) {
            ++st.positive;
            rows.midline[k] = '+';
        }
    }
    if (st.alignedColumns > 0) {
        const double d = static_cast<double>(st.alignedColumns);
        st.percentIdentity = 100.0 * static_cast<double>(st.identical) / d;
        st.percentSimilarity = 100.0 * static_cast<double>(st.positive) / d;
    }
    return st;
}

int requireInt(const nlohmann::json& j, const char* key, const std::string& where) {
    const auto it = j.find(key);
    if (it == j.end() || !it->is_number()) {
        throw Error::parse(where + ": missing required numeric field \"" + std::string(key) + "\"");
    }
    return it->get<int>();
}

double numberOr(const nlohmann::json& j, const char* key, double fallback) {
    const auto it = j.find(key);
    return (it != j.end() && it->is_number()) ? it->get<double>() : fallback;
}

std::string stringOr(const nlohmann::json& j, const char* key) {
    const auto it = j.find(key);
    return (it != j.end() && it->is_string()) ? it->get<std::string>() : std::string{};
}

}  // namespace

// ------------------------------------------------------------- matrix loading

std::size_t SubstitutionMatrix::index(char c) {
    const auto u = static_cast<unsigned char>(std::toupper(static_cast<unsigned char>(c)));
    return (u < kSlots) ? u : static_cast<std::size_t>('X');
}

int SubstitutionMatrix::score(char a, char b) const {
    const std::size_t ia = index(a);
    const std::size_t ib = index(b);
    const bool ka = known_[ia];
    const bool kb = known_[ib];
    if (ka && kb) return table_[ia * kSlots + ib];
    // Unknown symbol: fall back to the matrix's own X row, as BLAST does, so a
    // stray letter costs the published mismatch rather than a guess.
    if (hasX_) {
        const std::size_t x = static_cast<std::size_t>('X');
        if (ka) return table_[ia * kSlots + x];
        if (kb) return table_[x * kSlots + ib];
        return table_[x * kSlots + x];
    }
    return unknownScore_;
}

const KarlinAltschul* SubstitutionMatrix::statisticsFor(int gapOpen, int gapExtend) const {
    for (const auto& s : stats_) {
        if (s.gapOpen == gapOpen && s.gapExtend == gapExtend) return &s;
    }
    return nullptr;
}

SubstitutionMatrix parseSubstitutionMatrix(const nlohmann::json& j) {
    if (!j.is_object()) throw Error::parse("substitution matrix: document is not an object");
    const std::string where = "substitution matrix";
    const int version = requireInt(j, "schemaVersion", where);
    if (version != 1) {
        throw Error::parse(where + ": unsupported schemaVersion " + std::to_string(version));
    }

    SubstitutionMatrix m;
    m.id_ = stringOr(j, "id");
    m.source_ = stringOr(j, "source");
    m.alphabet_ = stringOr(j, "alphabet");
    m.defaultGapOpen_ = static_cast<int>(numberOr(j, "defaultGapOpen", 11));
    m.defaultGapExtend_ = static_cast<int>(numberOr(j, "defaultGapExtend", 1));

    const auto scores = j.find("scores");
    if (scores == j.end() || !scores->is_object() || scores->empty()) {
        throw Error::parse(where + ": missing \"scores\" object");
    }
    for (auto row = scores->begin(); row != scores->end(); ++row) {
        if (row.key().size() != 1 || !row.value().is_object()) {
            throw Error::parse(where + ": row \"" + row.key() + "\" is not a single-symbol object");
        }
        const std::size_t ia = SubstitutionMatrix::index(row.key()[0]);
        m.known_[ia] = 1;
        for (auto col = row.value().begin(); col != row.value().end(); ++col) {
            if (col.key().size() != 1 || !col.value().is_number()) {
                throw Error::parse(where + ": entry [" + row.key() + "][" + col.key() +
                                   "] is not a numeric single-symbol score");
            }
            const std::size_t ib = SubstitutionMatrix::index(col.key()[0]);
            m.table_[ia * SubstitutionMatrix::kSlots + ib] = col.value().get<int>();
        }
    }
    m.hasX_ = m.known_[static_cast<std::size_t>('X')] != 0;
    m.unknownScore_ = static_cast<int>(numberOr(j, "unknownScore", -4));

    if (const auto st = j.find("statistics"); st != j.end() && st->is_array()) {
        for (const auto& s : *st) {
            KarlinAltschul k;
            k.gapOpen = requireInt(s, "gapOpen", where + ".statistics");
            k.gapExtend = requireInt(s, "gapExtend", where + ".statistics");
            k.lambda = numberOr(s, "lambda", 0.0);
            k.K = numberOr(s, "K", 0.0);
            k.H = numberOr(s, "H", 0.0);
            k.alpha = numberOr(s, "alpha", 0.0);
            k.beta = numberOr(s, "beta", 0.0);
            k.source = stringOr(s, "source");
            if (k.lambda <= 0.0 || k.K <= 0.0) {
                throw Error::parse(where + ".statistics: lambda and K must be positive");
            }
            m.stats_.push_back(std::move(k));
        }
    }
    return m;
}

SubstitutionMatrix loadSubstitutionMatrix(const std::filesystem::path& file) {
    std::ifstream in(file, std::ios::binary);
    if (!in) throw Error::io("cannot open substitution matrix: " + file.string());
    nlohmann::json j;
    try {
        in >> j;
    } catch (const nlohmann::json::exception& e) {
        throw Error::parse(file.string() + ": " + e.what());
    }
    try {
        return parseSubstitutionMatrix(j);
    } catch (const Error& e) {
        throw Error::parse(file.string() + ": " + e.message);
    }
}

// ----------------------------------------------------------------- alignment

GlobalAlignment alignGlobal(std::string_view a, std::string_view b,
                            const SubstitutionMatrix& m, GapCost gaps) {
    if (gaps.open < 0 || gaps.extend <= 0) {
        throw Error::invalidArgument("alignGlobal: gap open must be >= 0 and extend > 0");
    }
    GlobalAlignment out;
    if (a.empty() || b.empty()) {
        // An end-to-end alignment against an empty sequence is all gaps; there is
        // no aligned column, so identity and similarity stay 0 rather than 100.
        out.rows.a.assign(a);
        out.rows.b.assign(b);
        out.rows.a.append(b.size(), '-');
        out.rows.b.insert(out.rows.b.begin(), a.size(), '-');
        const std::size_t len = a.size() + b.size();
        out.score = len == 0 ? 0 : -(gaps.open + static_cast<int>(len) * gaps.extend);
        out.stats = statsOf(out.rows, m);
        return out;
    }

    Dp dp(a.size(), b.size());
    fill(dp, a, b, m, gaps, /*local=*/false);
    const std::size_t last = dp.at(a.size(), b.size());
    unsigned char state = kM;
    int best = dp.m[last];
    if (dp.ix[last] > best) { best = dp.ix[last]; state = kIx; }
    if (dp.iy[last] > best) { best = dp.iy[last]; state = kIy; }
    out.score = best;
    std::size_t iStart = 0, jStart = 0;
    traceback(dp, a, b, a.size(), b.size(), state, /*local=*/false, out.rows, iStart, jStart);
    out.stats = statsOf(out.rows, m);
    return out;
}

LocalAlignment alignLocal(std::string_view a, std::string_view b,
                          const SubstitutionMatrix& m, GapCost gaps) {
    if (gaps.open < 0 || gaps.extend <= 0) {
        throw Error::invalidArgument("alignLocal: gap open must be >= 0 and extend > 0");
    }
    LocalAlignment out;
    if (a.empty() || b.empty()) return out;

    Dp dp(a.size(), b.size());
    fill(dp, a, b, m, gaps, /*local=*/true);

    int best = 0;
    std::size_t bi = 0, bj = 0;
    for (std::size_t i = 1; i <= a.size(); ++i) {
        for (std::size_t j = 1; j <= b.size(); ++j) {
            if (dp.m[dp.at(i, j)] > best) { best = dp.m[dp.at(i, j)]; bi = i; bj = j; }
        }
    }
    out.score = best;
    if (best <= 0) return out;   // no positively scoring segment exists

    std::size_t iStart = 0, jStart = 0;
    traceback(dp, a, b, bi, bj, kM, /*local=*/true, out.rows, iStart, jStart);
    out.aBegin = iStart;
    out.aEnd = bi;
    out.bBegin = jStart;
    out.bEnd = bj;
    out.stats = statsOf(out.rows, m);
    return out;
}

// -------------------------------------------------------------- KA statistics

long computeLengthAdjustment(const KarlinAltschul& ka, std::size_t queryLength, double dbLength,
                             long dbNumSeqs) {
    if (ka.lambda <= 0.0 || ka.K <= 0.0 || dbNumSeqs <= 0) return 0;
    const double m = static_cast<double>(queryLength);
    const double n = dbLength;
    const double N = static_cast<double>(dbNumSeqs);
    const double alphaDLambda = ka.alpha / ka.lambda;
    const double logK = std::log(ka.K);

    // Upper bound: the largest ell for which both effective lengths stay
    // positive, taken from the quadratic (m - ell)(n - N*ell) = max(m,n)/K.
    double ellMax = 0.0;
    {
        const double aq = N;
        const double mb = m * N + n;
        const double c = n * m - std::max(m, n) / ka.K;
        if (c < 0.0) return 0;
        ellMax = 2.0 * c / (mb + std::sqrt(mb * mb - 4.0 * aq * c));
    }

    double ellMin = 0.0;
    double ellNext = 0.0;
    bool converged = false;
    constexpr int kMaxIters = 20;
    for (int i = 1; i <= kMaxIters; ++i) {
        const double ell = ellNext;
        const double ss = (m - ell) * (n - N * ell);
        const double ellBar = alphaDLambda * (logK + std::log(ss)) + ka.beta;
        if (ellBar >= ell) {
            ellMin = ell;
            if (ellBar - ellMin <= 1.0) { converged = true; break; }
            if (ellMin == ellMax) break;
        } else {
            ellMax = ell;
        }
        if (ellMin <= ellBar && ellBar <= ellMax) ellNext = ellBar;
        else if (i == 1)                          ellNext = ellMax;
        else                                      ellNext = 0.5 * (ellMin + ellMax);
    }

    long adjustment = static_cast<long>(ellMin);
    if (converged) {
        const double ell = std::ceil(ellMin);
        if (ell <= ellMax) {
            const double ellBar =
                alphaDLambda * (logK + std::log((m - ell) * (n - N * ell))) + ka.beta;
            if (ellBar >= ell) adjustment = static_cast<long>(ell);
        }
    }
    return adjustment < 0 ? 0 : adjustment;
}

Significance evalueOf(const LocalAlignment& hit, const KarlinAltschul& ka,
                      std::size_t queryLength, double dbLength, long dbNumSeqs) {
    Significance s;
    if (ka.lambda <= 0.0 || ka.K <= 0.0) return s;
    s.lengthAdjustment = computeLengthAdjustment(ka, queryLength, dbLength, dbNumSeqs);
    const double N = static_cast<double>(std::max<long>(dbNumSeqs, 1));
    s.effectiveQueryLength =
        std::max(1.0, static_cast<double>(queryLength) - static_cast<double>(s.lengthAdjustment));
    s.effectiveDbLength = std::max(1.0, dbLength - N * static_cast<double>(s.lengthAdjustment));
    const double score = static_cast<double>(hit.score);
    s.bitScore = (ka.lambda * score - std::log(ka.K)) / std::log(2.0);
    s.evalue = ka.K * s.effectiveQueryLength * s.effectiveDbLength * std::exp(-ka.lambda * score);
    return s;
}

}  // namespace biocad::bio
