// chem/Formula.cpp - see chem/Formula.h for the data source, the provenance rule
// and the safety scope. This file hard-codes no isotope mass and no atomic
// weight: everything numeric comes from assets/packs/descriptors/isotopes.json.
#include "chem/Formula.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <numeric>
#include <sstream>
#include <unordered_map>

#include <nlohmann/json.hpp>

#include "core/Assets.h"

namespace biocad::chem {
namespace {

constexpr int kSchemaVersion = 1;

// CODATA 2018 electron rest mass in unified atomic mass units. An ion is the
// neutral molecule minus (cation) or plus (anion) this many electron masses, so
// m/z is not simply M/|z| - at high resolution that difference is measurable.
constexpr double kElectronMass = 0.00054857990907;

struct Iso {
    int    massNumber = 0;
    double mass = 0.0;
    double abundance = 0.0;
};

struct Elem {
    std::string      symbol;
    int              z = 0;
    double           standardWeight = 0.0;
    std::vector<Iso> isotopes;       // as shipped, ascending mass number
    int              base = -1;      // index of the most abundant isotope
};

struct Table {
    bool                               ok = false;
    std::string                        note;
    std::string                        source;
    std::unordered_map<int, Elem>      byZ;
    std::unordered_map<std::string, int> zBySymbol;

    const Elem* find(int z) const {
        auto it = byZ.find(z);
        return it == byZ.end() ? nullptr : &it->second;
    }
};

// biocad_chem cannot link biocad_packs (biocad_contracts already depends on
// biocad_chem, so that direction would be a cycle), so the pack is located
// through core::assetRoot(). The environment override lets a test or a harness
// point at the in-tree copy without a packaged layout.
std::filesystem::path findPack() {
    std::error_code ec;
    if (const char* env = std::getenv("BIOCAD_DESCRIPTOR_DIR")) {
        const auto p = std::filesystem::path(env) / "isotopes.json";
        if (std::filesystem::is_regular_file(p, ec)) return p;
    }
    const auto dir = core::assetDir("packs/descriptors");
    if (!dir.empty()) {
        const auto p = dir / "isotopes.json";
        if (std::filesystem::is_regular_file(p, ec)) return p;
    }
    return {};
}

Table loadTable() {
    Table t;
    const auto path = findPack();
    if (path.empty()) {
        t.note = "assets/packs/descriptors/isotopes.json not found (set BIOCAD_DESCRIPTOR_DIR "
                 "to its directory)";
        return t;
    }
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        t.note = path.string() + ": cannot open";
        return t;
    }
    nlohmann::json j;
    try {
        in >> j;
    } catch (const std::exception& e) {
        t.note = path.string() + ": " + e.what();
        return t;
    }
    const int version = j.value("schemaVersion", 0);
    if (version != kSchemaVersion) {
        t.note = path.string() + ": unsupported schemaVersion " + std::to_string(version);
        return t;
    }
    t.source = j.value("source", std::string{});
    if (!j.contains("elements") || !j["elements"].is_array()) {
        t.note = path.string() + ": no elements array";
        return t;
    }
    for (const auto& je : j["elements"]) {
        Elem e;
        e.symbol = je.value("symbol", std::string{});
        e.z = je.value("z", 0);
        e.standardWeight = je.value("standardAtomicWeight", 0.0);
        if (e.symbol.empty() || e.z <= 0 || e.standardWeight <= 0.0) {
            t.note = path.string() + ": incomplete element record near '" + e.symbol + "'";
            return t;
        }
        if (!je.contains("isotopes") || !je["isotopes"].is_array() || je["isotopes"].empty()) {
            t.note = path.string() + ": element " + e.symbol + " has no isotopes";
            return t;
        }
        for (const auto& ji : je["isotopes"]) {
            Iso iso;
            iso.massNumber = ji.value("massNumber", 0);
            iso.mass = ji.value("relativeAtomicMass", 0.0);
            iso.abundance = ji.value("isotopicComposition", 0.0);
            if (iso.massNumber <= 0 || iso.mass <= 0.0) {
                t.note = path.string() + ": element " + e.symbol + " has a malformed isotope";
                return t;
            }
            e.isotopes.push_back(iso);
        }
        std::sort(e.isotopes.begin(), e.isotopes.end(),
                  [](const Iso& a, const Iso& b) { return a.massNumber < b.massNumber; });
        e.base = 0;
        for (std::size_t i = 1; i < e.isotopes.size(); ++i)
            if (e.isotopes[i].abundance > e.isotopes[e.base].abundance) e.base = static_cast<int>(i);
        t.zBySymbol[e.symbol] = e.z;
        t.byZ.emplace(e.z, std::move(e));
    }
    t.ok = true;
    t.note = path.string();
    return t;
}

const Table& table() {
    static const Table t = loadTable();
    return t;
}

std::string citation() {
    static const std::string s = [] {
        const Table& t = table();
        if (!t.source.empty()) return t.source;
        return std::string("NIST Standard Reference Database 144: Atomic Weights and Isotopic "
                           "Compositions with Relative Atomic Masses (isotope table unavailable)");
    }();
    return s;
}

// Conventional valences for the rings-plus-double-bond count. Elements absent
// here (noble gases and the transition metals) have no single conventional
// valence, so an RDBE that included them would be arbitrary; they are excluded
// and named instead.
int conventionalValence(int z) {
    switch (z) {
        case 1: case 3: case 9: case 11: case 17: case 19: case 35: case 53: return 1;
        case 8: case 12: case 16: case 20: case 34: return 2;
        case 5: case 7: case 13: case 15: case 33: return 3;
        case 6: case 14: return 4;
        default: return 0;
    }
}

// --- formula parsing -------------------------------------------------------

using Counts = std::map<std::pair<int, int>, long long>;

void addCount(Counts& c, int z, int massNumber, long long n) { c[{z, massNumber}] += n; }

bool readInt(std::string_view s, std::size_t& i, long long& out) {
    const std::size_t start = i;
    long long v = 0;
    while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i]))) {
        v = v * 10 + (s[i] - '0');
        ++i;
    }
    if (i == start) return false;
    out = v;
    return true;
}

// Longest symbol match at `i` that the table knows: two letters before one, so
// "Cl" never reads as "C" + "l".
int readSymbol(std::string_view s, std::size_t& i, std::string& symbol) {
    const Table& t = table();
    if (i >= s.size() || !std::isupper(static_cast<unsigned char>(s[i]))) return 0;
    std::size_t len = 1;
    while (i + len < s.size() && std::islower(static_cast<unsigned char>(s[i + len]))) ++len;
    for (std::size_t l = len; l >= 1; --l) {
        const std::string cand(s.substr(i, l));
        auto it = t.zBySymbol.find(cand);
        if (it != t.zBySymbol.end()) {
            symbol = cand;
            i += l;
            return it->second;
        }
    }
    symbol.assign(s.substr(i, len));
    return 0;
}

char closerFor(char open) { return open == '(' ? ')' : (open == '[' ? ']' : '}'); }

bool parseSeq(std::string_view s, std::size_t& i, Counts& out, long long mult,
              std::vector<std::string>& warnings, bool digitsAreMultiplier);

// Is s[i] the start of an isotope label the table recognises - "[13C]" or a bare
// "13C"? A bare digit run is a label only when <digits><symbol> is a real
// isotope, otherwise it is a segment multiplier like the 5 of ".5H2O".
bool tryIsotope(std::string_view s, std::size_t& i, Counts& out, long long mult,
                std::vector<std::string>& warnings) {
    const Table& t = table();
    std::size_t j = i;
    bool bracketed = false;
    if (j < s.size() && s[j] == '[') {
        bracketed = true;
        ++j;
    }
    long long a = 0;
    if (!readInt(s, j, a)) return false;
    std::string symbol;
    const int z = readSymbol(s, j, symbol);
    if (z == 0) return false;
    if (bracketed) {
        if (j >= s.size() || s[j] != ']') return false;
        ++j;
    }
    const Elem* e = t.find(z);
    bool known = false;
    if (e)
        for (const auto& iso : e->isotopes)
            if (iso.massNumber == static_cast<int>(a)) known = true;
    if (!known) {
        if (!bracketed) return false;  // a bare "5H" is a multiplier, not an isotope
        warnings.push_back("isotope " + std::to_string(a) + symbol + " is not in the NIST table");
        return false;
    }
    long long n = 1;
    readInt(s, j, n);
    addCount(out, z, static_cast<int>(a), n * mult);
    i = j;
    return true;
}

bool parseSeq(std::string_view s, std::size_t& i, Counts& out, long long mult,
              std::vector<std::string>& warnings, bool digitsAreMultiplier) {
    bool any = false;
    while (i < s.size()) {
        const char c = s[i];
        if (c == ')' || c == ']' || c == '}') return any;
        if (c == '(' || c == '{' || c == '[') {
            if (c == '[' && tryIsotope(s, i, out, mult, warnings)) {
                any = true;
                continue;
            }
            const char close = closerFor(c);
            ++i;
            Counts sub;
            parseSeq(s, i, sub, 1, warnings, false);
            if (i >= s.size() || s[i] != close) {
                warnings.push_back("unbalanced bracket");
                return any;
            }
            ++i;
            long long k = 1;
            readInt(s, i, k);
            for (const auto& [key, n] : sub) addCount(out, key.first, key.second, n * k * mult);
            any = any || !sub.empty();
            continue;
        }
        if (std::isdigit(static_cast<unsigned char>(c))) {
            // A digit run is ambiguous: the 5 of ".5H2O" is a multiplier and the
            // 13 of "13C" is an isotope label. At the start of a segment that
            // follows a hydrate dot the multiplier reading wins, because that is
            // the only place a multiplier is conventionally written; everywhere
            // else a digit run that names a real isotope is a label.
            if (!(digitsAreMultiplier && i == 0) && tryIsotope(s, i, out, mult, warnings)) {
                any = true;
                continue;
            }
            // A digit run where an element was expected multiplies the rest of
            // this segment, which is how a hydrate's ".5H2O" is written.
            long long k = 1;
            readInt(s, i, k);
            const bool sub = parseSeq(s, i, out, mult * k, warnings, false);
            any = any || sub;
            continue;
        }
        if (std::isupper(static_cast<unsigned char>(c))) {
            std::string symbol;
            const int z = readSymbol(s, i, symbol);
            if (z == 0) {
                warnings.push_back("unknown element symbol '" + symbol + "'");
                i += symbol.empty() ? 1 : symbol.size();
                continue;
            }
            long long n = 1;
            readInt(s, i, n);
            addCount(out, z, 0, n * mult);
            any = true;
            continue;
        }
        warnings.push_back(std::string("unexpected character '") + c + "'");
        ++i;
    }
    return any;
}

}  // namespace

std::optional<ParsedFormula> parseFormula(std::string_view text) {
    ParsedFormula f;
    // Whitespace is significant for exactly one thing - separating a charge
    // magnitude from an atom count, as in "SO4 2-" - so it is trimmed here and
    // removed from the body only after the charge has been taken off the end.
    std::string s(text);
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.pop_back();
    if (s.empty()) return std::nullopt;

    // Trailing charge: "+", "--", "3+", "+2", "^2-", " 2-".
    std::size_t end = s.size();
    if (s[end - 1] == '+' || s[end - 1] == '-') {
        const char sign = s[end - 1];
        std::size_t p = end;
        while (p > 0 && s[p - 1] == sign) --p;
        long long magnitude = static_cast<long long>(end - p);
        if (magnitude == 1) {
            std::size_t q = p;
            while (q > 0 && std::isdigit(static_cast<unsigned char>(s[q - 1]))) --q;
            // Digits immediately before the sign are the magnitude only when
            // something separates them from the formula body: "SO4 2-" and
            // "SO4^2-" are a charge, whereas the 4 of "NH4+" is an atom count.
            const bool separated =
                q > 0 && (s[q - 1] == '^' || std::isspace(static_cast<unsigned char>(s[q - 1])));
            if (q < p && (separated || q == 0)) {
                magnitude = std::stoll(s.substr(q, p - q));
                p = q;
            }
        }
        f.charge = static_cast<int>(sign == '+' ? magnitude : -magnitude);
        end = p;
    } else if (std::isdigit(static_cast<unsigned char>(s[end - 1]))) {
        std::size_t p = end;
        while (p > 0 && std::isdigit(static_cast<unsigned char>(s[p - 1]))) --p;
        if (p > 0 && (s[p - 1] == '+' || s[p - 1] == '-')) {
            const long long magnitude = std::stoll(s.substr(p, end - p));
            f.charge = static_cast<int>(s[p - 1] == '+' ? magnitude : -magnitude);
            end = p - 1;
        }
    }
    while (end > 0 && (s[end - 1] == '^' || std::isspace(static_cast<unsigned char>(s[end - 1]))))
        --end;
    s.resize(end);
    {
        std::string body;
        body.reserve(s.size());
        for (const char c : s)
            if (!std::isspace(static_cast<unsigned char>(c))) body.push_back(c);
        s = std::move(body);
    }
    if (s.empty()) return std::nullopt;

    if (!table().ok) f.warnings.push_back(std::string("isotope table unavailable: ") + table().note);

    // Hydrate dots split the formula into independently multiplied segments.
    Counts counts;
    bool any = false;
    std::size_t seg = 0;
    while (seg <= s.size()) {
        std::size_t stop = seg;
        int depth = 0;
        while (stop < s.size()) {
            const char c = s[stop];
            if (c == '(' || c == '[' || c == '{') ++depth;
            else if (c == ')' || c == ']' || c == '}') --depth;
            else if (depth == 0 && (c == '.' || c == '*')) break;
            ++stop;
        }
        const std::string_view part(s.data() + seg, stop - seg);
        if (!part.empty()) {
            std::size_t i = 0;
            any = parseSeq(part, i, counts, 1, f.warnings, seg != 0) || any;
        }
        if (stop >= s.size()) break;
        seg = stop + 1;
    }
    if (!any && counts.empty()) return std::nullopt;

    for (const auto& [key, n] : counts) {
        if (n == 0) continue;
        if (n < 0) f.warnings.push_back("negative atom count for element Z=" + std::to_string(key.first));
        f.terms.push_back({key.first, key.second, static_cast<int>(n)});
    }
    f.ok = table().ok && f.warnings.empty() && !f.terms.empty();
    return f;
}

std::string hillFormula(const ParsedFormula& f) {
    struct Entry {
        std::string label;
        std::string sortKey;
        int         rank = 2;  // 0 = carbon, 1 = hydrogen, 2 = everything else
        int         count = 0;
    };
    std::vector<Entry> entries;
    for (const auto& t : f.terms) {
        const Elem* e = table().find(t.z);
        const std::string symbol = e ? e->symbol : ("Z" + std::to_string(t.z));
        Entry en;
        en.label = t.massNumber ? "[" + std::to_string(t.massNumber) + symbol + "]" : symbol;
        en.sortKey = symbol + (t.massNumber ? std::to_string(t.massNumber) : std::string{});
        en.rank = t.z == 6 ? 0 : (t.z == 1 ? 1 : 2);
        en.count = t.count;
        entries.push_back(std::move(en));
    }
    std::sort(entries.begin(), entries.end(), [](const Entry& a, const Entry& b) {
        if (a.rank != b.rank) return a.rank < b.rank;
        return a.sortKey < b.sortKey;
    });
    std::string out;
    for (const auto& en : entries) {
        out += en.label;
        if (en.count != 1) out += std::to_string(en.count);
    }
    if (f.charge > 0) out += f.charge == 1 ? "+" : std::to_string(f.charge) + "+";
    if (f.charge < 0) out += f.charge == -1 ? "-" : std::to_string(-f.charge) + "-";
    return out;
}

namespace {

const Iso* isotopeOf(const Elem& e, int massNumber) {
    for (const auto& iso : e.isotopes)
        if (iso.massNumber == massNumber) return &iso;
    return nullptr;
}

}  // namespace

double monoisotopicMass(const ParsedFormula& f) {
    double m = 0.0;
    for (const auto& t : f.terms) {
        const Elem* e = table().find(t.z);
        if (!e) return 0.0;
        const Iso* iso = t.massNumber ? isotopeOf(*e, t.massNumber) : &e->isotopes[e->base];
        if (!iso) return 0.0;
        m += iso->mass * t.count;
    }
    return m;
}

double averageMass(const ParsedFormula& f) {
    double m = 0.0;
    for (const auto& t : f.terms) {
        const Elem* e = table().find(t.z);
        if (!e) return 0.0;
        if (t.massNumber) {
            // A labelled position is one isotope, not an isotopic average.
            const Iso* iso = isotopeOf(*e, t.massNumber);
            if (!iso) return 0.0;
            m += iso->mass * t.count;
        } else {
            m += e->standardWeight * t.count;
        }
    }
    return m;
}

int electronCount(const ParsedFormula& f) {
    long long n = 0;
    for (const auto& t : f.terms) n += static_cast<long long>(t.z) * t.count;
    return static_cast<int>(n - f.charge);
}

double ringPlusDoubleBondEquivalents(const ParsedFormula& f, std::vector<std::string>* warnings) {
    double rdbe = 1.0;
    for (const auto& t : f.terms) {
        const int v = conventionalValence(t.z);
        if (v == 0) {
            if (warnings) {
                const Elem* e = table().find(t.z);
                warnings->push_back(std::string("excluded from the rings+double-bond count (no "
                                                "conventional valence): ") +
                                    (e ? e->symbol : ("Z" + std::to_string(t.z))));
            }
            continue;
        }
        rdbe += t.count * (v - 2) / 2.0;
    }
    return rdbe;
}

FormulaMass toFormulaMass(const ParsedFormula& f) {
    FormulaMass out;
    out.formula = hillFormula(f);
    out.charge = f.charge;
    out.warnings = f.warnings;
    for (const auto& t : f.terms) {
        const Elem* e = table().find(t.z);
        const std::string symbol = e ? e->symbol : ("Z" + std::to_string(t.z));
        out.elements.push_back({t.massNumber ? "[" + std::to_string(t.massNumber) + symbol + "]"
                                             : symbol,
                                t.count});
    }
    out.electrons = electronCount(f);
    out.unsaturation = ringPlusDoubleBondEquivalents(f, &out.warnings);

    if (!table().ok) {
        const std::string why = std::string("isotope table: ") + table().note;
        out.monoisotopic = notComputed(why);
        out.average = notComputed(why);
        out.mz = notComputed(why);
        return out;
    }
    // Isotope masses and terrestrial abundances are measured constants, so a sum
    // of them is Measured, not a prediction.
    const double mono = monoisotopicMass(f);
    const double avg = averageMass(f);
    out.monoisotopic = makeQuantity(mono, "Da", 0.0, Provenance::Measured, citation());
    out.average = makeQuantity(avg, "Da", 0.0, Provenance::Measured, citation());
    if (f.charge == 0) {
        out.mz = notComputed("charge is zero");
    } else {
        // An ion is the neutral molecule minus (cation) or plus (anion) charge
        // electrons, so the electron mass is removed/added before dividing.
        const double mz = (mono - f.charge * kElectronMass) / std::abs(f.charge);
        out.mz = makeQuantity(mz, "Da", 0.0, Provenance::Measured,
                              citation() + "; electron mass 0.00054857990907 Da (CODATA 2018)");
    }
    return out;
}

FormulaMass toFormulaMass(std::string_view text) {
    const auto parsed = parseFormula(text);
    if (!parsed) {
        FormulaMass out;
        out.formula = std::string(text);
        out.warnings.push_back("no atoms could be parsed from '" + std::string(text) + "'");
        out.monoisotopic = notComputed("a parseable formula");
        out.average = notComputed("a parseable formula");
        out.mz = notComputed("a parseable formula");
        return out;
    }
    return toFormulaMass(*parsed);
}

namespace {

// One aggregated envelope peak during the convolution: `probMass` is the
// abundance-weighted exact mass sum, so the reported mass of a nominal peak is
// the centroid of the combinations that fall in it, not one arbitrary member.
struct Peak {
    int    shift = 0;
    double prob = 0.0;
    double probMass = 0.0;
};

std::vector<Peak> convolve(const std::vector<Peak>& a, const std::vector<Peak>& b,
                           double minIntensity) {
    std::map<int, Peak> acc;
    for (const auto& pa : a) {
        const double massA = pa.prob > 0 ? pa.probMass / pa.prob : 0.0;
        for (const auto& pb : b) {
            const double massB = pb.prob > 0 ? pb.probMass / pb.prob : 0.0;
            const double prob = pa.prob * pb.prob;
            if (prob <= 0.0) continue;
            Peak& p = acc[pa.shift + pb.shift];
            p.shift = pa.shift + pb.shift;
            p.prob += prob;
            p.probMass += prob * (massA + massB);
        }
    }
    double maxProb = 0.0;
    for (const auto& [_, p] : acc) maxProb = std::max(maxProb, p.prob);
    std::vector<Peak> out;
    for (const auto& [_, p] : acc)
        if (p.prob >= minIntensity * maxProb) out.push_back(p);
    return out;
}

}  // namespace

IsotopeEnvelope isotopeEnvelope(const ParsedFormula& f, double minIntensity) {
    IsotopeEnvelope env;
    env.formula = hillFormula(f);
    env.prunedBelow = minIntensity;
    env.source = citation();
    if (!table().ok || f.terms.empty()) return env;

    std::vector<Peak> total{{0, 1.0, 0.0}};
    for (const auto& t : f.terms) {
        const Elem* e = table().find(t.z);
        if (!e) return {};
        std::vector<Peak> one;
        if (t.massNumber) {
            const Iso* iso = isotopeOf(*e, t.massNumber);
            if (!iso) return {};
            one.push_back({0, 1.0, iso->mass});
        } else {
            const Iso& base = e->isotopes[e->base];
            for (const auto& iso : e->isotopes) {
                if (iso.abundance <= 0.0) continue;
                one.push_back({iso.massNumber - base.massNumber, iso.abundance,
                               iso.abundance * iso.mass});
            }
        }
        // Convolving the single-atom distribution `count` times is the sparse
        // form of the multinomial, and pruning inside the loop keeps the peak
        // list from growing combinatorially.
        for (int k = 0; k < t.count; ++k) total = convolve(total, one, minIntensity);
    }

    double maxProb = 0.0;
    for (const auto& p : total) maxProb = std::max(maxProb, p.prob);
    if (maxProb <= 0.0) return env;
    for (const auto& p : total) {
        IsotopePeak peak;
        peak.mass = p.prob > 0 ? p.probMass / p.prob : 0.0;
        peak.intensity = p.prob / maxProb;
        peak.nominalShift = p.shift;
        env.peaks.push_back(peak);
    }
    std::sort(env.peaks.begin(), env.peaks.end(),
              [](const IsotopePeak& a, const IsotopePeak& b) { return a.mass < b.mass; });
    return env;
}

namespace {

struct Bound {
    int    z = 0;
    int    lo = 0;
    int    hi = 0;
    double mass = 0.0;  // monoisotopic mass of one atom
};

}  // namespace

std::vector<ParsedFormula> findFormulas(
    double targetMass, double toleranceDa,
    const std::map<std::string, std::pair<int, int>>& elementBounds) {
    std::vector<ParsedFormula> out;
    if (!table().ok || elementBounds.empty() || toleranceDa < 0.0) return out;

    std::vector<Bound> bounds;
    for (const auto& [symbol, range] : elementBounds) {
        auto it = table().zBySymbol.find(symbol);
        if (it == table().zBySymbol.end()) continue;
        const Elem* e = table().find(it->second);
        if (!e) continue;
        Bound b;
        b.z = e->z;
        b.lo = std::max(0, range.first);
        b.hi = std::max(b.lo, range.second);
        b.mass = e->isotopes[e->base].mass;
        bounds.push_back(b);
    }
    if (bounds.empty()) return out;
    // Heaviest element first: the mass bound prunes hardest at the top.
    std::sort(bounds.begin(), bounds.end(),
              [](const Bound& a, const Bound& b) { return a.mass > b.mass; });

    // Minimum mass still to be contributed by the remaining elements' lower
    // bounds, used to abandon a branch that can no longer reach the target.
    std::vector<double> minTail(bounds.size() + 1, 0.0);
    for (std::size_t i = bounds.size(); i-- > 0;)
        minTail[i] = minTail[i + 1] + bounds[i].lo * bounds[i].mass;

    std::vector<int> counts(bounds.size(), 0);
    const auto emit = [&](double mass) {
        ParsedFormula f;
        for (std::size_t i = 0; i < bounds.size(); ++i)
            if (counts[i] > 0) f.terms.push_back({bounds[i].z, 0, counts[i]});
        if (f.terms.empty()) return;
        std::sort(f.terms.begin(), f.terms.end(), [](const FormulaTerm& a, const FormulaTerm& b) {
            return a.z == b.z ? a.massNumber < b.massNumber : a.z < b.z;
        });
        f.ok = true;
        (void)mass;
        out.push_back(std::move(f));
    };

    const auto recurse = [&](auto&& self, std::size_t i, double mass) -> void {
        if (i == bounds.size()) {
            if (std::abs(mass - targetMass) <= toleranceDa) emit(mass);
            return;
        }
        for (int n = bounds[i].lo; n <= bounds[i].hi; ++n) {
            const double m = mass + n * bounds[i].mass;
            if (m + minTail[i + 1] > targetMass + toleranceDa) break;
            counts[i] = n;
            self(self, i + 1, m);
        }
        counts[i] = 0;
    };
    recurse(recurse, 0, 0.0);

    std::sort(out.begin(), out.end(), [&](const ParsedFormula& a, const ParsedFormula& b) {
        return std::abs(monoisotopicMass(a) - targetMass) <
               std::abs(monoisotopicMass(b) - targetMass);
    });
    return out;
}

// --- equation balancing ----------------------------------------------------

namespace {

// Exact rational arithmetic, because the whole point of balancing is that the
// coefficients are integers: a floating-point elimination would turn 1/3 into
// something that has to be rounded, and rounding is how a wrong equation gets
// reported as balanced.
struct Rat {
    long long n = 0;
    long long d = 1;
};

Rat makeRat(long long n, long long d) {
    if (d == 0) return {0, 1};
    if (d < 0) {
        n = -n;
        d = -d;
    }
    const long long g = std::gcd(std::abs(n), d);
    return {g ? n / g : n, g ? d / g : d};
}

Rat add(Rat a, Rat b) { return makeRat(a.n * b.d + b.n * a.d, a.d * b.d); }
Rat mul(Rat a, Rat b) { return makeRat(a.n * b.n, a.d * b.d); }
Rat div(Rat a, Rat b) { return makeRat(a.n * b.d, a.d * b.n); }
bool isZero(Rat a) { return a.n == 0; }

}  // namespace

BalancedEquation balanceEquation(const std::vector<std::string>& reactants,
                                 const std::vector<std::string>& products,
                                 const std::vector<double>&      reactantGrams) {
    BalancedEquation out;
    out.reactants = reactants;
    out.products = products;
    out.theoreticalYield = notComputed("reactant amounts in grams");
    out.atomEconomy = notComputed("a balanced equation");

    const std::size_t nr = reactants.size();
    const std::size_t np = products.size();
    if (nr == 0 || np == 0) {
        out.warnings.push_back("an equation needs at least one reactant and one product");
        return out;
    }

    std::vector<ParsedFormula> species;
    species.reserve(nr + np);
    for (const auto& s : reactants) {
        const auto p = parseFormula(s);
        if (!p || p->terms.empty()) {
            out.warnings.push_back("cannot parse reactant '" + s + "'");
            return out;
        }
        species.push_back(*p);
    }
    for (const auto& s : products) {
        const auto p = parseFormula(s);
        if (!p || p->terms.empty()) {
            out.warnings.push_back("cannot parse product '" + s + "'");
            return out;
        }
        species.push_back(*p);
    }

    // Rows: one per distinct (element, isotope label) plus one for charge, so a
    // labelled atom and a charged species are both conserved explicitly.
    std::map<std::pair<int, int>, std::size_t> rowOf;
    for (const auto& sp : species)
        for (const auto& t : sp.terms) rowOf.emplace(std::make_pair(t.z, t.massNumber), rowOf.size());
    const std::size_t nCols = species.size();
    const std::size_t nRows = rowOf.size() + 1;
    std::vector<std::vector<Rat>> a(nRows, std::vector<Rat>(nCols, Rat{0, 1}));
    for (std::size_t c = 0; c < nCols; ++c) {
        const long long sign = c < nr ? 1 : -1;
        for (const auto& t : species[c].terms) {
            const std::size_t r = rowOf.at({t.z, t.massNumber});
            a[r][c] = add(a[r][c], makeRat(sign * t.count, 1));
        }
        a[nRows - 1][c] = add(a[nRows - 1][c], makeRat(sign * species[c].charge, 1));
    }

    // Reduced row echelon form over the rationals.
    std::vector<std::size_t> pivotCol;
    std::size_t row = 0;
    for (std::size_t col = 0; col < nCols && row < nRows; ++col) {
        std::size_t sel = nRows;
        for (std::size_t r = row; r < nRows; ++r)
            if (!isZero(a[r][col])) {
                sel = r;
                break;
            }
        if (sel == nRows) continue;
        std::swap(a[row], a[sel]);
        const Rat p = a[row][col];
        for (std::size_t c = col; c < nCols; ++c) a[row][c] = div(a[row][c], p);
        for (std::size_t r = 0; r < nRows; ++r) {
            if (r == row || isZero(a[r][col])) continue;
            const Rat factor = a[r][col];
            for (std::size_t c = col; c < nCols; ++c)
                a[r][c] = add(a[r][c], mul(factor, makeRat(-a[row][c].n, a[row][c].d)));
        }
        pivotCol.push_back(col);
        ++row;
    }

    std::vector<std::size_t> freeCols;
    for (std::size_t c = 0; c < nCols; ++c)
        if (std::find(pivotCol.begin(), pivotCol.end(), c) == pivotCol.end()) freeCols.push_back(c);
    if (freeCols.size() != 1) {
        out.warnings.push_back(
            freeCols.empty()
                ? "no non-trivial solution: the element-conservation matrix has full column rank, "
                  "so this equation cannot be balanced as written"
                : "the null space has dimension " + std::to_string(freeCols.size()) +
                      ", so the coefficients are not unique; split the equation into independent "
                      "reactions");
        return out;
    }

    const std::size_t freeCol = freeCols[0];
    std::vector<Rat> x(nCols, Rat{0, 1});
    x[freeCol] = Rat{1, 1};
    for (std::size_t r = 0; r < pivotCol.size(); ++r)
        x[pivotCol[r]] = makeRat(-a[r][freeCol].n, a[r][freeCol].d);

    long long lcm = 1;
    for (const auto& v : x) lcm = std::lcm(lcm, v.d);
    std::vector<long long> coeff(nCols, 0);
    for (std::size_t c = 0; c < nCols; ++c) coeff[c] = x[c].n * (lcm / x[c].d);
    long long g = 0;
    for (const auto v : coeff) g = std::gcd(g, std::abs(v));
    if (g > 1)
        for (auto& v : coeff) v /= g;
    // The free coefficient's sign is arbitrary; the physical solution is the one
    // with positive coefficients.
    if (std::count_if(coeff.begin(), coeff.end(), [](long long v) { return v < 0; }) >
        static_cast<long long>(nCols) / 2)
        for (auto& v : coeff) v = -v;

    for (std::size_t c = 0; c < nCols; ++c) {
        if (coeff[c] > 0) continue;
        out.warnings.push_back("coefficient " + std::to_string(coeff[c]) + " for '" +
                               (c < nr ? reactants[c] : products[c - nr]) +
                               "' is not positive, so this equation cannot be balanced as written");
        return out;
    }

    for (std::size_t c = 0; c < nr; ++c) out.reactantCoefficients.push_back(static_cast<int>(coeff[c]));
    for (std::size_t c = nr; c < nCols; ++c) out.productCoefficients.push_back(static_cast<int>(coeff[c]));
    out.balanced = true;

    const std::string yieldSource =
        std::string("stoichiometry of the balanced equation; molar masses from ") + citation();
    std::vector<double> molar(nCols, 0.0);
    for (std::size_t c = 0; c < nCols; ++c) molar[c] = averageMass(species[c]);

    double reactantMass = 0.0;
    for (std::size_t c = 0; c < nr; ++c) reactantMass += coeff[c] * molar[c];
    const double productMass = coeff[nr] * molar[nr];
    if (reactantMass > 0.0)
        out.atomEconomy = makeQuantity(100.0 * productMass / reactantMass, "%", 0.0,
                                       Provenance::Measured, yieldSource);

    if (reactantGrams.empty()) return out;
    if (reactantGrams.size() != nr) {
        out.warnings.push_back("reactantGrams has " + std::to_string(reactantGrams.size()) +
                               " entries for " + std::to_string(nr) + " reactants; amounts ignored");
        return out;
    }
    double limitingRatio = -1.0;
    std::size_t limiting = nr;
    for (std::size_t c = 0; c < nr; ++c) {
        if (reactantGrams[c] <= 0.0 || molar[c] <= 0.0) continue;
        const double ratio = (reactantGrams[c] / molar[c]) / static_cast<double>(coeff[c]);
        if (limiting == nr || ratio < limitingRatio) {
            limitingRatio = ratio;
            limiting = c;
        }
    }
    if (limiting == nr) {
        out.warnings.push_back("no positive reactant amount was supplied");
        return out;
    }
    out.limitingReagent = reactants[limiting];
    out.theoreticalYield = makeQuantity(limitingRatio * coeff[nr] * molar[nr], "g", 0.0,
                                        Provenance::Measured, yieldSource);
    return out;
}

const char* formulaCitation() {
    static const std::string s = citation();
    return s.c_str();
}

bool isotopeTableOk() { return table().ok; }

const char* isotopeTableNote() { return table().note.c_str(); }

}  // namespace biocad::chem
