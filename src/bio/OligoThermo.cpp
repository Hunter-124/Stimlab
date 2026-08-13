#include "bio/OligoThermo.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <mutex>

#include "bio/NucSeq.h"
#include "core/Assets.h"
#include "core/Error.h"

namespace biocad::bio {

namespace {

std::string toDna(std::string_view seq) {
    std::string out;
    out.reserve(seq.size());
    for (char c : seq) {
        char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if (u == 'U') u = 'T';   // the DNA/DNA parameter set is the only one in the pack
        out.push_back(u);
    }
    return out;
}

// The first non-ACGT symbol, or 0. Ambiguity is refused rather than averaged.
char firstDegenerate(const std::string& dna) {
    for (char c : dna) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return c;
    }
    return 0;
}

bool watsonCrick(char a, char b) {
    return (a == 'A' && b == 'T') || (a == 'T' && b == 'A') || (a == 'G' && b == 'C') ||
           (a == 'C' && b == 'G');
}

}  // namespace

const NnParameters::Term* NnParameters::step(char a, char b) const {
    for (const auto& e : steps_) {
        if (e.a == a && e.b == b) return &e.term;
    }
    return nullptr;
}

NnParameters::Term NnParameters::initiation(char terminalBase) const {
    return terminalBase == 'G' || terminalBase == 'C' ? initGC_ : initAT_;
}

NnParameters parseNnParameters(const nlohmann::json& j) {
    const int version = j.value("schemaVersion", 0);
    if (version != 1) {
        throw Error::parse("nn-thermodynamics pack: unsupported schemaVersion " +
                           std::to_string(version));
    }
    NnParameters p;
    p.source_ = j.value("title", std::string{}) + " - " + j.value("description", std::string{});
    p.gasConstant_ = j.value("gasConstant", 0.0);
    p.referenceK_ = j.value("referenceTemperatureK", 0.0);
    if (p.gasConstant_ <= 0 || p.referenceK_ <= 0) {
        throw Error::parse("nn-thermodynamics pack: gasConstant and referenceTemperatureK are "
                           "required and must be positive");
    }
    const auto salt = j.find("saltCorrection");
    if (salt == j.end()) throw Error::parse("nn-thermodynamics pack: no saltCorrection block");
    p.saltCoefficient_ = salt->value("coefficient", 0.0);
    p.saltMethod_ = salt->value("method", std::string{}) + ": " +
                    salt->value("formula", std::string{}) + " [" +
                    salt->value("source", std::string{}) + "]";
    p.magnesiumNote_ = salt->value("magnesium", std::string{});
    if (p.saltCoefficient_ <= 0) {
        throw Error::parse("nn-thermodynamics pack: salt coefficient must be positive");
    }

    const auto nn = j.find("nearestNeighbours");
    if (nn == j.end() || !nn->is_array()) {
        throw Error::parse("nn-thermodynamics pack: no nearestNeighbours");
    }
    for (const auto& e : *nn) {
        // Keys are written in duplex form, "AA/TT", i.e. the top-strand step and
        // its complement. Both orientations are registered from the one entry
        // because reading a duplex 5'->3' on either strand must give the same
        // number; that identity is what makes the ten parameters cover sixteen
        // steps.
        const std::string duplex = e.value("duplex", std::string{});
        if (duplex.size() != 5 || duplex[2] != '/') {
            throw Error::parse("nn-thermodynamics pack: malformed duplex key '" + duplex + "'");
        }
        NnParameters::Term term;
        term.dH = e.value("deltaH", 0.0);
        term.dS = e.value("deltaS", 0.0);
        if (term.dS == 0.0) {
            throw Error::parse("nn-thermodynamics pack: " + duplex + " has no deltaS");
        }
        const char a = duplex[0], b = duplex[1];
        const char ca = complementBase(b), cb = complementBase(a);
        if (!p.step(a, b)) p.steps_.push_back({a, b, term});
        if (!p.step(ca, cb)) p.steps_.push_back({ca, cb, term});
    }
    if (p.steps_.size() != 16) {
        throw Error::parse("nn-thermodynamics pack: the ten unique parameters must expand to 16 "
                           "steps, got " + std::to_string(p.steps_.size()));
    }
    const auto init = j.find("initiation");
    if (init == j.end()) throw Error::parse("nn-thermodynamics pack: no initiation block");
    p.initGC_.dH = init->at("terminalGC").value("deltaH", 0.0);
    p.initGC_.dS = init->at("terminalGC").value("deltaS", 0.0);
    p.initAT_.dH = init->at("terminalAT").value("deltaH", 0.0);
    p.initAT_.dS = init->at("terminalAT").value("deltaS", 0.0);
    const auto sym = j.find("symmetry");
    if (sym == j.end()) throw Error::parse("nn-thermodynamics pack: no symmetry correction");
    p.symmetry_.dH = sym->value("deltaH", 0.0);
    p.symmetry_.dS = sym->value("deltaS", 0.0);
    return p;
}

NnParameters loadNnParameters(const std::filesystem::path& file) {
    std::ifstream in(file, std::ios::binary);
    if (!in) throw Error::io("cannot open nn-thermodynamics pack: " + file.string());
    nlohmann::json j;
    try {
        in >> j;
    } catch (const nlohmann::json::exception& e) {
        throw Error::parse(file.string() + ": " + e.what());
    }
    try {
        return parseNnParameters(j);
    } catch (const Error& e) {
        throw Error::parse(file.string() + ": " + e.message);
    }
}

const NnParameters& builtinNnParameters() {
    static NnParameters params;
    static std::once_flag once;
    static std::string failure;
    std::call_once(once, [] {
        const auto dir = core::assetDir("packs/nucleic");
        if (dir.empty()) {
            failure = "asset root not found";
            return;
        }
        try {
            params = loadNnParameters(dir / "nn-thermodynamics.json");
        } catch (const Error& e) {
            failure = e.message;
        }
    });
    if (params.empty()) throw Error::io("NN thermodynamics unavailable: " + failure);
    return params;
}

bool isSelfComplementary(std::string_view seq) {
    const std::string dna = toDna(seq);
    return !dna.empty() && reverseComplement(dna) == dna;
}

OligoThermo oligoThermo(std::string_view seq, const NnParameters& params,
                        const ThermoOptions& options) {
    OligoThermo out;
    out.sequence = toDna(seq);
    out.naMolar = options.naMolar;
    out.mgMolar = options.mgMolar;
    out.oligoMolar = options.oligoMolar;
    out.dntpMolar = options.dntpMolar;
    out.gcPercent = gcPercent(out.sequence);

    const std::string& dna = out.sequence;
    const auto refuse = [&](const std::string& what) {
        out.deltaH = notComputed(what);
        out.deltaS = notComputed(what);
        out.deltaG37 = notComputed(what);
        out.tm = notComputed(what);
        out.assumptions.push_back("no value was computed: " + what);
        return out;
    };
    if (dna.size() < 2) return refuse("a duplex of at least two base pairs");
    if (const char bad = firstDegenerate(dna)) {
        return refuse(std::string("an unambiguous sequence; '") + bad +
                      "' has no nearest-neighbour parameters and averaging over its expansions "
                      "would be an invention");
    }
    if (options.oligoMolar <= 0.0) return refuse("a positive total strand concentration");
    if (options.naMolar <= 0.0) return refuse("a positive monovalent salt concentration");

    double dH = 0.0, dS = 0.0;
    for (std::size_t i = 0; i + 1 < dna.size(); ++i) {
        const NnParameters::Term* t = params.step(dna[i], dna[i + 1]);
        if (!t) return refuse(std::string("a parameter for step ") + dna[i] + dna[i + 1]);
        dH += t->dH;
        dS += t->dS;
    }
    // One initiation term per end, selected by that end's base pair. The A/T term
    // IS the terminal-AT penalty of the unified set.
    for (char end : {dna.front(), dna.back()}) {
        const auto term = params.initiation(end);
        dH += term.dH;
        dS += term.dS;
    }
    const bool selfComp = isSelfComplementary(dna);
    if (selfComp) {
        const auto term = params.symmetry();
        dH += term.dH;
        dS += term.dS;
    }
    const double x = selfComp ? 1.0 : 4.0;
    const double dSsalt =
        dS + params.saltCoefficient() * static_cast<double>(dna.size() - 1) *
                 std::log(options.naMolar);
    const double tmK = dH * 1000.0 /
                       (dSsalt + params.gasConstant() * std::log(options.oligoMolar / x));
    const double dG37 = dH - params.referenceTemperatureK() * dS / 1000.0;

    const std::string source =
        "SantaLucia unified nearest-neighbour parameter set, from "
        "assets/packs/nucleic/nn-thermodynamics.json; " + params.saltMethod();
    out.deltaH = makeQuantity(dH, "kcal/mol", 0.0, Provenance::Predicted, source);
    out.deltaS = makeQuantity(dS, "cal/(mol*K)", 0.0, Provenance::Predicted,
                              source + " (1 M Na+ standard state)");
    out.deltaG37 = makeQuantity(dG37, "kcal/mol", 0.0, Provenance::Predicted,
                                source + " (dG at 37 C, 1 M Na+ standard state)");
    out.tm = makeQuantity(tmK - 273.15, "degC", 0.0, Provenance::Predicted, source);

    out.assumptions.push_back(
        std::string("Tm = dH / (dS + R ln(Ct/x)) - 273.15 with x = ") + (selfComp ? "1" : "4") +
        (selfComp ? " because the sequence is self-complementary"
                  : " because the sequence is not self-complementary and the two strands are "
                    "assumed present at equal concentration"));
    out.assumptions.push_back("Ct is the total strand concentration, " +
                              std::to_string(options.oligoMolar) + " mol/L");
    out.assumptions.push_back("monovalent salt correction applied to dS: " + params.saltMethod());
    out.assumptions.push_back("dH, dS and dG37 are the 1 M Na+ standard-state values; only Tm "
                              "carries the salt correction");
    if (options.mgMolar > 0.0) {
        out.assumptions.push_back("Mg2+ was supplied (" + std::to_string(options.mgMolar) +
                                  " mol/L) and is UNUSED. " + params.magnesiumNote());
    }
    if (options.dntpMolar > 0.0) {
        out.assumptions.push_back("dNTPs were supplied (" + std::to_string(options.dntpMolar) +
                                  " mol/L) and are UNUSED: dNTP chelation of Mg2+ only matters "
                                  "through a divalent correction, which is not applied");
    }
    out.assumptions.push_back("no error bar is asserted on Tm: the pack carries no published "
                              "prediction RMSE that could be cited here, so the Quantity error "
                              "is 0, meaning unavailable rather than zero");
    return out;
}

OligoThermo oligoThermo(std::string_view seq, const ThermoOptions& options) {
    return oligoThermo(seq, builtinNnParameters(), options);
}

// ------------------------------------------------------- structure scanning

namespace {

// One contiguous run of Watson-Crick pairs, scored with the NN table: the helix
// initiation terms for its two ends plus one step term per adjacent pair. A run
// shorter than two pairs has no step and is not a helix.
double runDeltaG(const std::string& top, const std::string& bottom, std::size_t begin,
                 std::size_t length, const NnParameters& params) {
    if (length < 2) return 0.0;
    double dG = 0.0;
    for (std::size_t k = 0; k + 1 < length; ++k) {
        const NnParameters::Term* t = params.step(top[begin + k], top[begin + k + 1]);
        if (!t) return 0.0;
        dG += NnParameters::deltaG(*t, params.referenceTemperatureK());
    }
    for (char end : {top[begin], top[begin + length - 1]}) {
        dG += NnParameters::deltaG(params.initiation(end), params.referenceTemperatureK());
    }
    (void)bottom;
    return dG;
}

std::string renderDuplex(const std::string& a, const std::string& bReversed, int offset) {
    // Two lines: the top strand 5'->3' and the other strand 3'->5' under it, with
    // the second shifted by `offset` columns. This is a rendering, not a model.
    const int pad = offset < 0 ? -offset : 0;
    std::string top = "5'-" + std::string(static_cast<std::size_t>(pad), ' ') + a + "-3'";
    std::string bottom = "3'-" + std::string(static_cast<std::size_t>(offset > 0 ? offset : 0), ' ') +
                         bReversed + "-5'";
    return top + "\n" + bottom;
}

struct Best { double dG = 0.0; int offset = 0; int position = 0; std::string alignment; };

// Scores every register shift of `b` (read 3'->5') against `a` (5'->3') and keeps
// the most stable contiguous helix found in each. O(n*m), which is the O(n^2) the
// header promises for equal lengths.
std::vector<Best> scanRegisters(const std::string& a, const std::string& b,
                                const NnParameters& params, double maxDeltaG37) {
    std::vector<Best> out;
    const int n = static_cast<int>(a.size());
    const int m = static_cast<int>(b.size());
    std::string bRev(b.rbegin(), b.rend());
    for (int offset = -(m - 1); offset <= n - 1; ++offset) {
        // Column i of `a` pairs with bRev[i - offset].
        std::string paired(a.size(), '.');
        for (int i = 0; i < n; ++i) {
            const int j = i - offset;
            if (j < 0 || j >= m) continue;
            if (watsonCrick(a[static_cast<std::size_t>(i)], bRev[static_cast<std::size_t>(j)])) {
                paired[static_cast<std::size_t>(i)] = '|';
            }
        }
        Best best;
        for (int i = 0; i < n;) {
            if (paired[static_cast<std::size_t>(i)] != '|') {
                ++i;
                continue;
            }
            int run = 0;
            while (i + run < n && paired[static_cast<std::size_t>(i + run)] == '|') ++run;
            const double dG = runDeltaG(a, b, static_cast<std::size_t>(i),
                                        static_cast<std::size_t>(run), params);
            if (dG < best.dG) {
                best.dG = dG;
                best.offset = offset;
                best.position = i;
            }
            i += run;
        }
        if (best.dG <= maxDeltaG37 && best.dG < 0.0) {
            best.alignment = renderDuplex(a, bRev, offset);
            out.push_back(best);
        }
    }
    std::sort(out.begin(), out.end(), [](const Best& x, const Best& y) { return x.dG < y.dG; });
    return out;
}

SecondaryStructure make(const std::string& kind, const Best& b, const std::string& note) {
    SecondaryStructure s;
    s.kind = kind;
    s.position = b.position;
    s.alignment = b.alignment;
    s.deltaG37 = makeQuantity(b.dG, "kcal/mol", 0.0, Provenance::Predicted, note);
    return s;
}

const char* kStructureNote =
    "SantaLucia unified nearest-neighbour terms for the most stable contiguous Watson-Crick "
    "helix in this register, at 1 M Na+ standard state and 37 C. No loop-closure, bulge, "
    "internal-loop or mismatch penalty is applied because none is transcribed in the pack, so "
    "this is the helix's dG37 and is more negative than the whole structure's would be.";

}  // namespace

std::vector<SecondaryStructure> selfDimers(std::string_view seq, const NnParameters& params,
                                           double maxDeltaG37) {
    const std::string dna = toDna(seq);
    if (firstDegenerate(dna)) return {};
    std::vector<SecondaryStructure> out;
    for (const auto& b : scanRegisters(dna, dna, params, maxDeltaG37)) {
        out.push_back(make("self-dimer", b, kStructureNote));
    }
    return out;
}

std::vector<SecondaryStructure> heteroDimers(std::string_view a, std::string_view b,
                                             const NnParameters& params, double maxDeltaG37) {
    const std::string x = toDna(a), y = toDna(b);
    if (firstDegenerate(x) || firstDegenerate(y)) return {};
    std::vector<SecondaryStructure> out;
    for (const auto& hit : scanRegisters(x, y, params, maxDeltaG37)) {
        out.push_back(make("hetero-dimer", hit, kStructureNote));
    }
    return out;
}

std::vector<SecondaryStructure> hairpins(std::string_view seq, const NnParameters& params,
                                         double maxDeltaG37, int minLoop) {
    const std::string dna = toDna(seq);
    std::vector<SecondaryStructure> out;
    if (firstDegenerate(dna) || minLoop < 3) return out;
    const int n = static_cast<int>(dna.size());
    // Every (stem end, loop size) pair, extending the stem outward-in while the
    // bases pair. O(n^2) pairs, each extended in O(n) at worst, which is the same
    // cost class as the dimer scan for the lengths oligos actually have.
    for (int i = 0; i < n; ++i) {
        for (int j = i + minLoop + 1; j < n; ++j) {
            if (!watsonCrick(dna[static_cast<std::size_t>(i)], dna[static_cast<std::size_t>(j)])) {
                continue;
            }
            int stem = 0;
            while (i - stem >= 0 && j + stem < n &&
                   watsonCrick(dna[static_cast<std::size_t>(i - stem)],
                               dna[static_cast<std::size_t>(j + stem)])) {
                ++stem;
            }
            if (stem < 2) continue;
            const int begin = i - stem + 1;
            const double dG = runDeltaG(dna, dna, static_cast<std::size_t>(begin),
                                        static_cast<std::size_t>(stem), params);
            if (dG >= 0.0 || dG > maxDeltaG37) continue;
            Best best;
            best.dG = dG;
            best.position = begin;
            std::string loop = dna.substr(static_cast<std::size_t>(i + 1),
                                          static_cast<std::size_t>(j - i - 1));
            best.alignment = "5'-" + dna.substr(static_cast<std::size_t>(begin),
                                                static_cast<std::size_t>(stem)) +
                             "  loop " + loop + "\n3'-" +
                             std::string(dna.rbegin() + (n - j - stem), dna.rbegin() + (n - j)) +
                             "  (stem " + std::to_string(stem) + " bp)";
            out.push_back(make("hairpin", best, kStructureNote));
        }
    }
    std::sort(out.begin(), out.end(),
              [](const SecondaryStructure& p, const SecondaryStructure& q) {
                  return p.deltaG37.value < q.deltaG37.value;
              });
    return out;
}

std::vector<SecondaryStructure> hairpins(std::string_view seq, double maxDeltaG37) {
    return hairpins(seq, builtinNnParameters(), maxDeltaG37);
}
std::vector<SecondaryStructure> selfDimers(std::string_view seq, double maxDeltaG37) {
    return selfDimers(seq, builtinNnParameters(), maxDeltaG37);
}
std::vector<SecondaryStructure> heteroDimers(std::string_view a, std::string_view b,
                                             double maxDeltaG37) {
    return heteroDimers(a, b, builtinNnParameters(), maxDeltaG37);
}

}  // namespace biocad::bio
