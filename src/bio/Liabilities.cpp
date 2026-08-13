#include "bio/Liabilities.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <regex>
#include <set>

#include <nlohmann/json.hpp>

#include "bio/Imgt.h"
#include "bio/Score.h"
#include "chem/Formula.h"
#include "core/Assets.h"

namespace biocad::bio {
namespace {

using nlohmann::json;

std::string readFile(const std::filesystem::path& p) {
    std::ifstream in(p, std::ios::binary);
    if (!in) return {};
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

// Residue compositions (the residue, i.e. the amino acid minus one water). These
// are COMPOSITIONS, not masses: every number a user sees is monoisotopicMass() or
// averageMass() of the assembled formula over the NIST SRD 144 table.
struct ResidueComposition {
    int c = 0, h = 0, n = 0, o = 0, s = 0;
};

const std::map<char, ResidueComposition>& residueCompositions() {
    static const std::map<char, ResidueComposition> t = {
        {'G', {2, 3, 1, 1, 0}},   {'A', {3, 5, 1, 1, 0}},   {'S', {3, 5, 1, 2, 0}},
        {'P', {5, 7, 1, 1, 0}},   {'V', {5, 9, 1, 1, 0}},   {'T', {4, 7, 1, 2, 0}},
        {'C', {3, 5, 1, 1, 1}},   {'L', {6, 11, 1, 1, 0}},  {'I', {6, 11, 1, 1, 0}},
        {'N', {4, 6, 2, 2, 0}},   {'D', {4, 5, 1, 3, 0}},   {'Q', {5, 8, 2, 2, 0}},
        {'K', {6, 12, 2, 1, 0}},  {'E', {5, 7, 1, 3, 0}},   {'M', {5, 9, 1, 1, 1}},
        {'H', {6, 7, 3, 1, 0}},   {'F', {9, 9, 1, 1, 0}},   {'R', {6, 12, 4, 1, 0}},
        {'Y', {9, 9, 1, 2, 0}},   {'W', {11, 10, 2, 1, 0}},
    };
    return t;
}

// Monosaccharide residue compositions, for the glycoform ladder.
const std::map<std::string, ResidueComposition>& glycanUnits() {
    static const std::map<std::string, ResidueComposition> t = {
        {"HexNAc", {8, 13, 1, 5, 0}},
        {"Hex", {6, 10, 0, 5, 0}},
        {"Fuc", {6, 10, 0, 4, 0}},
        {"NeuAc", {11, 17, 1, 8, 0}},
    };
    return t;
}

std::string formulaOf(const ResidueComposition& r) {
    std::string s;
    auto add = [&s](const char* sym, int n) {
        if (n <= 0) return;
        s += sym;
        if (n != 1) s += std::to_string(n);
    };
    add("C", r.c);
    add("H", r.h);
    add("N", r.n);
    add("O", r.o);
    add("S", r.s);
    return s;
}

// Composition of a polypeptide: residues plus one water.
std::optional<ResidueComposition> compose(const std::string& seq) {
    ResidueComposition t;
    const auto& tbl = residueCompositions();
    for (char c : seq) {
        const auto it = tbl.find(c);
        if (it == tbl.end()) return std::nullopt;
        t.c += it->second.c;
        t.h += it->second.h;
        t.n += it->second.n;
        t.o += it->second.o;
        t.s += it->second.s;
    }
    t.h += 2;   // the terminating water
    t.o += 1;
    return t;
}

double monoOf(const std::string& formula) {
    const auto p = chem::parseFormula(formula);
    return p ? chem::monoisotopicMass(*p) : 0.0;
}
double avgOf(const std::string& formula) {
    const auto p = chem::parseFormula(formula);
    return p ? chem::averageMass(*p) : 0.0;
}

// Mass differences that are elsewhere written as literals. Evaluated once against
// the NIST table so a corrected table corrects them too.
double hydrogenMass() { return monoOf("H"); }
double waterMono() { return monoOf("H2O"); }
double protonMass() {
    // The proton is a hydrogen atom minus one electron, which toFormulaMass()
    // already accounts for when a charge is written on the formula.
    const auto p = chem::parseFormula("H+");
    if (!p) return 0.0;
    return chem::toFormulaMass(*p).mz.value;
}
double deamidationDelta() { return monoOf("O") - monoOf("NH"); }
double carbon13Spacing() { return monoOf("[13C]") - monoOf("C"); }

// ------------------------------------------------------------------- packs

LiabilityPack loadLiabilityPack() {
    LiabilityPack p;
    const auto dir = core::assetDir("packs/biologics");
    if (dir.empty()) {
        p.errors.push_back("no asset tree was found, so liabilities.json could not be read");
        return p;
    }
    const auto file = dir / "liabilities.json";
    const std::string text = readFile(file);
    if (text.empty()) {
        p.errors.push_back("could not read " + file.string());
        return p;
    }
    json j;
    try {
        j = json::parse(text);
    } catch (const std::exception& e) {
        p.errors.push_back(file.string() + ": " + e.what());
        return p;
    }
    if (j.value("schemaVersion", 0) != 1) {
        p.errors.push_back(file.string() + ": unsupported schemaVersion");
        return p;
    }
    p.exposureThreshold = j.value("exposureThreshold", 0.20);
    p.exposureNote = j.value("exposureNote", "");
    for (const auto& r : j["rules"]) {
        LiabilityRule rule;
        rule.ruleId = r.value("ruleId", "");
        rule.motif = r.value("motif", "");
        rule.label = r.value("label", "");
        rule.pattern = r.value("pattern", "");
        rule.risk = r.value("risk", 0);
        rule.requiresExposure = r.value("requiresExposure", false);
        rule.oddCountOnly = r.value("oddCountOnly", false);
        rule.citation = r.value("citation", "");
        // A rule with no citation is exactly the uncitable claim this tree removes.
        if (rule.citation.empty() || rule.pattern.empty()) {
            p.errors.push_back(file.string() + ": rule '" + rule.ruleId +
                               "' is missing a pattern or a citation and was dropped");
            continue;
        }
        p.rules.push_back(std::move(rule));
    }
    p.ok = !p.rules.empty();
    return p;
}

struct Descriptors {
    bool ok = false;
    std::map<char, double> hydropathy;
    std::map<char, std::map<char, double>> diwv;
    std::map<std::string, double> posPk, negPk, cTermPk, nTermPk;
    std::map<std::string, std::string> citations;
    std::vector<std::string> errors;
};

Descriptors loadDescriptors() {
    Descriptors d;
    const auto dir = core::assetDir("packs/biologics");
    if (dir.empty()) {
        d.errors.push_back("no asset tree was found, so descriptors.json could not be read");
        return d;
    }
    const auto file = dir / "descriptors.json";
    const std::string text = readFile(file);
    if (text.empty()) {
        d.errors.push_back("could not read " + file.string());
        return d;
    }
    json j;
    try {
        j = json::parse(text);
    } catch (const std::exception& e) {
        d.errors.push_back(file.string() + ": " + e.what());
        return d;
    }
    for (const auto& [k, v] : j["hydropathy"].items())
        if (!k.empty()) d.hydropathy[k[0]] = v.get<double>();
    for (const auto& [a, row] : j["instabilityDiwv"].items())
        for (const auto& [b, v] : row.items())
            if (!a.empty() && !b.empty()) d.diwv[a[0]][b[0]] = v.get<double>();
    for (const auto& [k, v] : j["pKa"]["positive"].items()) d.posPk[k] = v.get<double>();
    for (const auto& [k, v] : j["pKa"]["negative"].items()) d.negPk[k] = v.get<double>();
    for (const auto& [k, v] : j["pKa"]["cTerminal"].items()) d.cTermPk[k] = v.get<double>();
    for (const auto& [k, v] : j["pKa"]["nTerminal"].items()) d.nTermPk[k] = v.get<double>();
    for (const auto& [k, v] : j["citations"].items()) d.citations[k] = v.get<std::string>();
    d.ok = !d.hydropathy.empty() && !d.diwv.empty() && !d.posPk.empty();
    if (!d.ok) d.errors.push_back(file.string() + ": incomplete scale tables");
    return d;
}

const Descriptors& desc() {
    static const Descriptors d = loadDescriptors();
    return d;
}

std::string citation(const char* key) {
    const auto it = desc().citations.find(key);
    return it == desc().citations.end() ? std::string("assets/packs/biologics/descriptors.json")
                                       : it->second;
}

int countOf(const std::string& s, char c) {
    return static_cast<int>(std::count(s.begin(), s.end(), c));
}

}  // namespace

const LiabilityPack& liabilityPack() {
    static const LiabilityPack p = loadLiabilityPack();
    return p;
}

const DescriptorPack& descriptorPack() {
    static DescriptorPack p = [] {
        DescriptorPack d;
        d.ok = desc().ok;
        d.errors = desc().errors;
        return d;
    }();
    return p;
}

// ---------------------------------------------------------------------------
// Liability scan
// ---------------------------------------------------------------------------

std::vector<SequenceLiability> scanLiabilities(const AbDomain& domain, const Structure* structure,
                                               const std::string& chainId) {
    std::vector<SequenceLiability> out;
    const LiabilityPack& pack = liabilityPack();
    if (!pack.ok) return out;

    const std::string& seq = domain.sequence;
    // IMGT position and region per sequence index, when the domain is numbered. An
    // unnumbered domain still gets scanned - the motif is in the sequence either
    // way - but its flags carry no IMGT position rather than a guessed one.
    std::map<int, const NumberedResidue*> byIndex;
    for (const auto& r : domain.residues) byIndex[r.sequenceIndex] = &r;

    std::vector<double> rel;
    if (structure && !chainId.empty()) rel = relativeSasaFor(domain, *structure, chainId);
    std::map<int, double> relByIndex;
    for (std::size_t i = 0; i < rel.size() && i < domain.residues.size(); ++i)
        if (rel[i] >= 0.0) relByIndex[domain.residues[i].sequenceIndex] = rel[i];

    for (const auto& rule : pack.rules) {
        if (rule.oddCountOnly && countOf(seq, 'C') % 2 == 0) continue;
        std::regex re;
        try {
            re = std::regex(rule.pattern, std::regex::ECMAScript);
        } catch (const std::regex_error&) {
            continue;
        }
        for (auto it = std::sregex_iterator(seq.begin(), seq.end(), re);
             it != std::sregex_iterator(); ++it) {
            const int idx = static_cast<int>(it->position(0));
            SequenceLiability f;
            f.ruleId = rule.ruleId;
            f.motif = rule.motif;
            f.label = rule.label;
            f.sequenceIndex = idx;
            f.citation = rule.citation;
            const auto nr = byIndex.find(idx);
            if (nr != byIndex.end()) {
                f.imgtPosition = positionLabel(*nr->second);
                f.region = nr->second->region;
            }
            const auto rv = relByIndex.find(idx);
            if (rv != relByIndex.end()) {
                f.exposureKnown = true;
                f.relativeSasa = rv->second;
            }
            // An exposure-dependent rule with no coordinates is still REPORTED,
            // with exposureKnown false. Dropping it would hide a real site; calling
            // it exposed would invent one.
            if (rule.requiresExposure && f.exposureKnown &&
                f.relativeSasa < pack.exposureThreshold)
                continue;
            out.push_back(std::move(f));
        }
    }
    std::sort(out.begin(), out.end(), [](const SequenceLiability& a, const SequenceLiability& b) {
        return a.sequenceIndex < b.sequenceIndex;
    });
    return out;
}

// ---------------------------------------------------------------------------
// Developability descriptors
// ---------------------------------------------------------------------------

double netCharge(const std::string& seq, double ph) {
    const Descriptors& d = desc();
    if (!d.ok || seq.empty()) return 0.0;
    auto pk = [&](const std::map<std::string, double>& t, const std::string& k, double dflt) {
        const auto it = t.find(k);
        return it == t.end() ? dflt : it->second;
    };
    double q = 0.0;
    // Termini first: the Bjellqvist set gives residue-specific terminal pKa values.
    {
        const std::string n(1, seq.front()), c(1, seq.back());
        const double nPk = pk(d.nTermPk, n, pk(d.posPk, "Nterm", 7.5));
        const double cPk = pk(d.cTermPk, c, pk(d.negPk, "Cterm", 3.55));
        q += 1.0 / (1.0 + std::pow(10.0, ph - nPk));
        q -= 1.0 / (1.0 + std::pow(10.0, cPk - ph));
    }
    for (const auto& [res, p] : d.posPk) {
        if (res == "Nterm") continue;
        const int n = countOf(seq, res[0]);
        if (n) q += n / (1.0 + std::pow(10.0, ph - p));
    }
    for (const auto& [res, p] : d.negPk) {
        if (res == "Cterm") continue;
        const int n = countOf(seq, res[0]);
        if (n) q -= n / (1.0 + std::pow(10.0, p - ph));
    }
    return q;
}

std::vector<std::pair<double, double>> netChargeCurve(const std::string& seq) {
    std::vector<std::pair<double, double>> out;
    for (int i = 20; i <= 120; ++i) {
        const double ph = i / 10.0;
        out.emplace_back(ph, netCharge(seq, ph));
    }
    return out;
}

namespace {

double isoelectricPoint(const std::string& seq) {
    double lo = 2.0, hi = 13.0;
    for (int i = 0; i < 60; ++i) {
        const double mid = 0.5 * (lo + hi);
        (netCharge(seq, mid) > 0 ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

// Residues in and around the CDRs, which is the region every TAP-like metric is
// defined over. The definition used is stated in the report, because "CDR vicinity"
// is not a standard with one meaning.
std::set<std::string> cdrVicinity(const Structure& s, const std::vector<std::string>& chainIds,
                                  const std::vector<AbDomain>& domains, double radius) {
    std::set<std::string> keys;
    const Model* m = s.model(1);
    if (!m) return keys;
    std::vector<const Atom*> cdrAtoms;
    for (std::size_t k = 0; k < chainIds.size() && k < domains.size(); ++k) {
        const Chain* ch = nullptr;
        for (const auto& c : m->chains)
            if (c.id == chainIds[k]) ch = &c;
        if (!ch || !domains[k].numbered) continue;
        std::set<int> cdrIdx;
        for (const auto& r : domains[k].residues)
            if (r.region.rfind("CDR", 0) == 0) cdrIdx.insert(r.sequenceIndex);
        int poly = 0;
        for (const auto& res : ch->residues) {
            if (res.oneLetter() == 'X') continue;
            if (cdrIdx.count(poly))
                for (const auto& a : res.atoms) cdrAtoms.push_back(&a);
            ++poly;
        }
    }
    if (cdrAtoms.empty()) return keys;
    for (const auto& c : m->chains) {
        for (const auto& res : c.residues) {
            if (res.oneLetter() == 'X') continue;
            bool near = false;
            for (const auto& a : res.atoms) {
                for (const Atom* b : cdrAtoms) {
                    const double dx = a.x - b->x, dy = a.y - b->y, dz = a.z - b->z;
                    if (dx * dx + dy * dy + dz * dz <= radius * radius) {
                        near = true;
                        break;
                    }
                }
                if (near) break;
            }
            if (near)
                keys.insert(c.id + ":" + std::to_string(res.authSeqId) +
                            (res.insertionCode == ' ' ? "" : std::string(1, res.insertionCode)));
        }
    }
    return keys;
}

}  // namespace

DevelopabilityReport developability(const DevelopabilityInput& in) {
    DevelopabilityReport rep;
    std::string all;
    for (const auto& c : in.chains) all += c;
    if (!descriptorPack().ok || all.empty()) {
        const std::string why = all.empty() ? "a sequence"
                                            : "assets/packs/biologics/descriptors.json";
        rep.isoelectricPoint = notComputed(why);
        rep.netChargeAtPh74 = notComputed(why);
        rep.extinctionCoefficient280 = notComputed(why);
        rep.gravy = notComputed(why);
        rep.aliphaticIndex = notComputed(why);
        rep.instabilityIndex = notComputed(why);
        rep.fvChargeSymmetry = notComputed(why);
        rep.tapPsh = notComputed(why);
        rep.tapPpc = notComputed(why);
        rep.tapPnc = notComputed(why);
        rep.tapSfvcsp = notComputed(why);
        for (const auto& e : descriptorPack().errors) rep.warnings.push_back(e);
        return rep;
    }
    const Descriptors& d = desc();

    rep.isoelectricPoint =
        makeQuantity(isoelectricPoint(all), "pH", 0.0, Provenance::Measured, citation("pka"));
    rep.netChargeAtPh74 =
        makeQuantity(netCharge(all, 7.4), "e", 0.0, Provenance::Measured, citation("pka"));

    // eps280. Cystines are PAIRS: an odd cysteine count means one is free, and the
    // free one contributes nothing at 280 nm.
    const int trp = countOf(all, 'W'), tyr = countOf(all, 'Y'), cys = countOf(all, 'C');
    const int cystine = cys / 2;
    rep.extinctionCoefficient280 =
        makeQuantity(5500.0 * trp + 1490.0 * tyr + 125.0 * cystine, "M^-1 cm^-1", 0.0,
                     Provenance::Measured,
                     citation("extinction") + "; assumes all " + std::to_string(cystine) +
                         " cysteine pair(s) are oxidised");

    double gravy = 0;
    int gravyN = 0;
    for (char c : all) {
        const auto it = d.hydropathy.find(c);
        if (it == d.hydropathy.end()) continue;
        gravy += it->second;
        ++gravyN;
    }
    rep.gravy = gravyN ? makeQuantity(gravy / gravyN, "kcal/mol per residue (Kyte-Doolittle)", 0.0,
                                      Provenance::Measured, citation("hydropathy"))
                       : notComputed("any residue with a hydropathy value");

    const double L = static_cast<double>(all.size());
    const double aliphatic = 100.0 * (countOf(all, 'A') / L + 2.9 * countOf(all, 'V') / L +
                                      3.9 * (countOf(all, 'I') + countOf(all, 'L')) / L);
    rep.aliphaticIndex =
        makeQuantity(aliphatic, "", 0.0, Provenance::Measured, citation("aliphatic"));

    // Instability index. Computed per chain, because a dipeptide spanning the joint
    // between two separate chains does not exist.
    {
        double sum = 0;
        double n = 0;
        bool complete = true;
        for (const auto& chain : in.chains) {
            for (std::size_t i = 0; i + 1 < chain.size(); ++i) {
                const auto row = d.diwv.find(chain[i]);
                if (row == d.diwv.end()) {
                    complete = false;
                    continue;
                }
                const auto cell = row->second.find(chain[i + 1]);
                if (cell == row->second.end()) {
                    complete = false;
                    continue;
                }
                sum += cell->second;
            }
            n += static_cast<double>(chain.size());
        }
        rep.instabilityIndex = n > 0 ? makeQuantity(10.0 / n * sum, "", 0.0, Provenance::Measured,
                                                    citation("instability"))
                                     : notComputed("a sequence");
        if (!complete)
            rep.warnings.push_back(
                "the instability index skipped dipeptides containing a non-standard residue, so it "
                "is computed over fewer than L-1 pairs");
    }

    // Fv charge symmetry parameter: the PRODUCT of the two Fv chain net charges at
    // pH 5.5, which is only defined when exactly two chains were supplied.
    if (in.chains.size() == 2) {
        const double qh = netCharge(in.chains[0], 5.5), ql = netCharge(in.chains[1], 5.5);
        rep.fvChargeSymmetry = makeQuantity(qh * ql, "e^2", 0.0, Provenance::Measured,
                                           citation("fvcsp") + "; VH q=" + std::to_string(qh) +
                                               ", VL q=" + std::to_string(ql) + " at pH 5.5");
        rep.assumptions.push_back(
            "Fv charge symmetry treats the first supplied chain as VH and the second as VL, at "
            "pH 5.5.");
    } else {
        rep.fvChargeSymmetry = notComputed("exactly two Fv chains (VH then VL)");
    }

    rep.structureOrigin = in.structureOrigin;
    if (in.structureOrigin.empty() || !in.structure) {
        const char* why =
            "a structure AND its stated origin (crystal entry, or homology model plus protocol)";
        rep.tapPsh = notComputed(why);
        rep.tapPpc = notComputed(why);
        rep.tapPnc = notComputed(why);
        rep.tapSfvcsp = notComputed(why);
        rep.assumptions.push_back(
            "The TAP-style metrics are NotComputed. They are only interpretable next to the origin "
            "of the coordinates: the published thresholds were derived on homology models, and "
            "applying them as a verdict to a crystal structure - or the reverse - is the error "
            "this refusal prevents.");
        return rep;
    }

    // With a structure AND an origin: computed under the protocol printed in every
    // source string. This is NOT the published TAP protocol, so the published
    // thresholds must not be applied to it and no verdict is emitted.
    const double radius = 4.0;
    std::vector<AbDomain> domains;
    for (const auto& c : in.chains) domains.push_back(numberDomain(c));
    const auto vicinity = cdrVicinity(*in.structure, in.chainIds, domains, radius);
    const SasaResult sr = sasa(*in.structure);
    const std::string protocol =
        "in-house protocol, NOT the published TAP protocol (Raybould et al. PNAS 116:4025, 2019): "
        "CDR vicinity = any residue with a heavy atom within " + std::to_string(radius) +
        " A of a CDR heavy atom; hydrophobicity = Kyte-Doolittle; weights = relative side-chain "
        "SASA. Structure origin: " + in.structureOrigin +
        ". The published TAP thresholds are NOT applicable to these numbers and none is shown.";

    double psh = 0, ppc = 0, pnc = 0;
    for (const auto& pr : sr.perResidue) {
        const std::string key =
            pr.key.chainId + ":" + std::to_string(pr.key.authSeqId) +
            (pr.key.insertionCode == ' ' ? "" : std::string(1, pr.key.insertionCode));
        if (!vicinity.count(key)) continue;
        const double rel =
            pr.relative.provenance == Provenance::Measured ? pr.relative.value : 0.0;
        char one = 'X';
        {
            Residue tmp;
            tmp.name = pr.residueName;
            one = tmp.oneLetter();
        }
        const auto h = d.hydropathy.find(one);
        if (h != d.hydropathy.end() && h->second > 0) psh += h->second * rel;
        if (one == 'K' || one == 'R') ppc += rel;
        if (one == 'D' || one == 'E') pnc += rel;
    }
    rep.tapPsh = makeQuantity(psh, "", 0.0, Provenance::Model, protocol);
    rep.tapPpc = makeQuantity(ppc, "", 0.0, Provenance::Model, protocol);
    rep.tapPnc = makeQuantity(pnc, "", 0.0, Provenance::Model, protocol);
    rep.tapSfvcsp = rep.fvChargeSymmetry;
    rep.assumptions.push_back(protocol);
    return rep;
}

// ---------------------------------------------------------------------------
// Mass ladders
// ---------------------------------------------------------------------------

std::string peptideFormula(const std::string& sequence) {
    const auto comp = compose(sequence);
    if (!comp) return {};
    return formulaOf(*comp);
}

std::vector<Glycan> defaultGlycoforms() {
    // The IgG core-fucosylated ladder. Compositions, not masses.
    return {{"G0F", 4, 3, 1, 0}, {"G1F", 4, 4, 1, 0}, {"G2F", 4, 5, 1, 0},
            {"G0", 4, 3, 0, 0},  {"G2F+NeuAc", 4, 5, 1, 1}};
}

namespace {

// Assembles a species: chains, disulfides (each removes two hydrogens), glycans,
// and the pyroGlu / C-terminal-Lys variants.
struct Assembly {
    ResidueComposition comp;
    bool ok = false;
    int  disulfides = 0;
};

Assembly assemble(const std::vector<std::string>& chains, int disulfides,
                  const Glycan* glycan, int glycanCopies, bool pyroGlu, bool clipLys) {
    Assembly a;
    for (auto seq : chains) {
        if (pyroGlu && !seq.empty() && (seq.front() == 'Q' || seq.front() == 'E')) {
            // Pyroglutamate is a cyclisation: Gln loses NH3, Glu loses H2O. Both are
            // expressed as composition changes rather than as a mass literal.
            const auto c = compose(seq);
            if (!c) return a;
            a.comp.c += c->c;
            a.comp.h += c->h - (seq.front() == 'Q' ? 3 : 2);
            a.comp.n += c->n - (seq.front() == 'Q' ? 1 : 0);
            a.comp.o += c->o - (seq.front() == 'Q' ? 0 : 1);
            a.comp.s += c->s;
            continue;
        }
        if (clipLys && !seq.empty() && seq.back() == 'K') seq.pop_back();
        const auto c = compose(seq);
        if (!c) return a;
        a.comp.c += c->c;
        a.comp.h += c->h;
        a.comp.n += c->n;
        a.comp.o += c->o;
        a.comp.s += c->s;
    }
    a.comp.h -= 2 * disulfides;   // one H lost per participating cysteine
    a.disulfides = disulfides;
    if (glycan) {
        const auto& u = glycanUnits();
        auto addUnit = [&](const char* name, int n) {
            const auto it = u.find(name);
            if (it == u.end()) return;
            a.comp.c += it->second.c * n;
            a.comp.h += it->second.h * n;
            a.comp.n += it->second.n * n;
            a.comp.o += it->second.o * n;
            a.comp.s += it->second.s * n;
        };
        for (int k = 0; k < glycanCopies; ++k) {
            addUnit("HexNAc", glycan->hexNAc);
            addUnit("Hex", glycan->hex);
            addUnit("Fuc", glycan->fuc);
            addUnit("NeuAc", glycan->neuAc);
        }
    }
    a.ok = true;
    return a;
}

MassLadderEntry entryOf(const std::string& species, const Assembly& a, const std::string& note) {
    MassLadderEntry e;
    e.species = species;
    e.disulfides = a.disulfides;
    e.note = note;
    if (!a.ok) {
        e.average = notComputed("a sequence containing only the 20 standard residues");
        e.monoisotopic = e.average;
        return e;
    }
    const std::string f = formulaOf(a.comp);
    const double mono = monoOf(f), avg = avgOf(f);
    const std::string src = std::string(chem::formulaCitation()) + "; formula " + f;
    e.average = makeQuantity(avg, "Da", 0.0, Provenance::Measured, src);
    e.monoisotopic = makeQuantity(mono, "Da", 0.0, Provenance::Measured, src);
    return e;
}

}  // namespace

MassLadder massLadder(const MassLadderInput& in) {
    MassLadder out;
    if (!chem::isotopeTableOk()) {
        out.requiredResolvingPower = notComputed(chem::isotopeTableNote());
        out.assumptions.push_back(std::string("no mass could be computed: ") +
                                  chem::isotopeTableNote());
        return out;
    }
    std::vector<std::string> all = in.heavyChains;
    all.insert(all.end(), in.lightChains.begin(), in.lightChains.end());
    const int totalSs = in.interchainDisulfides + in.intrachainDisulfides;

    const std::vector<Glycan>& glycans =
        in.glycoforms.empty() ? defaultGlycoforms() : in.glycoforms;
    const int glycanSites = static_cast<int>(in.heavyChains.size());

    // Above ~10 kDa the monoisotopic peak is not observable on a typical instrument:
    // the isotope envelope of a 150 kDa protein has thousands of peaks and the
    // all-12C one is a vanishing fraction of the total ion current. The average mass
    // is the number that matches the measurement, so it leads.
    const std::string above10k =
        "Above ~10 kDa the AVERAGE mass is the reported one: the all-light-isotope "
        "(monoisotopic) peak of a protein this size carries a negligible fraction of the "
        "envelope and is not resolved on a routine instrument. The monoisotopic column is kept "
        "for peptide-level work.";

    out.entries.push_back(entryOf("intact (aglycosylated)",
                                  assemble(all, totalSs, nullptr, 0, in.pyroGlutamate,
                                           in.cTerminalLysClipped),
                                  above10k));
    out.entries.push_back(
        entryOf("reduced (all disulfides broken)",
                assemble(all, 0, nullptr, 0, in.pyroGlutamate, in.cTerminalLysClipped),
                "Reduction adds one hydrogen per cysteine: " + std::to_string(2 * totalSs) +
                    " H relative to the intact species."));
    out.entries.push_back(entryOf("deglycosylated (PNGase F, no glycan)",
                                  assemble(all, totalSs, nullptr, 0, in.pyroGlutamate,
                                           in.cTerminalLysClipped),
                                  "Identical to the aglycosylated intact mass by construction; "
                                  "PNGase F also converts the sequon Asn to Asp (+0.9840 Da per "
                                  "site), which is NOT applied here."));
    if (!in.heavyChains.empty())
        out.entries.push_back(entryOf("reduced heavy chain",
                                      assemble({in.heavyChains.front()}, 0, nullptr, 0,
                                               in.pyroGlutamate, in.cTerminalLysClipped),
                                      "One heavy chain, no disulfides."));
    if (!in.lightChains.empty())
        out.entries.push_back(entryOf("reduced light chain",
                                      assemble({in.lightChains.front()}, 0, nullptr, 0, false,
                                               false),
                                      "One light chain, no disulfides."));
    if (!in.heavyChains.empty() && !in.lightChains.empty()) {
        // Fab and Fc are DOMAIN subsets and this build does not know where the
        // hinge is in an arbitrary sequence, so they are reported as the chain
        // combinations they are, with that stated.
        out.entries.push_back(
            entryOf("Fab-like (one heavy chain + one light chain)",
                    assemble({in.heavyChains.front(), in.lightChains.front()},
                             in.intrachainDisulfides / 2 + 1, nullptr, 0, in.pyroGlutamate,
                             in.cTerminalLysClipped),
                    "This is one HC plus one LC, not a papain Fab: the hinge cleavage site is not "
                    "inferred from sequence, so a true Fab/Fc split needs the user to supply the "
                    "two fragment sequences."));
    }
    if (in.includeGlycoforms && glycanSites > 0) {
        for (const auto& g : glycans) {
            out.entries.push_back(entryOf(
                "intact + " + g.name + " on " + std::to_string(glycanSites) + " site(s)",
                assemble(all, totalSs, &g, glycanSites, in.pyroGlutamate, in.cTerminalLysClipped),
                "Glycan composition HexNAc" + std::to_string(g.hexNAc) + "Hex" +
                    std::to_string(g.hex) + "Fuc" + std::to_string(g.fuc) + "NeuAc" +
                    std::to_string(g.neuAc) + ", one copy per heavy chain."));
        }
    }
    if (!in.pyroGlutamate)
        out.entries.push_back(entryOf("intact, N-terminal pyroGlu on every heavy chain",
                                      assemble(all, totalSs, nullptr, 0, true,
                                               in.cTerminalLysClipped),
                                      "Gln loses NH3, Glu loses H2O; the composition change, not a "
                                      "mass literal, is what is applied."));
    if (!in.cTerminalLysClipped)
        out.entries.push_back(entryOf("intact, C-terminal Lys clipped",
                                      assemble(all, totalSs, nullptr, 0, in.pyroGlutamate, true),
                                      "One Lys residue removed from every chain that ends in K."));

    // The resolving power that separates a deamidation from an isotope peak. Both
    // deltas come from the same NIST table as the masses.
    const double dDeamid = deamidationDelta(), d13c = carbon13Spacing();
    double mass = 0;
    for (const auto& e : out.entries)
        if (e.species.rfind("intact", 0) == 0 && e.average.provenance == Provenance::Measured) {
            mass = e.average.value;
            break;
        }
    if (mass > 0) {
        const double dm = std::abs(d13c - dDeamid);
        out.requiredResolvingPower = makeQuantity(
            mass / dm, "m/dm (FWHM-equivalent)", 0.0, Provenance::Measured,
            "R = m/dm with m the intact average mass and dm = |(13C-12C) - (O-NH)| = " +
                std::to_string(dm) + " Da, i.e. the gap between the 13C spacing (" +
                std::to_string(d13c) + " Da) and a deamidation (" + std::to_string(dDeamid) +
                " Da). Both deltas are computed from " + chem::formulaCitation() + ".");
        out.assumptions.push_back(
            "A +1 Da shift on an intact antibody is NOT evidence of deamidation: at this mass the "
            "isotope spacing and the deamidation shift differ by " + std::to_string(dm) +
            " Da, so telling them apart needs a resolving power of " +
            std::to_string(static_cast<long long>(mass / dm)) + ".");
    } else {
        out.requiredResolvingPower = notComputed("an intact average mass");
    }
    out.assumptions.push_back(above10k);
    out.assumptions.push_back(
        "Disulfides are applied as " + std::to_string(totalSs) +
        " bond(s) = " + std::to_string(2 * totalSs) + " hydrogens removed from the formula, which "
        "is 2 x " + std::to_string(hydrogenMass()) + " Da per bond from the NIST table.");
    return out;
}

// ---------------------------------------------------------------------------
// Peptide mapping
// ---------------------------------------------------------------------------

std::vector<std::string> proteaseNames() {
    return {"trypsin", "lysc", "gluc", "aspn", "chymotrypsin"};
}

namespace {

// Cleavage sites as indices AFTER which the chain is cut.
std::vector<int> cleavageSites(const std::string& s, const std::string& protease,
                               std::vector<std::string>& warnings) {
    std::vector<int> sites;
    const bool aspn = protease == "aspn";
    for (std::size_t i = 0; i < s.size(); ++i) {
        const char c = s[i];
        const char next = i + 1 < s.size() ? s[i + 1] : '\0';
        bool cut = false;
        if (protease == "trypsin") cut = (c == 'K' || c == 'R') && next != 'P';
        else if (protease == "lysc") cut = c == 'K';
        else if (protease == "gluc") cut = c == 'E';
        else if (protease == "chymotrypsin")
            cut = (c == 'F' || c == 'W' || c == 'Y' || c == 'L') && next != 'P';
        else if (aspn) cut = next == 'D';   // Asp-N cuts N-terminal to Asp
        if (cut && i + 1 < s.size()) sites.push_back(static_cast<int>(i));
    }
    if (protease == "gluc")
        warnings.push_back(
            "Glu-C is applied at Glu only. In phosphate buffer it also cleaves after Asp, which "
            "would change the peptide list; that variant is not assumed here.");
    return sites;
}

}  // namespace

PeptideMap digest(const std::string& chain, const DigestOptions& opts) {
    PeptideMap map;
    map.protease = opts.protease;
    if (chain.empty()) {
        map.warnings.push_back("no sequence was supplied");
        return map;
    }
    // One temporary, bound to a local: begin() and end() of two separate temporaries
    // are iterators into two different vectors, which is undefined behaviour.
    const std::vector<std::string> known = proteaseNames();
    if (std::find(known.begin(), known.end(), opts.protease) == known.end()) {
        map.warnings.push_back("unknown protease '" + opts.protease + "'; known: trypsin, lysc, "
                               "gluc, aspn, chymotrypsin");
        return map;
    }
    const bool massesOk = chem::isotopeTableOk();
    if (!massesOk) map.warnings.push_back(chem::isotopeTableNote());

    std::vector<int> bounds{-1};
    for (int s : cleavageSites(chain, opts.protease, map.warnings)) bounds.push_back(s);
    bounds.push_back(static_cast<int>(chain.size()) - 1);

    const double proton = protonMass(), water = waterMono();
    std::vector<char> covered(chain.size(), 0);

    for (std::size_t i = 0; i + 1 < bounds.size(); ++i) {
        for (int mc = 0; mc <= opts.maxMissedCleavages; ++mc) {
            const std::size_t j = i + 1 + static_cast<std::size_t>(mc);
            if (j >= bounds.size()) break;
            const int begin = bounds[i] + 1, end = bounds[j];
            const int len = end - begin + 1;
            if (len < opts.minLength) continue;
            PeptideFragment f;
            f.sequence = chain.substr(static_cast<std::size_t>(begin),
                                      static_cast<std::size_t>(len));
            f.begin = begin;
            f.end = end;
            f.missedCleavages = mc;
            const std::string formula = peptideFormula(f.sequence);
            if (formula.empty() || !massesOk) {
                f.monoisotopic = notComputed("a peptide of standard residues and the NIST table");
                f.average = f.monoisotopic;
            } else {
                const std::string src = std::string(chem::formulaCitation()) + "; formula " +
                                        formula;
                f.monoisotopic =
                    makeQuantity(monoOf(formula), "Da", 0.0, Provenance::Measured, src);
                f.average = makeQuantity(avgOf(formula), "Da", 0.0, Provenance::Measured, src);
            }
            for (const auto& m : opts.mods) {
                bool applies = m.residues.empty()
                                   ? begin == 0
                                   : f.sequence.find_first_of(m.residues) != std::string::npos;
                if (!applies) continue;
                const double delta = (m.addFormula.empty() ? 0.0 : monoOf(m.addFormula)) -
                                     (m.removeFormula.empty() ? 0.0 : monoOf(m.removeFormula));
                f.modifications.push_back(m.name + " " + (delta >= 0 ? "+" : "") +
                                          std::to_string(delta) + " Da");
            }
            if (opts.bAndYIons && massesOk && !formula.empty()) {
                // b_i = sum of the first i residues + proton; y_i = the last i
                // residues + water + proton. Singly charged.
                double running = 0;
                for (std::size_t k = 0; k + 1 < f.sequence.size(); ++k) {
                    const auto comp = residueCompositions().find(f.sequence[k]);
                    if (comp == residueCompositions().end()) break;
                    running += monoOf(formulaOf(comp->second));
                    f.bIons.push_back(running + proton);
                }
                double rev = 0;
                for (std::size_t k = f.sequence.size(); k-- > 1;) {
                    const auto comp = residueCompositions().find(f.sequence[k]);
                    if (comp == residueCompositions().end()) break;
                    rev += monoOf(formulaOf(comp->second));
                    f.yIons.push_back(rev + water + proton);
                }
            }
            if (mc == 0)
                for (int k = begin; k <= end; ++k) covered[static_cast<std::size_t>(k)] = 1;
            map.peptides.push_back(std::move(f));
        }
    }
    const int cov = static_cast<int>(std::count(covered.begin(), covered.end(), 1));
    map.coveragePct = 100.0 * cov / static_cast<double>(chain.size());
    if (map.coveragePct < 100.0)
        map.warnings.push_back(
            "Sequence coverage is " + std::to_string(map.coveragePct) +
            "% at zero missed cleavages and minimum length " + std::to_string(opts.minLength) +
            ": the short peptides below that length are real, they are just not listed.");
    return map;
}

}  // namespace biocad::bio
