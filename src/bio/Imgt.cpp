#include "bio/Imgt.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <string_view>

#include <nlohmann/json.hpp>

#include "bio/Align.h"
#include "bio/Score.h"
#include "core/Assets.h"
#include "core/Error.h"

namespace biocad::bio {
namespace {

using nlohmann::json;

std::filesystem::path packDir() { return core::assetDir("packs/biologics"); }

// The one aligner. A missing matrix must not be replaced by an identity matrix:
// that would produce a plausible germline assignment from nothing.
const SubstitutionMatrix* blosum62() {
    struct Loaded {
        bool ok = false;
        SubstitutionMatrix m;
        Loaded() {
            const auto dir = core::assetDir("packs/matrices");
            if (dir.empty()) return;
            try {
                m = loadSubstitutionMatrix(dir / "blosum62.json");
                ok = true;
            } catch (const Error&) {
            } catch (const std::exception&) {
            }
        }
    };
    static const Loaded loaded;
    return loaded.ok ? &loaded.m : nullptr;
}

std::string readFile(const std::filesystem::path& p) {
    std::ifstream in(p, std::ios::binary);
    if (!in) return {};
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

// IMGT 118 is the first residue of FR4, and FR4 opens with the J-PHE/J-TRP
// followed by Gly-x-Gly. Locating it in the J reference by that motif (rather than
// by a fixed offset) is what lets a 12-residue light J and a 17-residue heavy J
// both report the same position.
int fr4OffsetIn(const std::string& s) {
    for (std::size_t i = 0; i + 3 < s.size(); ++i) {
        if ((s[i] == 'W' || s[i] == 'F') && s[i + 1] == 'G' && s[i + 3] == 'G')
            return static_cast<int>(i);
    }
    return -1;
}

ImgtReference loadReference() {
    ImgtReference r;
    const auto dir = packDir();
    if (dir.empty()) {
        r.errors.push_back("no asset tree was found, so assets/packs/biologics/ could not be read");
        return r;
    }
    const auto file = dir / "imgt-reference.json";
    const std::string text = readFile(file);
    if (text.empty()) {
        r.errors.push_back("could not read " + file.string());
        return r;
    }
    json j;
    try {
        j = json::parse(text);
    } catch (const std::exception& e) {
        r.errors.push_back(file.string() + ": " + e.what());
        return r;
    }
    if (j.value("schemaVersion", 0) != 1) {
        r.errors.push_back(file.string() + ": unsupported schemaVersion " +
                           std::to_string(j.value("schemaVersion", 0)));
        return r;
    }
    r.source = j.value("source", "");
    r.licence = j.value("licence", "");
    r.citation = j.value("citation", "");
    for (const auto& e : j["regions"])
        r.regions.push_back({e.value("name", ""), e.value("first", 0), e.value("last", 0)});
    for (const auto& e : j["anchors"])
        r.anchors.push_back({e.value("position", 0), e.value("accept", ""), e.value("label", "")});
    for (const auto& e : j["vhhHallmark"])
        r.vhhHallmark.push_back({e.value("position", 0), e.value("accept", ""), ""});
    for (const auto& e : j["cdrLayout"]) {
        ImgtCdrLayout l;
        l.region = e.value("region", "");
        l.deletionOrder = e.value("deletionOrder", std::vector<int>{});
        l.insertionAnchors = e.value("insertionAnchors", std::vector<int>{});
        r.cdrLayouts.push_back(std::move(l));
    }
    for (const auto& e : j["vGenes"]) {
        ImgtVGene g;
        g.gene = e.value("gene", "");
        g.geneClass = e.value("geneClass", "");
        g.gapped = e.value("gapped", "");
        g.ungapped.reserve(g.gapped.size());
        g.column.reserve(g.gapped.size());
        for (std::size_t i = 0; i < g.gapped.size(); ++i) {
            if (g.gapped[i] == '.') continue;
            g.ungapped.push_back(g.gapped[i]);
            g.column.push_back(static_cast<int>(i) + 1);
        }
        if (g.ungapped.empty()) continue;
        r.vGenes.push_back(std::move(g));
    }
    for (const auto& e : j["jGenes"]) {
        ImgtJGene g;
        g.gene = e.value("gene", "");
        g.geneClass = e.value("geneClass", "");
        g.sequence = e.value("sequence", "");
        g.fr4Offset = fr4OffsetIn(g.sequence);
        if (g.sequence.empty()) continue;
        r.jGenes.push_back(std::move(g));
    }
    if (r.vGenes.empty() || r.jGenes.empty()) {
        r.errors.push_back(file.string() + ": no usable V or J reference entries");
        return r;
    }
    r.ok = true;
    return r;
}

double bitScoreOf(int rawScore, const SubstitutionMatrix& m) {
    const KarlinAltschul* ka = m.statisticsFor(11, 1);
    if (!ka) return 0.0;
    return (ka->lambda * rawScore - std::log(ka->K)) / std::log(2.0);
}

// Query index -> reference index, from a local alignment's gapped rows. -1 where
// the query residue is an insertion relative to the reference.
struct IndexMap {
    std::vector<int> refOf;      // size = query length
    int              score = 0;
    double           bits = 0;
};

IndexMap mapThrough(const LocalAlignment& al, std::size_t queryLen) {
    IndexMap m;
    m.refOf.assign(queryLen, -1);
    m.score = al.score;
    std::size_t qi = al.aBegin, ri = al.bBegin;
    for (std::size_t c = 0; c < al.rows.a.size(); ++c) {
        const char a = al.rows.a[c], b = al.rows.b[c];
        if (a != '-' && b != '-') {
            if (qi < queryLen) m.refOf[qi] = static_cast<int>(ri);
            ++qi;
            ++ri;
        } else if (a != '-') {
            ++qi;
        } else {
            ++ri;
        }
    }
    return m;
}

std::string setOf(const std::string& gene) {
    const auto star = gene.find('*');
    return star == std::string::npos ? gene : gene.substr(0, star);
}

// Slots of a region, with germline-absent slots and the IMGT deletion order used
// to decide which slot goes empty first, and the two central anchors used to place
// insertion codes when the loop is longer than the region.
struct Slot {
    int position = 0;
    int insertion = 0;   // 0 = none, else 1..n rendered as ".1", ".2"
};

// `trimFromEnd` says which edge a short region loses: FR1 is missing its start
// (an expression construct truncated at the N terminus), FR4 is missing its end
// (a light-chain FR4 is 10 residues, so IMGT 128 simply does not exist), and a
// CDR loses neither because `deletionOrder` covers every one of its slots.
std::vector<Slot> fitRegion(int first, int last, const std::set<int>& germlinePresent,
                            const std::vector<int>& deletionOrder,
                            const std::vector<int>& insertionAnchors, int n,
                            bool trimFromEnd = false) {
    std::vector<int> slots;
    for (int p = first; p <= last; ++p) slots.push_back(p);
    const int capacity = static_cast<int>(slots.size());

    if (n < capacity) {
        std::vector<int> removalOrder;
        // 1. Slots the chosen germline itself does not use.
        for (int p = first; p <= last; ++p)
            if (!germlinePresent.empty() && !germlinePresent.count(p)) removalOrder.push_back(p);
        // 2. The IMGT centre-outward order for this region.
        for (int p : deletionOrder)
            if (p >= first && p <= last &&
                std::find(removalOrder.begin(), removalOrder.end(), p) == removalOrder.end())
                removalOrder.push_back(p);
        // 3. Anything left, from whichever edge this region loses residues at.
        if (trimFromEnd) {
            for (int p = last; p >= first; --p)
                if (std::find(removalOrder.begin(), removalOrder.end(), p) == removalOrder.end())
                    removalOrder.push_back(p);
        } else {
            for (int p = first; p <= last; ++p)
                if (std::find(removalOrder.begin(), removalOrder.end(), p) == removalOrder.end())
                    removalOrder.push_back(p);
        }
        for (int p : removalOrder) {
            if (static_cast<int>(slots.size()) <= n) break;
            slots.erase(std::remove(slots.begin(), slots.end(), p), slots.end());
        }
    }

    std::vector<Slot> out;
    out.reserve(std::max(n, 0));
    if (n <= capacity) {
        for (int i = 0; i < n && i < static_cast<int>(slots.size()); ++i)
            out.push_back({slots[static_cast<std::size_t>(i)], 0});
        return out;
    }

    // Longer than the region: the extra residues become insertion codes on the two
    // central anchors, left ascending and right descending, so every neighbouring
    // position keeps the number it had.
    const int extra = n - capacity;
    int leftAnchor = insertionAnchors.size() > 0 ? insertionAnchors[0] : first;
    int rightAnchor = insertionAnchors.size() > 1 ? insertionAnchors[1] : leftAnchor;
    int leftCount = (extra + 1) / 2, rightCount = extra / 2;
    for (int p : slots) {
        out.push_back({p, 0});
        if (p == leftAnchor)
            for (int k = 1; k <= leftCount; ++k) out.push_back({leftAnchor, k});
        if (p == rightAnchor - 1 && rightAnchor != leftAnchor)
            for (int k = rightCount; k >= 1; --k) out.push_back({rightAnchor, k});
    }
    // A degenerate anchor pair (both equal, or off the slot list) must not silently
    // drop residues: append what is left onto the left anchor.
    int k = leftCount;
    while (static_cast<int>(out.size()) < n) out.push_back({leftAnchor, ++k});
    return out;
}

const ImgtCdrLayout* layoutFor(const ImgtReference& r, const std::string& region) {
    for (const auto& l : r.cdrLayouts)
        if (l.region == region) return &l;
    return nullptr;
}


// ------------------------------------------------------------------ scheme maps

struct Cdr3Rule {
    std::map<int, std::string> leading;
    int                        cdr3First = 0;
    std::vector<std::string>   head;
    std::string                insertionBase;
    std::vector<std::string>   tail;
    std::string                fr4First;
};

struct SchemeTable {
    std::string                scheme;
    std::string                geneClass;
    std::map<int, std::string> positions;
    Cdr3Rule                   cdr3;
};

struct SchemeMaps {
    bool                     ok = false;
    std::vector<SchemeTable> tables;
    // scheme+"/"+geneClass (or scheme+"/*") -> reason it cannot be produced.
    std::map<std::string, std::string> unavailable;
    std::vector<std::string> notes;
    std::vector<std::string> errors;
};

SchemeMaps loadSchemeMaps() {
    SchemeMaps m;
    const auto dir = packDir();
    if (dir.empty()) {
        m.errors.push_back("no asset tree was found, so scheme-maps.json could not be read");
        return m;
    }
    const auto file = dir / "scheme-maps.json";
    const std::string text = readFile(file);
    if (text.empty()) {
        m.errors.push_back("could not read " + file.string());
        return m;
    }
    json j;
    try {
        j = json::parse(text);
    } catch (const std::exception& e) {
        m.errors.push_back(file.string() + ": " + e.what());
        return m;
    }
    if (j.value("schemaVersion", 0) != 1) {
        m.errors.push_back(file.string() + ": unsupported schemaVersion");
        return m;
    }
    for (const auto& n : j.value("notes", std::vector<std::string>{})) m.notes.push_back(n);
    for (const auto& t : j["tables"]) {
        SchemeTable st;
        st.scheme = t.value("scheme", "");
        st.geneClass = t.value("geneClass", "");
        for (const auto& [k, v] : t["positions"].items())
            st.positions[std::stoi(k)] = v.get<std::string>();
        const auto& c = t["cdr3Rule"];
        for (const auto& [k, v] : c["leading"].items())
            st.cdr3.leading[std::stoi(k)] = v.get<std::string>();
        st.cdr3.cdr3First = c.value("cdr3First", 0);
        st.cdr3.head = c.value("head", std::vector<std::string>{});
        st.cdr3.insertionBase = c.value("insertionBase", "");
        st.cdr3.tail = c.value("tail", std::vector<std::string>{});
        st.cdr3.fr4First = c.value("fr4First", "");
        m.tables.push_back(std::move(st));
    }
    for (const auto& u : j.value("unavailable", json::array()))
        m.unavailable[u.value("scheme", "") + "/" + u.value("geneClass", "*")] =
            u.value("reason", "");
    m.ok = !m.tables.empty();
    return m;
}

const SchemeMaps& schemeMaps() {
    static const SchemeMaps m = loadSchemeMaps();
    return m;
}

const char* schemeId(NumberingScheme s) {
    switch (s) {
        case NumberingScheme::Imgt:    return "imgt";
        case NumberingScheme::Kabat:   return "kabat";
        case NumberingScheme::Chothia: return "chothia";
        case NumberingScheme::Martin:  return "martin";
        case NumberingScheme::Aho:     return "aho";
    }
    return "imgt";
}

std::string geneClassOf(AbChainType c) {
    switch (c) {
        case AbChainType::HeavyVh:
        case AbChainType::Vhh:          return "IGHV";
        case AbChainType::LightVKappa:  return "IGKV";
        case AbChainType::LightVLambda: return "IGLV";
        case AbChainType::TcrBeta:      return "TRBV";
        case AbChainType::TcrAlpha:     return "TRAV";
        case AbChainType::Unknown:      break;
    }
    return "";
}

}  // namespace

const ImgtReference& imgtReference() {
    static const ImgtReference r = loadReference();
    return r;
}

std::string positionLabel(const NumberedResidue& r) {
    std::string s = std::to_string(r.position);
    if (!r.insertionCode.empty()) s += "." + r.insertionCode;
    return s;
}

AbDomain numberDomain(const std::string& sequence, const NumberingOptions& opts) {
    AbDomain d;
    d.scheme = NumberingScheme::Imgt;

    std::string q;
    for (char c : sequence) {
        const char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if (u >= 'A' && u <= 'Z') q.push_back(u);
    }
    d.sequence = q;

    const ImgtReference& ref = imgtReference();
    if (!ref.ok) {
        for (const auto& e : ref.errors)
            d.anchorFailures.push_back("IMGT reference unavailable: " + e);
        return d;
    }
    const SubstitutionMatrix* mat = blosum62();
    if (!mat) {
        d.anchorFailures.push_back(
            "assets/packs/matrices/blosum62.json could not be loaded, so no germline alignment "
            "was possible");
        return d;
    }
    if (q.size() < 60) {
        d.anchorFailures.push_back("the sequence is " + std::to_string(q.size()) +
                                   " residues; a V-DOMAIN is ~95-130 and nothing shorter than 60 "
                                   "is numbered");
        return d;
    }

    // --- best and runner-up germline V ------------------------------------------
    const ImgtVGene* best = nullptr;
    IndexMap bestMap;
    double bestBits = 0;
    std::string bestSet, runnerSet;
    double runnerBits = 0;
    for (const auto& g : ref.vGenes) {
        const LocalAlignment al = alignLocal(q, g.ungapped, *mat, GapCost{11, 1});
        const double bits = bitScoreOf(al.score, *mat);
        if (bits > bestBits) {
            if (setOf(g.gene) != bestSet) {
                runnerSet = bestSet;
                runnerBits = bestBits;
            }
            best = &g;
            bestBits = bits;
            bestSet = setOf(g.gene);
            bestMap = mapThrough(al, q.size());
        } else if (bits > runnerBits && setOf(g.gene) != bestSet) {
            runnerSet = setOf(g.gene);
            runnerBits = bits;
        }
    }
    if (!best) {
        d.anchorFailures.push_back("no germline V reference could be aligned");
        return d;
    }
    d.closestGermlineSet = bestSet;
    d.runnerUpGermlineSet = runnerSet;
    d.bestBitScore = bestBits;
    d.runnerUpBitScore = runnerBits;
    d.vGene = best->gene;

    // Chain type, from the winning gene class plus the camelid FR2 hallmark.
    if (best->geneClass == "IGHV") d.chain = AbChainType::HeavyVh;
    else if (best->geneClass == "IGKV") d.chain = AbChainType::LightVKappa;
    else if (best->geneClass == "IGLV") d.chain = AbChainType::LightVLambda;
    else if (best->geneClass == "TRBV") d.chain = AbChainType::TcrBeta;
    else if (best->geneClass == "TRAV") d.chain = AbChainType::TcrAlpha;

    // Position of every query residue under the winning germline's IMGT columns.
    std::vector<int> imgtOf(q.size(), 0);
    for (std::size_t i = 0; i < q.size(); ++i) {
        const int ri = bestMap.refOf[i];
        if (ri >= 0 && ri < static_cast<int>(best->column.size()))
            imgtOf[i] = best->column[static_cast<std::size_t>(ri)];
    }
    auto indexOfPosition = [&](int p) -> int {
        for (std::size_t i = 0; i < imgtOf.size(); ++i)
            if (imgtOf[i] == p) return static_cast<int>(i);
        return -1;
    };

    if (bestBits < opts.minVBitScore) {
        d.anchorFailures.push_back(
            "the best germline V reference (" + bestSet + ") scored only " +
            std::to_string(static_cast<int>(bestBits)) + " bits, below the " +
            std::to_string(static_cast<int>(opts.minVBitScore)) +
            "-bit floor; this is not a V-DOMAIN the reference pack covers");
        return d;
    }

    // A T-cell receptor is a V-set domain and WILL satisfy the anchors, so it has
    // to be refused by identity, before the anchor test, or it would number
    // cleanly as an antibody.
    if (d.chain == AbChainType::TcrBeta || d.chain == AbChainType::TcrAlpha) {
        d.anchorFailures.push_back(
            "the closest germline set is " + bestSet + " (" + best->geneClass +
            "), a T-cell receptor V gene, not an immunoglobulin; this module numbers "
            "immunoglobulin V-DOMAINs only and refuses to apply IG CDR definitions to a TCR");
        return d;
    }

    // --- anchors ---------------------------------------------------------------
    int idx23 = -1, idx41 = -1, idx89 = -1, idx104 = -1, idx118 = -1;
    for (const auto& a : ref.anchors) {
        if (a.position == 118) continue;   // located through the J reference below
        const int i = indexOfPosition(a.position);
        if (i < 0) {
            d.anchorFailures.push_back("IMGT " + std::to_string(a.position) + " (" + a.label +
                                       ") has no aligned residue");
            continue;
        }
        if (a.accept.find(q[static_cast<std::size_t>(i)]) == std::string::npos) {
            d.anchorFailures.push_back("IMGT " + std::to_string(a.position) + " (" + a.label +
                                       ") is " + std::string(1, q[static_cast<std::size_t>(i)]) +
                                       ", expected one of " + a.accept);
            continue;
        }
        if (a.position == 23) idx23 = i;
        else if (a.position == 41) idx41 = i;
        else if (a.position == 89) idx89 = i;
        else if (a.position == 104) idx104 = i;
    }

    // IMGT 118 comes from the J reference: it is the J-PHE/J-TRP that opens FR4.
    // The matched J also fixes how LONG FR4 is - 11 residues (118-128) for a heavy
    // chain, 10 for a light one - which is what stops the first constant-domain
    // residue being swallowed as IMGT 128.
    int fr4Length = 11;
    {
        double jBits = 0;
        for (const auto& g : ref.jGenes) {
            if (g.fr4Offset < 0) continue;
            if (best->geneClass == "IGHV" && g.geneClass != "IGHJ") continue;
            if (best->geneClass == "IGKV" && g.geneClass != "IGKJ") continue;
            if (best->geneClass == "IGLV" && g.geneClass != "IGLJ") continue;
            const LocalAlignment al = alignLocal(q, g.sequence, *mat, GapCost{11, 1});
            const double bits = bitScoreOf(al.score, *mat);
            if (bits <= jBits) continue;
            const IndexMap jm = mapThrough(al, q.size());
            int cand = -1;
            for (std::size_t i = 0; i < q.size(); ++i)
                if (jm.refOf[i] == g.fr4Offset) cand = static_cast<int>(i);
            if (cand < 0) continue;
            jBits = bits;
            idx118 = cand;
            fr4Length = static_cast<int>(g.sequence.size()) - g.fr4Offset;
            d.jGene = g.gene;
        }
        const std::string accept118 = "FW";
        if (idx118 < 0) {
            d.anchorFailures.push_back(
                "IMGT 118 (J-PHE or J-TRP) could not be located: no J reference aligned onto an "
                "FR4 start in this sequence");
        } else if (accept118.find(q[static_cast<std::size_t>(idx118)]) == std::string::npos) {
            d.anchorFailures.push_back("IMGT 118 (J-PHE or J-TRP) is " +
                                       std::string(1, q[static_cast<std::size_t>(idx118)]) +
                                       ", expected F or W");
            idx118 = -1;
        }
    }

    if (!d.anchorFailures.empty()) return d;
    if (!(idx23 < idx41 && idx41 < idx89 && idx89 < idx104 && idx104 < idx118)) {
        d.anchorFailures.push_back("the five anchors are not in sequence order, so the alignment "
                                   "register is wrong");
        return d;
    }

    // Camelid hallmark: reported, never rendered as a species.
    if (d.chain == AbChainType::HeavyVh) {
        int hits = 0, tested = 0;
        for (const auto& h : ref.vhhHallmark) {
            const int i = indexOfPosition(h.position);
            if (i < 0) continue;
            ++tested;
            if (h.accept.find(q[static_cast<std::size_t>(i)]) != std::string::npos) ++hits;
        }
        if (tested >= 3 && hits >= 3) {
            d.chain = AbChainType::Vhh;
            d.warnings.push_back(
                "FR2 carries " + std::to_string(hits) + " of " + std::to_string(tested) +
                " camelid VHH hallmark residues (IMGT 42/49/50/52), so this is treated as a "
                "single-domain VHH. That is a sequence hallmark, NOT a species identification.");
        }
    }

    // --- region-by-region assignment -------------------------------------------
    std::set<int> present;
    for (int p : best->column) present.insert(p);
    auto presentIn = [&](int a, int b) {
        std::set<int> s;
        for (int p : present)
            if (p >= a && p <= b) s.insert(p);
        return s;
    };
    auto countPresent = [&](int a, int b) { return static_cast<int>(presentIn(a, b).size()); };

    // FR3 runs 66..104 and is the only framework whose germline gap pattern varies
    // (73 always, 81/82 in the light classes), so its start is derived from the 89
    // anchor and the winning germline's own slot count rather than a fixed offset.
    const int fr3Slots66to89 = countPresent(66, 89);
    const int idx66 = idx89 - (fr3Slots66to89 - 1);
    const int idx39 = idx41 - 2;
    const int idx55 = idx41 + 14;

    if (idx39 < 0 || idx66 <= idx55 || idx23 + 4 > idx39) {
        d.anchorFailures.push_back(
            "the framework spacing implied by the anchors is impossible (FR2 or FR3 would have a "
            "negative length), so no numbering is produced");
        return d;
    }

    struct Block {
        std::string name;
        int first, last;     // IMGT position range
        int begin, end;      // query index range, inclusive
    };
    const std::vector<Block> blocks = {
        {"FR1", 1, 26, 0, idx23 + 3},
        {"CDR1", 27, 38, idx23 + 4, idx39 - 1},
        {"FR2", 39, 55, idx39, idx55},
        {"CDR2", 56, 65, idx55 + 1, idx66 - 1},
        {"FR3", 66, 104, idx66, idx104},
        {"CDR3", 105, 117, idx104 + 1, idx118 - 1},
        {"FR4", 118, 128, idx118,
         std::min<int>(static_cast<int>(q.size()) - 1, idx118 + fr4Length - 1)},
    };

    for (const auto& b : blocks) {
        const int n = b.end - b.begin + 1;
        if (n < 0) {
            d.anchorFailures.push_back(b.name + " would have a negative length");
            return d;
        }
        const ImgtCdrLayout* lay = layoutFor(ref, b.name);
        const std::vector<Slot> slots =
            fitRegion(b.first, b.last, presentIn(b.first, b.last),
                      lay ? lay->deletionOrder : std::vector<int>{},
                      lay ? lay->insertionAnchors : std::vector<int>{}, n, b.name == "FR4");
        for (int k = 0; k < n; ++k) {
            NumberedResidue nr;
            const Slot s = slots[static_cast<std::size_t>(k)];
            nr.position = s.position;
            if (s.insertion > 0) nr.insertionCode = std::to_string(s.insertion);
            nr.sequenceIndex = b.begin + k;
            nr.aminoAcid = q[static_cast<std::size_t>(nr.sequenceIndex)];
            nr.region = b.name;
            d.residues.push_back(std::move(nr));
        }
    }

    const int trailing = static_cast<int>(q.size()) - 1 - blocks.back().end;
    if (trailing > 0 && opts.warnOnTrailingResidues)
        d.warnings.push_back(std::to_string(trailing) +
                             " residue(s) follow IMGT 128 and were NOT numbered: they are a "
                             "constant domain or a tag, not part of the V-DOMAIN.");

    for (const char* reg : {"CDR1", "CDR2", "CDR3"}) {
        int n = 0;
        for (const auto& r : d.residues)
            if (r.region == reg) ++n;
        d.cdrLengths.push_back(n);
    }
    d.numbered = true;
    return d;
}

AbDomain convertScheme(const AbDomain& domain, NumberingScheme to) {
    AbDomain out = domain;
    out.scheme = to;
    if (!domain.numbered) return out;

    if (to == NumberingScheme::Imgt && domain.scheme == NumberingScheme::Imgt) return out;
    if (to == NumberingScheme::Imgt) {
        // The IMGT positions were never discarded: a converted domain keeps its
        // region labels and sequence indices, so going back is a renumber from the
        // canonical reference, not an inverse table.
        const AbDomain again = numberDomain(domain.sequence);
        out.residues = again.residues;
        out.numbered = again.numbered;
        out.anchorFailures = again.anchorFailures;
        return out;
    }

    const SchemeMaps& maps = schemeMaps();
    const std::string cls = geneClassOf(domain.chain);
    if (!maps.ok) {
        out.numbered = false;
        out.residues.clear();
        for (const auto& e : maps.errors) out.warnings.push_back("scheme table unavailable: " + e);
        return out;
    }
    const std::string id = schemeId(to);
    for (const auto& key : {id + "/" + cls, id + "/*"}) {
        const auto it = maps.unavailable.find(key);
        if (it == maps.unavailable.end()) continue;
        out.numbered = false;
        out.residues.clear();
        out.warnings.push_back("the " + id + " table for " + cls +
                               " is not shipped and this conversion is REFUSED: " + it->second +
                               ". A guessed correspondence would renumber every residue "
                               "plausibly and wrongly.");
        return out;
    }
    const SchemeTable* tbl = nullptr;
    for (const auto& t : maps.tables)
        if (t.scheme == id && t.geneClass == cls) tbl = &t;
    if (!tbl) {
        out.numbered = false;
        out.residues.clear();
        out.warnings.push_back("no " + id + " position table for gene class " + cls +
                               " is present in assets/packs/biologics/scheme-maps.json");
        return out;
    }

    // CDR3 is length-addressed, not position-addressed: Kabat numbers H3 95-102
    // with insertions at 100 and L3 89-97 with insertions at 95, so the labels
    // depend on how many residues the loop actually has.
    std::vector<std::string> cdr3Labels;
    {
        int n = 0;
        for (const auto& r : domain.residues)
            if (r.region == "CDR3" && r.position >= tbl->cdr3.cdr3First) ++n;
        const int head = static_cast<int>(tbl->cdr3.head.size());
        const int tail = static_cast<int>(tbl->cdr3.tail.size());
        for (int i = 0; i < n; ++i) {
            if (i < head && i < n - tail) {
                cdr3Labels.push_back(tbl->cdr3.head[static_cast<std::size_t>(i)]);
            } else if (i >= n - tail) {
                const int k = tail - (n - i);
                cdr3Labels.push_back(tbl->cdr3.tail[static_cast<std::size_t>(k)]);
            } else {
                cdr3Labels.push_back(tbl->cdr3.insertionBase +
                                     std::string(1, static_cast<char>('A' + (i - head))));
            }
        }
    }

    int fr4 = std::stoi(tbl->cdr3.fr4First);
    int cdr3Seen = 0;
    out.residues.clear();
    for (const auto& r : domain.residues) {
        NumberedResidue nr = r;
        std::string label;
        if (r.region == "FR4") {
            label = std::to_string(fr4 + (r.position - 118));
        } else if (r.region == "CDR3" || r.position >= 105) {
            const auto lead = tbl->cdr3.leading.find(r.position);
            if (lead != tbl->cdr3.leading.end() && r.insertionCode.empty()) {
                label = lead->second;
            } else if (cdr3Seen < static_cast<int>(cdr3Labels.size())) {
                label = cdr3Labels[static_cast<std::size_t>(cdr3Seen++)];
            }
        } else {
            const auto it = tbl->positions.find(r.position);
            if (it != tbl->positions.end()) label = it->second;
        }
        if (label.empty()) {
            out.warnings.push_back("IMGT " + positionLabel(r) +
                                   " has no published correspondence in the " + id +
                                   " table and is reported at its IMGT position");
            out.residues.push_back(nr);
            continue;
        }
        // Split "95A" into 95 + "A"; an insertion code is part of the identity.
        std::size_t cut = 0;
        while (cut < label.size() && std::isdigit(static_cast<unsigned char>(label[cut]))) ++cut;
        nr.position = std::stoi(label.substr(0, cut));
        nr.insertionCode = label.substr(cut);
        out.residues.push_back(nr);
    }
    if (id == "chothia")
        out.warnings.push_back(
            "Chothia here overrides Kabat for IMGT 1-104 only; the CDR3 and FR4 labels are the "
            "Kabat 1991 length rule, because the IMGT correspondence chart does not tabulate a "
            "Chothia insertion behaviour for them.");
    return out;
}

std::vector<double> relativeSasaFor(const AbDomain& domain, const Structure& s,
                                    const std::string& chainId) {
    std::vector<double> out;
    const Model* m = s.model(1);
    if (!m) return out;
    const Chain* chain = nullptr;
    for (const auto& c : m->chains)
        if (c.id == chainId) chain = &c;
    if (!chain) return out;

    const std::vector<char> seqv = sequenceOf(*chain);
    const std::string seq(seqv.begin(), seqv.end());
    const std::size_t offset = seq.find(domain.sequence.substr(0, 20));
    if (offset == std::string::npos) return out;

    const SasaResult sr = sasa(s);
    // Polymer residues of this chain, in order, paired with their relative SASA.
    std::vector<double> rel;
    for (const auto& r : chain->residues) {
        if (r.oneLetter() == 'X') continue;
        double v = -1.0;
        for (const auto& pr : sr.perResidue) {
            if (pr.key.chainId != chain->id || pr.key.authSeqId != r.authSeqId ||
                pr.key.insertionCode != r.insertionCode)
                continue;
            if (pr.relative.provenance == Provenance::Measured) v = pr.relative.value;
            break;
        }
        rel.push_back(v);
    }
    out.assign(domain.residues.size(), -1.0);
    for (std::size_t i = 0; i < domain.residues.size(); ++i) {
        const std::size_t k = offset + static_cast<std::size_t>(domain.residues[i].sequenceIndex);
        if (k < rel.size()) out[i] = rel[k];
    }
    return out;
}

}  // namespace biocad::bio
