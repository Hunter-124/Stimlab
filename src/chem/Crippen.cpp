// chem/Crippen.cpp - see chem/Crippen.h for the method, the citation, and the
// hydrogen mapping. This file contains no chemistry constants: every class,
// SMARTS and contribution is read from assets/packs/descriptors/crippen.json.
#include "chem/Crippen.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "chem/Smarts.h"
#include "core/Assets.h"

namespace biocad::chem {
namespace {

constexpr int kSchemaVersion = 1;

struct CrippenClass {
    std::string   type;
    std::string   smarts;
    double        logP = 0.0;
    double        mr = 0.0;
    SmartsPattern pattern;
};

// Contributions for the four hydrogen classes, kept beside the heavy-atom table
// because the pack cannot express them (chem::Molecule has no explicit H atoms).
// Values are the same published table, transcribed once.
struct HydrogenClass {
    const char* type;
    double      logP;
    double      mr;
};
constexpr HydrogenClass kH1{"H1", 0.1230, 1.057};
constexpr HydrogenClass kH2{"H2", -0.2677, 1.395};
constexpr HydrogenClass kH3{"H3", 0.2142, 0.9627};
constexpr HydrogenClass kH4{"H4", 0.2980, 1.805};

// biocad_chem cannot link biocad_packs (biocad_contracts already depends on
// biocad_chem, so that direction would be a cycle), so the descriptor pack is
// located through core::assetRoot(), the resolver both layers share. The
// environment override exists so a test or a harness can point at the in-tree
// copy without a packaged layout.
std::filesystem::path findPack() {
    std::error_code ec;
    if (const char* env = std::getenv("BIOCAD_DESCRIPTOR_DIR"); env && *env) {
        const auto p = std::filesystem::path(env) / "crippen.json";
        return std::filesystem::is_regular_file(p, ec) ? p : std::filesystem::path{};
    }
    const auto dir = core::assetDir("packs/descriptors");
    if (!dir.empty()) {
        const auto p = dir / "crippen.json";
        if (std::filesystem::is_regular_file(p, ec)) return p;
    }
    return {};
}

struct Table {
    bool                      ok = false;
    std::string               note;
    std::vector<CrippenClass> classes;
};

Table loadTable() {
    Table t;
    const auto path = findPack();
    if (path.empty()) {
        t.note = "assets/packs/descriptors/crippen.json not found (set BIOCAD_DESCRIPTOR_DIR "
                 "to its directory)";
        return t;
    }
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        t.note = path.string() + ": cannot open";
        return t;
    }
    std::ostringstream buf;
    buf << in.rdbuf();
    try {
        const auto j = nlohmann::json::parse(buf.str());
        const int version = j.at("schemaVersion").get<int>();
        if (version != kSchemaVersion) {
            t.note = path.string() + ": schemaVersion " + std::to_string(version) +
                     " is not supported (expected " + std::to_string(kSchemaVersion) + ")";
            return t;
        }
        for (const auto& c : j.at("classes")) {
            CrippenClass cc;
            cc.type = c.at("type").get<std::string>();
            cc.smarts = c.at("smarts").get<std::string>();
            cc.logP = c.at("logP").get<double>();
            cc.mr = c.at("mr").get<double>();
            // The pack stores the paper's pattern, whose FIRST atom is the atom
            // being classified. Wrapping it in a recursive query turns it into a
            // single-atom pattern anchored on that atom, so every atom that can
            // start the pattern is reported exactly once - matchAll's
            // distinct-by-atom-set rule would otherwise collapse two orientations
            // of the same set into one match and leave an atom untyped.
            std::string err;
            auto pat = parseSmarts("[$(" + cc.smarts + ")]", &err);
            if (!pat) {
                t.note = path.string() + ": class " + cc.type + " has an unparseable SMARTS \"" +
                         cc.smarts + "\": " + err;
                return t;
            }
            cc.pattern = std::move(*pat);
            t.classes.push_back(std::move(cc));
        }
    } catch (const std::exception& e) {
        t.note = path.string() + ": " + e.what();
        return t;
    }
    if (t.classes.empty()) {
        t.note = path.string() + ": no classes";
        return t;
    }
    t.ok = true;
    t.note = path.string();
    return t;
}

// Cached on success only. A failed load is retried on the next call so a process
// that fixes its working directory or sets BIOCAD_DESCRIPTOR_DIR is not poisoned
// by one early miss.
const Table& table() {
    static std::mutex mu;
    static Table cached;
    std::lock_guard<std::mutex> lock(mu);
    if (!cached.ok) cached = loadTable();
    return cached;
}

bool hasDoubleBondTo(const Molecule& m, int atom, const std::vector<int>& elements) {
    for (int bi : m.atoms[static_cast<std::size_t>(atom)].bonds) {
        const Bond& b = m.bonds[static_cast<std::size_t>(bi)];
        if (b.order != 2.0) continue;
        const int z = m.atoms[static_cast<std::size_t>(b.other(atom))].z;
        for (int e : elements)
            if (z == e) return true;
    }
    return false;
}

// See the mapping table in chem/Crippen.h. First match wins, in the reference
// table's own order (H2 rules, then H3, then H4).
const HydrogenClass& hydrogenClass(const Molecule& m, int i) {
    const Atom& a = m.atoms[static_cast<std::size_t>(i)];
    if (a.z == 6) return kH1;   // [#1][#6]
    if (a.z == 7) return kH3;   // [#1][#7]
    if (a.z != 8) return kH2;   // [#1][!C;!N;!O] - S-H, P-H, halide H, ...
    // Oxygen: look one atom further, as the reference patterns do.
    for (int nb : a.nbr) {
        const Atom& n = m.atoms[static_cast<std::size_t>(nb)];
        if (n.z == 6) {
            // [#1]O[CX4,c] is tested before the acid rule, so an sp3 or aromatic
            // carbon neighbour means H2 (alcohol, phenol).
            if (n.aromatic || (n.degree() + n.totalH()) == 4) return kH2;
            if (hasDoubleBondTo(m, nb, {6, 7, 8, 16})) return kH4;  // [#1]OC=[#6,#7,O,S]
            return kH2;
        }
        if (n.z == 7) return kH3;                    // [#1]O[#7]
        if (n.z == 8 || n.z == 16) return kH4;       // [#1]O[O,S]
        return kH2;                                  // [#1]O[!C;!N;!O;!S]
    }
    return kH2;  // water, [OH2]
}

}  // namespace

const char* crippenCitation() {
    return "Wildman-Crippen (J Chem Inf Comput Sci 1999;39:868-873), method RMS ~0.67 log units";
}

CrippenResult crippen(const Molecule& m) {
    CrippenResult out;
    const Table& t = table();
    if (!t.ok) {
        out.note = "Wildman-Crippen parameters unavailable: " + t.note;
        return out;
    }

    // Ring and aromaticity perception is ours, not the caller's: the same
    // molecule written Kekule or aromatic must give the same number.
    const PreparedMolecule prep = prepareMolecule(m);
    const Molecule& mol = prep.mol;
    const std::size_t n = mol.atoms.size();

    std::vector<int> assigned(n, -1);  // index into t.classes
    for (std::size_t ci = 0; ci < t.classes.size(); ++ci) {
        const auto hits = matchAll(t.classes[ci].pattern, mol, prep.rings, n + 1);
        for (const auto& hit : hits) {
            if (hit.atoms.empty()) continue;
            int& slot = assigned[static_cast<std::size_t>(hit.atoms[0])];
            if (slot < 0) slot = static_cast<int>(ci);
        }
    }

    out.atomTypes.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        const int ci = assigned[i];
        if (ci < 0) {
            out.atomTypes[i] = "?";
            out.unclassified.push_back(std::to_string(i) + ":" +
                                       std::to_string(mol.atoms[i].z));
            continue;
        }
        const CrippenClass& cc = t.classes[static_cast<std::size_t>(ci)];
        out.atomTypes[i] = cc.type;
        out.logP += cc.logP;
        out.molarRefractivity += cc.mr;

        const int h = mol.atoms[i].totalH();
        if (h > 0) {
            const HydrogenClass& hc = hydrogenClass(mol, static_cast<int>(i));
            out.logP += h * hc.logP;
            out.molarRefractivity += h * hc.mr;
        }
    }

    if (!out.unclassified.empty()) {
        out.note = "unclassified heavy atoms (index:Z): ";
        for (std::size_t k = 0; k < out.unclassified.size(); ++k) {
            if (k) out.note += ", ";
            out.note += out.unclassified[k];
        }
        return out;  // ok stays false: a partial sum is not the method's answer
    }
    out.ok = true;
    out.note = t.note;
    return out;
}

}  // namespace biocad::chem
