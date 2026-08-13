#include "chem/Analysis.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include <nlohmann/json.hpp>

#include "chem/Smarts.h"
#include "core/Assets.h"

namespace biocad::chem {

namespace {

std::uint64_t mix(std::uint64_t h, std::uint64_t v) {
    h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return h;
}

// The pack's `key` names a member of FunctionalGroups, so a rule can only ever
// set the flag it claims to define and an unknown key is an error instead of a
// silently ignored line.
using Flag = bool FunctionalGroups::*;

const std::unordered_map<std::string, Flag>& flagsByKey() {
    static const std::unordered_map<std::string, Flag> kFlags = {
        {"aromaticRing", &FunctionalGroups::aromaticRing},
        {"phenol", &FunctionalGroups::phenol},
        {"catechol", &FunctionalGroups::catechol},
        {"ester", &FunctionalGroups::ester},
        {"carboxylicAcid", &FunctionalGroups::carboxylicAcid},
        {"amide", &FunctionalGroups::amide},
        {"ketone", &FunctionalGroups::ketone},
        {"arylKetone", &FunctionalGroups::arylKetone},
        {"aldehyde", &FunctionalGroups::aldehyde},
        {"ether", &FunctionalGroups::ether},
        {"methylenedioxy", &FunctionalGroups::methylenedioxy},
        {"primaryAmine", &FunctionalGroups::primaryAmine},
        {"secondaryAmine", &FunctionalGroups::secondaryAmine},
        {"tertiaryAmine", &FunctionalGroups::tertiaryAmine},
        {"basicAmine", &FunctionalGroups::basicAmine},
        {"nitrile", &FunctionalGroups::nitrile},
        {"nitro", &FunctionalGroups::nitro},
        {"halogen", &FunctionalGroups::halogen},
        {"sulfoxide", &FunctionalGroups::sulfoxide},
        {"sulfone", &FunctionalGroups::sulfone},
        {"phenethylamine", &FunctionalGroups::phenethylamine},
        {"catecholamine", &FunctionalGroups::catecholamine},
        {"anilide", &FunctionalGroups::anilide},
        {"maoLabileAmine", &FunctionalGroups::maoLabileAmine},
    };
    return kFlags;
}

struct GroupRule {
    std::string   key;
    SmartsPattern pattern;
    Flag          flag = nullptr;
};

struct GroupRules {
    std::vector<GroupRule>   rules;
    std::vector<std::string> errors;
};

GroupRules loadGroupRules() {
    GroupRules out;
    const auto packs = core::assetDir("packs");
    if (packs.empty()) {
        out.errors.push_back(
            "functional-groups: no asset tree found - expected "
            "assets/packs/rules/functional-groups.json beside the executable. No functional "
            "group can be perceived.");
        return out;
    }
    const auto path = packs / "rules" / "functional-groups.json";
    std::ifstream in(path);
    if (!in) {
        out.errors.push_back("functional-groups: rule pack not found at " + path.string() +
                             " - no functional group can be perceived");
        return out;
    }

    nlohmann::json j;
    try {
        in >> j;
    } catch (const std::exception& e) {
        out.errors.push_back("functional-groups: " + path.string() + " is not valid JSON: " +
                             e.what());
        return out;
    }

    const int version = j.value("schemaVersion", 0);
    if (version != 1) {
        out.errors.push_back("functional-groups: unsupported schemaVersion " +
                             std::to_string(version) + " in " + path.string() +
                             " (this build understands 1)");
        return out;
    }
    const auto groups = j.find("groups");
    if (groups == j.end() || !groups->is_array()) {
        out.errors.push_back("functional-groups: " + path.string() +
                             " has no \"groups\" array");
        return out;
    }

    for (const auto& g : *groups) {
        const std::string key = g.value("key", std::string{});
        const std::string smarts = g.value("smarts", std::string{});
        if (key.empty()) {
            out.errors.push_back("functional-groups: a rule has no \"key\"");
            continue;
        }
        const auto flag = flagsByKey().find(key);
        if (flag == flagsByKey().end()) {
            out.errors.push_back("functional-groups: unknown group key \"" + key +
                                 "\" (no such flag in chem::FunctionalGroups)");
            continue;
        }
        if (std::any_of(out.rules.begin(), out.rules.end(),
                        [&](const GroupRule& r) { return r.key == key; })) {
            out.errors.push_back("functional-groups: duplicate group key \"" + key +
                                 "\" - which definition wins would be arbitrary");
            continue;
        }
        if (smarts.empty()) {
            out.errors.push_back("functional-groups: group \"" + key + "\" has no \"smarts\"");
            continue;
        }
        std::string error;
        auto pattern = parseSmarts(smarts, &error);
        if (!pattern) {
            out.errors.push_back("functional-groups: group \"" + key + "\" has malformed SMARTS \"" +
                                 smarts + "\": " + error);
            continue;
        }
        out.rules.push_back(GroupRule{key, std::move(*pattern), flag->second});
    }

    for (const auto& [key, flag] : flagsByKey()) {
        (void)flag;
        if (std::none_of(out.rules.begin(), out.rules.end(),
                         [&](const GroupRule& r) { return r.key == key; })) {
            out.errors.push_back("functional-groups: no rule defines \"" + key +
                                 "\" - that flag is always false in this run");
        }
    }
    return out;
}

// Parsed once per process: parsing 24 SMARTS patterns for every molecule would
// dominate the cost of perception, and the pack cannot change under a live run.
const GroupRules& groupRules() {
    static const GroupRules kRules = loadGroupRules();
    return kRules;
}

}  // namespace

const std::vector<std::string>& groupPackErrors() { return groupRules().errors; }

FunctionalGroups detectGroups(const Molecule& m) {
    FunctionalGroups g;
    if (m.atoms.empty()) return g;

    // Rings and aromaticity once for the whole rule sweep, and on a copy: callers
    // hand us whatever parseSmiles produced, aromatic patterns need perception,
    // and perception must not be a side effect of asking a question.
    const PreparedMolecule prepared = prepareMolecule(m);
    for (const GroupRule& rule : groupRules().rules) {
        if (matches(rule.pattern, prepared.mol, prepared.rings)) g.*(rule.flag) = true;
    }
    return g;
}

std::vector<std::uint32_t> morganFingerprint(const Molecule& m, int radius) {
    const int n = static_cast<int>(m.atoms.size());
    std::vector<std::uint64_t> inv(n);
    for (int i = 0; i < n; ++i) {
        const Atom& a = m.atoms[i];
        std::uint64_t h = 1469598103934665603ULL;
        h = mix(h, static_cast<std::uint64_t>(a.z));
        h = mix(h, static_cast<std::uint64_t>(a.degree()));
        h = mix(h, static_cast<std::uint64_t>(a.totalH()));
        h = mix(h, static_cast<std::uint64_t>(a.charge + 8));
        h = mix(h, static_cast<std::uint64_t>(a.aromatic ? 1 : 0));
        h = mix(h, static_cast<std::uint64_t>(a.inRing ? 1 : 0));
        inv[i] = h;
    }

    std::vector<std::uint32_t> bits;
    auto emit = [&](std::uint64_t v) { bits.push_back(static_cast<std::uint32_t>(v ^ (v >> 32))); };
    for (int i = 0; i < n; ++i) emit(inv[i]);

    for (int r = 0; r < radius; ++r) {
        std::vector<std::uint64_t> next(n);
        for (int i = 0; i < n; ++i) {
            std::vector<std::uint64_t> env;
            for (int bi : m.atoms[i].bonds) {
                const Bond& b = m.bonds[bi];
                env.push_back(mix(static_cast<std::uint64_t>(b.order * 2), inv[b.other(i)]));
            }
            std::sort(env.begin(), env.end());
            std::uint64_t h = mix(inv[i], static_cast<std::uint64_t>(r + 1));
            for (std::uint64_t e : env) h = mix(h, e);
            next[i] = h;
        }
        inv.swap(next);
        for (int i = 0; i < n; ++i) emit(inv[i]);
    }

    std::sort(bits.begin(), bits.end());
    bits.erase(std::unique(bits.begin(), bits.end()), bits.end());
    return bits;
}

double tanimoto(const std::vector<std::uint32_t>& a, const std::vector<std::uint32_t>& b) {
    if (a.empty() && b.empty()) return 1.0;
    int inter = 0;
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] == b[j]) { ++inter; ++i; ++j; }
        else if (a[i] < b[j]) ++i;
        else ++j;
    }
    const int uni = static_cast<int>(a.size() + b.size()) - inter;
    return uni > 0 ? static_cast<double>(inter) / uni : 0.0;
}

}  // namespace biocad::chem
