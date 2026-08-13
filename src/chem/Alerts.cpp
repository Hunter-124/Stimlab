#include "chem/Alerts.h"

#include <fstream>

#include <nlohmann/json.hpp>

#include "chem/Aromaticity.h"
#include "chem/Rings.h"
#include "chem/Smarts.h"
#include "core/Assets.h"

namespace biocad::chem {
namespace {

// A compiled rule: the parsed pattern plus the text the UI renders. The pattern
// is parsed once at load, not per screen: parsing is the expensive half and the
// pack does not change while the app runs.
struct Rule {
    std::string   key;
    std::string   label;
    std::string   mechanism;
    std::string   citation;
    bool          warn = false;
    SmartsPattern pattern;
};

struct RulePack {
    std::vector<Rule>        rules;
    std::vector<std::string> errors;
};

// Required-string reader. A rule missing any of its text fields is dropped with a
// named error: an alert without a citation is exactly the uncitable claim this
// project removes, so it must not survive into a hit.
bool readField(const nlohmann::json& j, const char* field, std::string& out,
               const std::string& where, std::vector<std::string>& errors) {
    const auto it = j.find(field);
    if (it == j.end() || !it->is_string() || it->get<std::string>().empty()) {
        errors.push_back(where + ": missing or empty \"" + field + "\"");
        return false;
    }
    out = it->get<std::string>();
    return true;
}

RulePack loadRulePack() {
    RulePack pack;
    const std::filesystem::path dir = core::assetDir("packs");
    if (dir.empty()) {
        // No asset tree at all: say so, rather than reporting an empty screen that
        // reads as "this compound raised no alerts".
        pack.errors.push_back("no asset tree found, so assets/packs/rules/"
                              "alerts-bioactivation.json could not be read - no alert "
                              "screening is performed");
        return pack;
    }
    const std::filesystem::path file = dir / "rules" / "alerts-bioactivation.json";

    std::ifstream in(file);
    if (!in) {
        pack.errors.push_back("structural-alert pack not found: " + file.string() +
                              " - no alert screening is performed");
        return pack;
    }

    nlohmann::json doc;
    try {
        in >> doc;
    } catch (const std::exception& e) {
        pack.errors.push_back(file.string() + ": " + e.what());
        return pack;
    }

    const int version = doc.value("schemaVersion", 0);
    if (version != 1) {
        // An unknown version is a load error, not a silent skip: a pack the user
        // just edited that quietly does nothing is indistinguishable from a bug.
        pack.errors.push_back(file.string() + ": unsupported schemaVersion " +
                              std::to_string(version) + " (this build understands 1)");
        return pack;
    }

    const auto alerts = doc.find("alerts");
    if (alerts == doc.end() || !alerts->is_array()) {
        pack.errors.push_back(file.string() + ": no \"alerts\" array");
        return pack;
    }

    for (std::size_t i = 0; i < alerts->size(); ++i) {
        const nlohmann::json& a = (*alerts)[i];
        const std::string where = file.filename().string() + " alert #" + std::to_string(i + 1);
        Rule r;
        std::string smarts, severity;
        if (!readField(a, "key", r.key, where, pack.errors)) continue;
        if (!readField(a, "smarts", smarts, where + " (" + r.key + ")", pack.errors)) continue;
        if (!readField(a, "label", r.label, where + " (" + r.key + ")", pack.errors)) continue;
        if (!readField(a, "mechanism", r.mechanism, where + " (" + r.key + ")", pack.errors)) continue;
        if (!readField(a, "citation", r.citation, where + " (" + r.key + ")", pack.errors)) continue;
        if (!readField(a, "severity", severity, where + " (" + r.key + ")", pack.errors)) continue;

        // Only two severities exist by design. There is no "danger" level: a
        // substructure match cannot support a toxicity verdict, so the vocabulary
        // does not offer one to author.
        if (severity == "warn") {
            r.warn = true;
        } else if (severity == "info") {
            r.warn = false;
        } else {
            pack.errors.push_back(where + " (" + r.key + "): severity \"" + severity +
                                  "\" is not \"info\" or \"warn\"");
            continue;
        }

        std::string perr;
        auto pat = parseSmarts(smarts, &perr);
        if (!pat) {
            pack.errors.push_back(where + " (" + r.key + "): bad SMARTS \"" + smarts + "\" - " +
                                  perr);
            continue;
        }
        r.pattern = std::move(*pat);
        pack.rules.push_back(std::move(r));
    }
    return pack;
}

const RulePack& rulePack() {
    static const RulePack pack = loadRulePack();
    return pack;
}

}  // namespace

std::vector<AlertHit> screenAlerts(const Molecule& mol) {
    std::vector<AlertHit> hits;
    const RulePack& pack = rulePack();
    if (pack.rules.empty() || mol.atoms.empty()) return hits;

    // Perceived once for the whole sweep rather than once per rule: ring
    // perception is the same work for every pattern.
    const PreparedMolecule prepared = prepareMolecule(mol);

    for (const Rule& r : pack.rules) {
        for (const Match& m : matchAll(r.pattern, prepared.mol, prepared.rings)) {
            AlertHit h;
            h.key = r.key;
            h.label = r.label;
            h.mechanism = r.mechanism;
            h.citation = r.citation;
            h.warn = r.warn;
            h.atoms = m.atoms;
            hits.push_back(std::move(h));
        }
    }
    return hits;
}

const std::vector<std::string>& alertPackErrors() { return rulePack().errors; }

std::size_t alertRuleCount() { return rulePack().rules.size(); }

}  // namespace biocad::chem
