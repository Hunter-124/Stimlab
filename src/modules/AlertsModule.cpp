#include "modules/AlertsModule.h"

#include "chem/Alerts.h"
#include "chem/Perceive.h"
#include "chem/Smiles.h"

namespace biocad {

AlertReport RealAlerts::screen(const Molecule& m) const {
    AlertReport r;
    r.moleculeId = m.id.empty() ? m.name : m.id;

    const auto parsed = chem::parsePerceived(m.smiles);
    if (!parsed) {
        // A structure that could not be parsed was not screened. Reporting "no
        // flags" here would be the one genuinely dishonest outcome available.
        r.summary = "The SMILES could not be parsed, so no substructure was screened. This is "
                    "not a finding about the compound.";
        return r;
    }

    if (chem::alertRuleCount() == 0) {
        std::string why = "The structural-alert pack could not be loaded, so nothing was "
                          "screened.";
        for (const auto& e : chem::alertPackErrors()) why += " " + e + ".";
        r.summary = why;
        return r;
    }

    for (const auto& h : chem::screenAlerts(*parsed)) {
        AlertFlag f;
        f.label = h.label;
        f.mechanism = h.mechanism;
        f.citation = h.citation;
        f.atomCount = static_cast<int>(h.atoms.size());
        // Two levels, by construction. There is no branch that can produce
        // Verdict::Danger or Verdict::Good from an alert.
        f.severity = h.warn ? Verdict::Warn : Verdict::Info;
        r.flags.push_back(std::move(f));
    }

    if (r.flags.empty()) {
        r.summary = "No alert in the built-in bioactivation pack matched (" +
                    std::to_string(chem::alertRuleCount()) +
                    " alerts screened). That is not a safety claim: the pack is a short, "
                    "in-house list of literature-associated motifs, and absence of a listed "
                    "motif says nothing about the motifs it does not list.";
        return r;
    }

    int warns = 0;
    for (const auto& f : r.flags) {
        if (f.severity == Verdict::Warn) ++warns;
    }
    r.summary = std::to_string(r.flags.size()) + " liability flag(s) raised (" +
                std::to_string(warns) + " warn, " + std::to_string(r.flags.size() - warns) +
                " info) out of " + std::to_string(chem::alertRuleCount()) +
                " alerts screened. Each flag means the matched substructure has been ASSOCIATED "
                "with reactive-metabolite formation in the literature; none of them is a "
                "statement that this compound is toxic, and many marketed drugs match several.";
    return r;
}

}  // namespace biocad
