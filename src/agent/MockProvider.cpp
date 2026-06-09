#include "agent/MockProvider.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <utility>

namespace stimlab::agent {
namespace {

std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

bool contains(const std::string& hay, const char* needle) {
    return hay.find(needle) != std::string::npos;
}

// Lowercased text of the most recent user message (empty if none).
std::string lastUserText(const LlmRequest& req) {
    for (auto it = req.messages.rbegin(); it != req.messages.rend(); ++it)
        if (it->role == ChatRole::User && !it->text.empty()) return toLower(it->text);
    return {};
}

bool toolOffered(const LlmRequest& req, const char* name) {
    for (const auto& t : req.tools)
        if (t.name == name) return true;
    return false;
}

// Offline mirror of the system prompt's hard boundary.
bool mentionsSynthesis(const std::string& q) {
    return contains(q, "synthes") || contains(q, "manufactur") || contains(q, "precursor") ||
           contains(q, "how do i make") || contains(q, "how to make") || contains(q, "cook ") ||
           contains(q, "reaction route") || contains(q, "route to make");
}

// Map a keyword to {panelId, one-line blurb}.
std::pair<std::string, std::string> routePanel(const std::string& q) {
    if (contains(q, "dock") || contains(q, "target") || contains(q, "binding") ||
        contains(q, "affinity") || contains(q, "dat") || contains(q, "sert"))
        return {"Docking",
                "the Docking panel - the target dropdown selects DAT/NET/SERT/TAAR1 and the score is "
                "predicted binding affinity (pharmacology, not a make-it signal)."};
    if (contains(q, "absor") || contains(q, "bioavail") || contains(q, "permeab") ||
        contains(q, "bbb") || contains(q, "pk") || contains(q, "pharmacokinet"))
        return {"Absorption",
                "the Absorption / PK panel - HIA, oral F%, Caco-2 permeability, BBB partition and "
                "P-gp efflux."};
    if (contains(q, "stab") || contains(q, "shelf") || contains(q, "degrad"))
        return {"Stability",
                "the Stability panel - resistance to hydrolysis/oxidation/photolysis/thermal/pH and a "
                "shelf-life estimate (this replaces any manufacturability score - out of scope)."};
    if (contains(q, "metabol") || contains(q, "admet") || contains(q, "cyp") ||
        contains(q, "herg") || contains(q, "toxic"))
        return {"Metabolism",
                "the Metabolism (ADMET) panel - metabolic routes, harmful metabolites, drug-drug "
                "interactions, hERG and safety flags."};
    if (contains(q, "similar") || contains(q, "tanimoto") || contains(q, "fingerprint"))
        return {"Similarity",
                "the Similarity panel - structural + pharmacophore similarity to the curated "
                "known-substance reference set."};
    if (contains(q, "legal") || contains(q, "analog") || contains(q, "schedule") ||
        contains(q, "controlled"))
        return {"Legal",
                "the Legal Analog panel - a substantial-similarity scorecard vs controlled references "
                "(illustrative, not legal advice)."};
    if (contains(q, "compare") || contains(q, "side by side") || contains(q, "side-by-side"))
        return {"Compare",
                "the Compare panel - up to three compounds charted across stability, absorption and "
                "ADMET."};
    if (contains(q, "smiles") || contains(q, "paste") || contains(q, "enter ") ||
        contains(q, "input"))
        return {"Input",
                "the Molecule Input panel - paste any SMILES to run the full analysis on an arbitrary "
                "structure."};
    if (contains(q, "deriv") || contains(q, "candidate") || contains(q, "tweak") ||
        contains(q, "modif"))
        return {"Analog",
                "the Analog Explorer - model a candidate derivative and see its nearest existing "
                "sample, legal-analog score and predicted byproducts (analysis only)."};
    if (contains(q, "structure") || contains(q, "3d") || contains(q, "viewer") ||
        contains(q, "propert"))
        return {"Structure",
                "the Structure Workbench - identity, physicochemical properties and a live 3D viewer."};
    if (contains(q, "librar") || contains(q, "browse") || contains(q, "pick") ||
        contains(q, "select"))
        return {"Library", "the Library - browse and select compounds from the default + imported set."};
    if (contains(q, "run") || contains(q, "history"))
        return {"Runs", "the Runs panel - history of analyses with status and summaries."};
    if (contains(q, "setting") || contains(q, "provider") || contains(q, "key") ||
        contains(q, "api"))
        return {"Settings", "the Settings panel - AI provider/key, compute mode and storage paths."};
    return {"Dashboard",
            "the Dashboard - library size, recent runs and a snapshot of the selected compound. "
            "StimLab predicts what a compound IS and DOES; it does not provide synthesis guidance."};
}

void streamOut(const StreamCallback& onText, const std::string& text) {
    if (!onText) return;
    // Emit in a few chunks so the streaming path is exercised like the real provider.
    const std::size_t chunks = 3;
    const std::size_t step = (text.size() + chunks - 1) / std::max<std::size_t>(1, chunks);
    for (std::size_t i = 0; i < text.size(); i += step) onText(text.substr(i, step));
}

}  // namespace

LlmResponse MockProvider::send(const LlmRequest& req, const StreamCallback& onText) const {
    LlmResponse r;
    r.message.role = ChatRole::Assistant;

    const ChatMessage* last = req.messages.empty() ? nullptr : &req.messages.back();

    // Second leg of a tool turn: we just got tool results -> close out.
    if (last && last->role == ChatRole::Tool) {
        std::string txt = "Done - I've taken you there and highlighted it. Anything else?";
        streamOut(onText, txt);
        r.message.text = txt;
        r.stop = StopReason::EndTurn;
        return r;
    }

    const std::string q = lastUserText(req);

    if (mentionsSynthesis(q)) {
        std::string txt =
            "I can't help with making, synthesizing, or manufacturing substances - that is outside "
            "what StimLab does. I can analyze what a compound IS and DOES: structure and properties, "
            "molecular stability, absorption/PK, ADMET/metabolism, target binding (docking), "
            "similarity to known substances, and legal-analog scoring.";
        streamOut(onText, txt);
        r.message.text = txt;
        r.stop = StopReason::Refusal;
        return r;
    }

    auto [panel, blurb] = routePanel(q);

    if (toolOffered(req, "highlight_panel")) {
        std::string say = "Sure - I'll take you to " + blurb;
        streamOut(onText, say);
        r.message.text = say;
        ToolCall c;
        c.id = "mock_call_1";
        c.name = "highlight_panel";
        c.arguments = {{"panel", panel}, {"explanation", blurb}};
        r.message.toolCalls.push_back(std::move(c));
        r.stop = StopReason::ToolUse;
        return r;
    }

    // No navigation tool available: answer in text only.
    std::string say = "You'll want " + blurb;
    streamOut(onText, say);
    r.message.text = say;
    r.stop = StopReason::EndTurn;
    return r;
}

}  // namespace stimlab::agent
