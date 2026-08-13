// agent/OfflineAssistant.h - the deterministic, rule-based assistant.
//
// This is a shipping feature, not a test double: it is the assistant the app
// uses when no API key is configured, and it needs no network and no key. It
// inspects the conversation, routes a keyword in the user's turn to a panel and
// (if a highlight_panel tool is offered) emits that tool call; after the tool
// result comes back it emits a short closing message. Synthesis-style asks are
// refused - an offline mirror of the system prompt's safety boundary - so that
// boundary holds even with no provider reachable.
//
// What it deliberately does NOT do, because a limited feature described honestly
// is worth more than a grand name: it does not call a language model; it does not
// generate prose beyond its fixed templates; and it never answers a chemistry
// question it cannot route to a panel - the routing IS the answer.
#pragma once

#include "contracts/ILlmProvider.h"

namespace biocad::agent {

class OfflineAssistant : public ILlmProvider {
public:
    std::string id() const override { return "offline"; }
    std::string displayName() const override { return "Offline assistant (no API key)"; }
    bool        ready() const override { return true; }  // always available

    LlmResponse send(const LlmRequest& req, const StreamCallback& onText) const override;
};

}  // namespace biocad::agent
