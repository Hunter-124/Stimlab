// agent/MockProvider.h - deterministic, offline ILlmProvider.
//
// Used by tests (a reproducible multi-tool-turn) and as the default assistant
// when no API key is configured. It inspects the conversation: on a user turn it
// routes a keyword to a panel and (if a highlight_panel tool is offered) emits a
// tool call; after the tool result comes back it emits a short closing message.
// Synthesis-style asks are refused - an offline mirror of the system prompt's
// safety boundary - so the safety behavior is testable without a network/key.
#pragma once

#include "contracts/ILlmProvider.h"

namespace stimlab::agent {

class MockProvider : public ILlmProvider {
public:
    std::string id() const override { return "mock"; }
    std::string displayName() const override { return "Offline assistant (no API key)"; }
    bool        ready() const override { return true; }  // always available

    LlmResponse send(const LlmRequest& req, const StreamCallback& onText) const override;
};

}  // namespace stimlab::agent
