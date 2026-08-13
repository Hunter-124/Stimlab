// agent/AnthropicProvider.h - real provider: Anthropic Messages API.
//
// Streams the Messages API over libcurl (SSE) and assembles text + tool_use
// blocks into an LlmResponse. ALL networking is compiled only when
// BIOCAD_HAVE_SCIENCE is defined (the `science` vcpkg feature provides curl,
// using Windows-native Schannel TLS). Without it the class still exists but
// ready() is false and send() returns a clear error, so biocad_ui links either
// way and the agent simply falls back to the offline MockProvider.
//
// SAFETY SCOPE: transport only - the boundary is the system prompt + tool set.
#pragma once

#include <mutex>
#include <string>

#include "contracts/ILlmProvider.h"

namespace biocad::agent {

class AnthropicProvider : public ILlmProvider {
public:
    AnthropicProvider() = default;

    void setApiKey(std::string key);          // decrypted key (never persisted here)
    [[nodiscard]] bool hasKey() const;

    std::string id() const override { return "anthropic"; }
    std::string displayName() const override { return "Anthropic (Claude)"; }
    bool        ready() const override;        // true iff built with curl AND key present
    LlmResponse send(const LlmRequest& req, const StreamCallback& onText) const override;

    // True iff this build was compiled with the libcurl transport (science feature).
    [[nodiscard]] static bool transportAvailable();

private:
    mutable std::mutex mu_;
    std::string apiKey_;
};

}  // namespace biocad::agent
