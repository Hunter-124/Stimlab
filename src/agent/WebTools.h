// agent/WebTools.h - keyless web search + fetch for the assistant (WP-K).
//
// web_search uses DuckDuckGo's HTML endpoint (no API key); web_fetch GETs a URL and
// extracts readable text. Results are cached + rate-limited under
// %APPDATA%/StimLab/cache/web. All networking is compiled only with the science
// feature (libcurl); without it the tools report unavailable and the agent degrades.
// The HTML parsers (parseDuckDuckGoHtml / htmlToText) are PURE and always compiled so
// they are unit-testable from fixtures without any network.
//
// SAFETY SCOPE: these fetch public web text the assistant reasons over (pharmacology,
// literature, legal status). They return third-party content verbatim; the agent's
// system prompt keeps it from turning that into synthesis guidance.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace stimlab::agent {

// True iff this build has networking (science feature / libcurl) compiled in.
bool webToolsAvailable();

struct WebHit {
    std::string title;
    std::string url;
    std::string snippet;
};

struct WebSearchResult {
    bool                ok = false;
    std::vector<WebHit> hits;
    std::string         error;
    bool                fromCache = false;
};

// Keyless search via DuckDuckGo HTML. Cached + rate-limited. ok=false (with an error)
// when networking is absent or the request fails.
WebSearchResult webSearch(const std::string& query, int maxResults = 6);

struct WebFetchResult {
    bool        ok = false;
    std::string text;      // HTML-stripped readable text, truncated
    std::string title;
    std::string finalUrl;
    std::string error;
    bool        fromCache = false;
};

// Fetch a URL and extract readable text (scripts/styles/tags stripped). Cached.
WebFetchResult webFetch(const std::string& url, std::size_t maxChars = 8000);

// ---- pure, network-free helpers (unit-testable from fixtures) ----------------
// Parse a DuckDuckGo HTML results page into structured hits (real URLs un-wrapped
// from DDG's /l/?uddg= redirect, titles/snippets HTML-stripped).
std::vector<WebHit> parseDuckDuckGoHtml(const std::string& html, int maxResults = 6);

// Strip <script>/<style>, remove tags, decode common entities, collapse whitespace,
// and truncate to maxChars.
std::string htmlToText(const std::string& html, std::size_t maxChars = 8000);

}  // namespace stimlab::agent
