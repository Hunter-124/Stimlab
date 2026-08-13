// agent/WebTools.h - keyless web search + fetch for the assistant (WP-K).
//
// web_search uses DuckDuckGo's HTML endpoint (no API key); web_fetch GETs a URL and
// extracts readable text. Results are cached + rate-limited under
// %APPDATA%/BioCAD/cache/web. All networking is compiled only with the science
// feature (libcurl); without it the tools report unavailable and the agent degrades.
// The HTML parsers (parseDuckDuckGoHtml / htmlToText) are PURE and always compiled so
// they are unit-testable from fixtures without any network.
//
// SAFETY SCOPE: these fetch public web text the assistant reasons over (pharmacology,
// literature, legal status). They return third-party content verbatim; the agent's
// system prompt keeps it from turning that into synthesis guidance.
#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace biocad::agent {

// True iff this build has networking (science feature / libcurl) compiled in.
bool webToolsAvailable();

// True iff this build can render JavaScript pages (BIOCAD_ENABLE_WEBVIEW2) AND the
// machine's Evergreen WebView2 Runtime is installed. False -> web_fetch uses curl only.
bool webViewRenderAvailable();

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
// When renderJs is true and this build has WebView2 (and the Evergreen Runtime is
// present), the page is loaded in a headless WebView2 so its JavaScript runs before
// extraction; on any failure it transparently falls back to the plain curl GET.
WebFetchResult webFetch(const std::string& url, std::size_t maxChars = 8000, bool renderJs = false);

// Load `url` in a HEADLESS WebView2 (full Chromium), let its JavaScript run, and return
// the rendered DOM HTML (document.documentElement.outerHTML). nullopt if this build has
// no WebView2, the Evergreen Runtime is absent, navigation fails, or it times out - in
// which case the caller falls back to the curl path. Compiled in only with
// BIOCAD_ENABLE_WEBVIEW2; otherwise a stub returns nullopt.
std::optional<std::string> webFetchRenderedHtml(const std::string& url, int timeoutMs = 20000);

// ---- pure, network-free helpers (unit-testable from fixtures) ----------------
// Parse a DuckDuckGo HTML results page into structured hits (real URLs un-wrapped
// from DDG's /l/?uddg= redirect, titles/snippets HTML-stripped).
std::vector<WebHit> parseDuckDuckGoHtml(const std::string& html, int maxResults = 6);

// Strip <script>/<style>, remove tags, decode common entities, collapse whitespace,
// and truncate to maxChars.
std::string htmlToText(const std::string& html, std::size_t maxChars = 8000);

}  // namespace biocad::agent
