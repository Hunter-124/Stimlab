#include "agent/WebTools.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <sstream>
#include <thread>

#include <nlohmann/json.hpp>

#include "core/AppPaths.h"
#include "core/Hash.h"
#include "core/Version.h"

#ifdef BIOCAD_HAVE_SCIENCE
#  include <curl/curl.h>
#endif

namespace biocad::agent {
namespace {

namespace fs = std::filesystem;

std::string lower(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

// Erase every <open ...></close> block (case-insensitive) - scripts/styles.
void removeBlocks(std::string& s, const std::string& open, const std::string& close) {
    const std::string ls = lower(s), lo = lower(open), lc = lower(close);
    std::string out;
    out.reserve(s.size());
    size_t i = 0;
    while (i < s.size()) {
        const size_t a = ls.find(lo, i);
        if (a == std::string::npos) {
            out.append(s, i, std::string::npos);
            break;
        }
        out.append(s, i, a - i);
        const size_t b = ls.find(lc, a);
        if (b == std::string::npos) break;  // unterminated -> drop the rest
        i = b + lc.size();
    }
    s.swap(out);
}

std::string decodeEntities(std::string s) {
    struct E { const char* from; const char* to; };
    static const E kE[] = {{"&nbsp;", " "}, {"&lt;", "<"},   {"&gt;", ">"},  {"&quot;", "\""},
                           {"&#39;", "'"},  {"&apos;", "'"}, {"&#x27;", "'"}, {"&#x2F;", "/"},
                           {"&amp;", "&"}};  // &amp; LAST so we don't double-decode
    for (const auto& e : kE) {
        std::string out;
        const std::string from = e.from;
        size_t i = 0, p;
        while ((p = s.find(from, i)) != std::string::npos) {
            out.append(s, i, p - i);
            out += e.to;
            i = p + from.size();
        }
        out.append(s, i, std::string::npos);
        s.swap(out);
    }
    return s;
}

std::string collapseWs(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    bool inWs = false;
    for (char c : s) {
        const bool ws = (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f' || c == '\v');
        if (ws) {
            inWs = true;
        } else {
            if (inWs && !out.empty()) out.push_back(' ');
            inWs = false;
            out.push_back(c);
        }
    }
    return out;
}

std::string urlDecode(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '%' && i + 2 < s.size() && std::isxdigit((unsigned char)s[i + 1]) &&
            std::isxdigit((unsigned char)s[i + 2])) {
            auto hex = [](char c) {
                c = static_cast<char>(std::tolower((unsigned char)c));
                return c <= '9' ? c - '0' : 10 + (c - 'a');
            };
            out.push_back(static_cast<char>(hex(s[i + 1]) * 16 + hex(s[i + 2])));
            i += 2;
        } else if (s[i] == '+') {
            out.push_back(' ');
        } else {
            out.push_back(s[i]);
        }
    }
    return out;
}

// Read the value of attr `name="..."` from an HTML open tag.
std::string tagAttr(const std::string& tag, const std::string& name) {
    const std::string key = name + "=\"";
    const size_t a = tag.find(key);
    if (a == std::string::npos) return {};
    const size_t start = a + key.size();
    const size_t end = tag.find('"', start);
    if (end == std::string::npos) return {};
    return tag.substr(start, end - start);
}

// Turn a DDG result href into a real URL (un-wrap /l/?uddg=<encoded>&...).
std::string unwrapDdgHref(const std::string& href) {
    const size_t u = href.find("uddg=");
    if (u != std::string::npos) {
        const size_t start = u + 5;
        const size_t end = href.find('&', start);
        return urlDecode(href.substr(start, end == std::string::npos ? std::string::npos : end - start));
    }
    if (href.rfind("//", 0) == 0) return "https:" + href;
    return href;
}

#ifdef BIOCAD_HAVE_SCIENCE
std::string extractTitle(const std::string& html) {
    const std::string lh = lower(html);
    const size_t a = lh.find("<title");
    if (a == std::string::npos) return {};
    const size_t gt = html.find('>', a);
    const size_t end = lh.find("</title>", gt == std::string::npos ? a : gt);
    if (gt == std::string::npos || end == std::string::npos) return {};
    return htmlToText(html.substr(gt + 1, end - gt - 1), 200);
}

size_t writeCb(char* ptr, size_t size, size_t nmemb, void* ud) {
    static_cast<std::string*>(ud)->append(ptr, size * nmemb);
    return size * nmemb;
}

// Process-wide minimum spacing between outbound requests (politeness / rate limit).
void rateLimit() {
    static std::mutex mu;
    static std::chrono::steady_clock::time_point last{};
    std::lock_guard<std::mutex> lk(mu);
    const auto now = std::chrono::steady_clock::now();
    const auto minGap = std::chrono::milliseconds(700);
    if (last.time_since_epoch().count() != 0 && (now - last) < minGap)
        std::this_thread::sleep_for(minGap - (now - last));
    last = std::chrono::steady_clock::now();
}

bool httpGet(const std::string& url, std::string& body, std::string& err, long& code) {
    static std::once_flag once;
    std::call_once(once, [] { curl_global_init(CURL_GLOBAL_DEFAULT); });
    rateLimit();
    CURL* c = curl_easy_init();
    if (!c) {
        err = "curl init failed";
        return false;
    }
    curl_easy_setopt(c, CURLOPT_URL, url.c_str());
    curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, writeCb);
    curl_easy_setopt(c, CURLOPT_WRITEDATA, &body);
    curl_easy_setopt(c, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(c, CURLOPT_MAXREDIRS, 5L);
    curl_easy_setopt(c, CURLOPT_TIMEOUT, 25L);
    curl_easy_setopt(c, CURLOPT_CONNECTTIMEOUT, 15L);
    curl_easy_setopt(c, CURLOPT_NOPROGRESS, 1L);
    curl_easy_setopt(c, CURLOPT_ACCEPT_ENCODING, "");
    const std::string ua =
        std::string("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 BioCAD/") +
        kBioCadVersionShort;
    curl_easy_setopt(c, CURLOPT_USERAGENT, ua.c_str());
    const CURLcode rc = curl_easy_perform(c);
    curl_easy_getinfo(c, CURLINFO_RESPONSE_CODE, &code);
    curl_easy_cleanup(c);
    if (rc != CURLE_OK) {
        err = curl_easy_strerror(rc);
        return false;
    }
    return true;
}

std::string urlEscape(const std::string& s) {
    CURL* c = curl_easy_init();
    char* e = c ? curl_easy_escape(c, s.c_str(), static_cast<int>(s.size())) : nullptr;
    std::string out = e ? e : s;
    if (e) curl_free(e);
    if (c) curl_easy_cleanup(c);
    return out;
}

fs::path webCacheDir() {
    const fs::path d = AppPaths::instance().cache() / "web";
    std::error_code ec;
    fs::create_directories(d, ec);
    return d;
}

bool freshEnough(const fs::path& p, std::chrono::hours ttl) {
    std::error_code ec;
    if (!fs::exists(p, ec)) return false;
    const auto mt = fs::last_write_time(p, ec);
    if (ec) return false;
    return (fs::file_time_type::clock::now() - mt) < ttl;
}

std::string readFile(const fs::path& p) {
    std::ifstream in(p, std::ios::binary);
    if (!in) return {};
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

void writeFile(const fs::path& p, const std::string& s) {
    std::ofstream out(p, std::ios::binary);
    out << s;
}
#endif  // BIOCAD_HAVE_SCIENCE

}  // namespace

std::string htmlToText(const std::string& html, std::size_t maxChars) {
    std::string s = html;
    removeBlocks(s, "<script", "</script>");
    removeBlocks(s, "<style", "</style>");
    std::string stripped;
    stripped.reserve(s.size());
    bool inTag = false;
    for (char c : s) {
        if (c == '<') {
            inTag = true;
        } else if (c == '>') {
            inTag = false;
            stripped.push_back(' ');
        } else if (!inTag) {
            stripped.push_back(c);
        }
    }
    std::string out = collapseWs(decodeEntities(stripped));
    // trim
    const size_t b = out.find_first_not_of(' ');
    const size_t e = out.find_last_not_of(' ');
    if (b == std::string::npos) return {};
    out = out.substr(b, e - b + 1);
    if (out.size() > maxChars) {
        out.resize(maxChars);
        out += " ...[truncated]";
    }
    return out;
}

std::vector<WebHit> parseDuckDuckGoHtml(const std::string& html, int maxResults) {
    std::vector<WebHit> hits;
    size_t pos = 0;
    while (static_cast<int>(hits.size()) < maxResults) {
        const size_t marker = html.find("result__a", pos);
        if (marker == std::string::npos) break;
        const size_t tagStart = html.rfind("<a", marker);
        const size_t tagEnd = html.find('>', marker);
        if (tagStart == std::string::npos || tagEnd == std::string::npos) {
            pos = marker + 9;
            continue;
        }
        const std::string openTag = html.substr(tagStart, tagEnd - tagStart);
        const std::string url = unwrapDdgHref(tagAttr(openTag, "href"));
        const size_t innerEnd = html.find("</a>", tagEnd);
        std::string title;
        if (innerEnd != std::string::npos)
            title = htmlToText(html.substr(tagEnd + 1, innerEnd - tagEnd - 1), 300);

        std::string snippet;
        const size_t snipFrom = (innerEnd == std::string::npos) ? tagEnd : innerEnd;
        const size_t snip = html.find("result__snippet", snipFrom);
        if (snip != std::string::npos) {
            const size_t st = html.find('>', snip);
            size_t se = html.find("</a>", st == std::string::npos ? snip : st);
            const size_t td = html.find("</td>", st == std::string::npos ? snip : st);
            if (td != std::string::npos && (se == std::string::npos || td < se)) se = td;
            if (st != std::string::npos && se != std::string::npos)
                snippet = htmlToText(html.substr(st + 1, se - st - 1), 400);
        }

        if (!url.empty() && !title.empty()) hits.push_back({title, url, snippet});
        pos = (innerEnd == std::string::npos) ? tagEnd + 1 : innerEnd + 4;
    }
    return hits;
}

bool webToolsAvailable() {
#ifdef BIOCAD_HAVE_SCIENCE
    return true;
#else
    return false;
#endif
}

WebSearchResult webSearch(const std::string& query, int maxResults) {
    WebSearchResult r;
#ifndef BIOCAD_HAVE_SCIENCE
    (void)maxResults;
    (void)query;
    r.error = "web search unavailable: this build has no networking (rebuild with the science feature).";
    return r;
#else
    if (query.empty()) {
        r.error = "empty query";
        return r;
    }
    if (maxResults < 1) maxResults = 1;
    if (maxResults > 12) maxResults = 12;
    const fs::path cache = webCacheDir() / (hashHex("search:" + query) + ".json");
    if (freshEnough(cache, std::chrono::hours(12))) {
        try {
            const auto j = nlohmann::json::parse(readFile(cache));
            for (const auto& h : j)
                r.hits.push_back(
                    {h.value("title", ""), h.value("url", ""), h.value("snippet", "")});
            if (static_cast<int>(r.hits.size()) > maxResults) r.hits.resize(maxResults);
            r.ok = true;
            r.fromCache = true;
            return r;
        } catch (...) {
        }
    }
    std::string body, err;
    long code = 0;
    const std::string url = "https://html.duckduckgo.com/html/?q=" + urlEscape(query);
    if (!httpGet(url, body, err, code)) {
        r.error = "search request failed: " + err;
        return r;
    }
    if (code >= 400) {
        r.error = "search returned HTTP " + std::to_string(code);
        return r;
    }
    r.hits = parseDuckDuckGoHtml(body, maxResults);
    r.ok = true;
    nlohmann::json j = nlohmann::json::array();
    for (const auto& h : r.hits)
        j.push_back({{"title", h.title}, {"url", h.url}, {"snippet", h.snippet}});
    writeFile(cache, j.dump());
    return r;
#endif
}

WebFetchResult webFetch(const std::string& url, std::size_t maxChars, bool renderJs) {
    WebFetchResult r;
    r.finalUrl = url;
#ifndef BIOCAD_HAVE_SCIENCE
    (void)maxChars;
    (void)renderJs;
    r.error = "web fetch unavailable: this build has no networking (rebuild with the science feature).";
    return r;
#else
    if (url.rfind("http", 0) != 0) {
        r.error = "only http(s) URLs are supported";
        return r;
    }
    // Rendered (JS) and raw fetches are cached under distinct keys so one never
    // overwrites the other.
    const std::string keyPrefix = renderJs ? "fetch-js:" : "fetch:";
    const fs::path cache = webCacheDir() / (hashHex(keyPrefix + url) + ".txt");
    if (freshEnough(cache, std::chrono::hours(12))) {
        r.text = readFile(cache);
        r.ok = true;
        r.fromCache = true;
        return r;
    }
#ifdef BIOCAD_HAVE_WEBVIEW2
    if (renderJs) {
        rateLimit();  // stay polite on the rendered path too
        if (auto html = webFetchRenderedHtml(url)) {
            r.title = extractTitle(*html);
            r.text = htmlToText(*html, maxChars);
            if (!r.text.empty()) {
                r.ok = true;
                writeFile(cache, r.text);
                return r;
            }
        }
        // no runtime / navigation failure / empty result -> fall through to curl below.
    }
#else
    (void)renderJs;
#endif
    std::string body, err;
    long code = 0;
    if (!httpGet(url, body, err, code)) {
        r.error = "fetch failed: " + err;
        return r;
    }
    if (code >= 400) {
        r.error = "fetch returned HTTP " + std::to_string(code);
        return r;
    }
    r.title = extractTitle(body);
    r.text = htmlToText(body, maxChars);
    r.ok = true;
    writeFile(cache, r.text);
    return r;
#endif
}

}  // namespace biocad::agent
