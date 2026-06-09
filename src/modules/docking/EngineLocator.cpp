#include "modules/docking/EngineLocator.h"

#include <array>
#include <atomic>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <string>
#include <system_error>
#include <vector>

#include "core/AppPaths.h"

#if defined(_WIN32)
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#endif

namespace stimlab::docking {

// Process-wide compute-mode preference (a single global user setting). Thread-safe.
namespace {
std::atomic<int> g_computeMode{static_cast<int>(ComputeMode::Auto)};
}  // namespace
void setComputeMode(ComputeMode m) {
    g_computeMode.store(static_cast<int>(m), std::memory_order_relaxed);
}
ComputeMode computeMode() {
    return static_cast<ComputeMode>(g_computeMode.load(std::memory_order_relaxed));
}

namespace {

namespace fs = std::filesystem;

// Official AutoDock Vina 1.2.5 Windows release asset. GitHub publishes no per-asset
// digest for this release, so kVinaSha256 is left empty by default; integrity can be
// pinned at runtime (env STIMLAB_VINA_SHA256 or a `vina.sha256` file in the engines
// dir) WITHOUT a rebuild. The published byte size IS known and is always enforced as
// a sanity check, so an HTML error page / truncated transfer is rejected even when no
// hash is configured. Provisioning stays best-effort and non-blocking either way.
constexpr const char* kVinaUrl =
    "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_win.exe";
constexpr const char* kVinaSha256 = "";          // optional compile-time pin
constexpr long long   kVinaSizeBytes = 1203712;  // published size of vina_1.2.5_win.exe
constexpr long long   kVinaMinSizeBytes = 200000;  // reject anything implausibly small

// Candidate file names per engine (most specific first).
std::vector<std::string> candidateNames(Engine e) {
    switch (e) {
        case Engine::Vina:
            return {"vina.exe", "vina_1.2.5_win.exe", "vina_1.2.5.exe",
                    "vina_1.2.3_win.exe", "vina_win.exe", "vina"};
        case Engine::Smina:
            return {"smina.exe", "smina_win.exe", "smina.static.exe", "smina"};
        case Engine::Obabel:
            return {"obabel.exe", "obabel"};
    }
    return {};
}

// Split a PATH-style environment string on ';' (Windows) into directories.
std::vector<fs::path> pathDirs() {
    std::vector<fs::path> dirs;
#if defined(_WIN32)
    if (const char* path = std::getenv("PATH"); path && *path) {
        std::string cur;
        for (const char* p = path; ; ++p) {
            if (*p == ';' || *p == '\0') {
                if (!cur.empty()) {
                    std::error_code ec;
                    fs::path d(cur);
                    if (fs::exists(d, ec)) dirs.push_back(d);
                }
                cur.clear();
                if (*p == '\0') break;
            } else {
                cur.push_back(*p);
            }
        }
    }
#endif
    return dirs;
}

std::string trimHash(std::string s) {
    // Keep only hex digits (Get-FileHash emits uppercase; a file may have a newline).
    std::string out;
    for (char c : s) if (std::isxdigit(static_cast<unsigned char>(c))) out.push_back(c);
    return out;
}

std::string readFirstLine(const fs::path& p) {
    std::error_code ec;
    if (!fs::exists(p, ec)) return {};
    std::ifstream in(p);
    std::string line;
    std::getline(in, line);
    return line;
}

}  // namespace

fs::path enginesDir() {
    return AppPaths::instance().runtime() / "engines";
}

std::string expectedVinaSha256() {
    // Runtime override precedence: env var, then a vina.sha256 file beside the binary,
    // then the compile-time pin. All normalised to bare uppercase hex.
    if (const char* env = std::getenv("STIMLAB_VINA_SHA256"); env && *env) {
        std::string h = trimHash(env);
        if (h.size() == 64) {
            for (auto& c : h) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            return h;
        }
    }
    std::string fileHash = trimHash(readFirstLine(enginesDir() / "vina.sha256"));
    if (fileHash.size() == 64) {
        for (auto& c : fileHash) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        return fileHash;
    }
    std::string pin = trimHash(kVinaSha256);
    if (pin.size() == 64) {
        for (auto& c : pin) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        return pin;
    }
    return {};
}

std::optional<fs::path> locateEngine(Engine e) {
    std::error_code ec;
    const auto names = candidateNames(e);

    // 1) runtime/engines (preferred, app-provisioned).
    const fs::path dir = enginesDir();
    for (const auto& n : names) {
        const fs::path cand = dir / n;
        if (fs::exists(cand, ec) && fs::is_regular_file(cand, ec)) return cand;
    }
    // Also accept any vina_*/smina_* file dropped into engines without an exact name.
    if (fs::exists(dir, ec) && fs::is_directory(dir, ec)) {
        const std::string stem = (e == Engine::Vina) ? "vina" :
                                 (e == Engine::Smina) ? "smina" : "obabel";
        for (fs::directory_iterator it(dir, ec), end; it != end && !ec; it.increment(ec)) {
            if (!it->is_regular_file(ec)) continue;
            std::string fn = it->path().filename().string();
            std::string low = fn;
            for (auto& c : low) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            if (low.rfind(stem, 0) == 0 &&
                (low.size() >= 4 && low.substr(low.size() - 4) == ".exe"))
                return it->path();
        }
    }

    // 2) PATH.
    for (const auto& d : pathDirs()) {
        for (const auto& n : names) {
            const fs::path cand = d / n;
            if (fs::exists(cand, ec) && fs::is_regular_file(cand, ec)) return cand;
        }
    }
    return std::nullopt;
}

bool engineAvailable(Engine e) { return locateEngine(e).has_value(); }

ProvisionResult ensureVina(bool allowDownload) {
    ProvisionResult r;

    // Already provisioned? Done.
    if (auto p = locateEngine(Engine::Vina)) {
        r.fetched = true;
        r.path = p->string();
        r.note = "vina already present at " + r.path;
        return r;
    }
    if (!allowDownload) {
        r.fetched = false;
        r.note = "vina not found; download not requested (locate-only).";
        return r;
    }

    std::error_code ec;
    const fs::path dir = enginesDir();
    fs::create_directories(dir, ec);
    if (ec) {
        r.note = "could not create engines dir: " + ec.message();
        return r;
    }
    const fs::path target = dir / "vina.exe";

#if defined(_WIN32)
    // Best-effort PowerShell download with a size sanity-check + optional SHA-256
    // verify. Quoted so spaces in the profile path are safe; -ErrorAction Stop turns
    // transfer failures into a non-zero exit we can detect. We NEVER block app
    // startup on this (callers run it off the UI thread / opportunistically), and any
    // failure leaves us in the locate-only state with a note. Exit codes: 0 = ok,
    // 2 = transfer error, 3 = hash mismatch, 4 = implausibly small (error page).
    const std::string expected = expectedVinaSha256();
    std::string ps;
    ps += "$ErrorActionPreference='Stop';";
    ps += "try{";
    ps += "[Net.ServicePointManager]::SecurityProtocol=[Net.SecurityProtocolType]::Tls12;";
    ps += "Invoke-WebRequest -UseBasicParsing -Uri '";
    ps += kVinaUrl;
    ps += "' -OutFile '";
    ps += target.string();
    ps += "';";
    // Reject an HTML error page / truncated transfer masquerading as the binary.
    ps += "if((Get-Item '";
    ps += target.string();
    ps += "').Length -lt ";
    ps += std::to_string(kVinaMinSizeBytes);
    ps += "){Remove-Item -Force '";
    ps += target.string();
    ps += "';exit 4};";
    if (expected.size() == 64) {
        ps += "$h=(Get-FileHash -Algorithm SHA256 '";
        ps += target.string();
        ps += "').Hash;";
        ps += "if($h -ne '";
        ps += expected;
        ps += "'){Remove-Item -Force '";
        ps += target.string();
        ps += "';exit 3};";
    }
    ps += "exit 0";
    ps += "}catch{exit 2}";

    std::string cmd = "powershell -NoProfile -NonInteractive -Command \"" + ps + "\"";
    const int code = std::system(cmd.c_str());

    if (code == 0 && fs::exists(target, ec) && fs::file_size(target, ec) >= kVinaMinSizeBytes) {
        r.fetched = true;
        r.path = target.string();
        const long long sz = static_cast<long long>(fs::file_size(target, ec));
        const bool sizeMatch = (sz == kVinaSizeBytes);
        if (expected.size() == 64)
            r.note = "downloaded + SHA-256 verified vina.exe";
        else if (sizeMatch)
            r.note = "downloaded vina.exe (size matches published " +
                     std::to_string(kVinaSizeBytes) + " B; no pinned SHA-256 - "
                     "set STIMLAB_VINA_SHA256 or runtime/engines/vina.sha256 to verify)";
        else
            r.note = "downloaded vina.exe (" + std::to_string(sz) +
                     " B; UNVERIFIED - no pinned SHA-256 and size differs from published)";
        return r;
    }
    r.fetched = false;
    r.note = code == 3   ? "vina download SHA-256 mismatch; discarded (check your pinned hash)."
             : code == 4 ? "vina download too small (likely an error page); discarded."
                         : "vina download failed (exit " + std::to_string(code) +
                               "); staying in descriptor-estimate fallback.";
    // Remove a partial/zero-byte file so a later locate() does not trust it.
    if (fs::exists(target, ec) && fs::file_size(target, ec) < kVinaMinSizeBytes)
        fs::remove(target, ec);
    return r;
#else
    (void)target;
    r.note = "provisioning only implemented on Windows.";
    return r;
#endif
}

ProvisionResult ensureObabel() {
    ProvisionResult r;
    if (auto p = locateEngine(Engine::Obabel)) {
        r.fetched = true;
        r.path = p->string();
        r.note = "obabel located at " + r.path + " (used for receptor protonation).";
    } else {
        r.fetched = false;
        r.note = "obabel.exe not found; receptor prep uses the built-in heavy-atom writer. "
                 "Drop obabel.exe into runtime/engines for full protonation.";
    }
    return r;
}

}  // namespace stimlab::docking
