#include "modules/docking/EngineLocator.h"

#include <array>
#include <cstdlib>
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
namespace {

namespace fs = std::filesystem;

// Official AutoDock Vina 1.2.5 Windows release asset. The SHA-256 is left empty
// here: if a build pins the exact published hash, set kVinaSha256 and verification
// becomes mandatory; with an empty hash the download is accepted but flagged
// "unverified" in the note. Either way provisioning stays best-effort/non-blocking.
constexpr const char* kVinaUrl =
    "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_win.exe";
constexpr const char* kVinaSha256 = "";  // optional pin; empty = accept unverified

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

}  // namespace

fs::path enginesDir() {
    return AppPaths::instance().runtime() / "engines";
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
    // Best-effort PowerShell download + optional SHA-256 verify. Quoted so spaces
    // in the profile path are safe; -ErrorAction Stop turns transfer failures into
    // a non-zero exit we can detect. We NEVER block app startup on this (callers
    // run it off the UI thread / opportunistically), and any failure leaves us in
    // the locate-only state with a note.
    std::string ps;
    ps += "$ErrorActionPreference='Stop';";
    ps += "try{";
    ps += "[Net.ServicePointManager]::SecurityProtocol=[Net.SecurityProtocolType]::Tls12;";
    ps += "Invoke-WebRequest -UseBasicParsing -Uri '";
    ps += kVinaUrl;
    ps += "' -OutFile '";
    ps += target.string();
    ps += "';";
    if (kVinaSha256 && kVinaSha256[0] != '\0') {
        ps += "$h=(Get-FileHash -Algorithm SHA256 '";
        ps += target.string();
        ps += "').Hash;";
        ps += "if($h -ne '";
        ps += kVinaSha256;
        ps += "'){Remove-Item -Force '";
        ps += target.string();
        ps += "';exit 3};";
    }
    ps += "exit 0";
    ps += "}catch{exit 2}";

    std::string cmd = "powershell -NoProfile -NonInteractive -Command \"" + ps + "\"";
    const int code = std::system(cmd.c_str());

    if (code == 0 && fs::exists(target, ec) && fs::file_size(target, ec) > 0) {
        r.fetched = true;
        r.path = target.string();
        r.note = (kVinaSha256 && kVinaSha256[0] != '\0')
                     ? "downloaded + SHA-256 verified vina.exe"
                     : "downloaded vina.exe (hash unverified - no pinned SHA-256)";
        return r;
    }
    r.fetched = false;
    r.note = "vina download failed (exit " + std::to_string(code) +
             "); staying in descriptor-estimate fallback.";
    // Remove a partial/zero-byte file so a later locate() does not trust it.
    if (fs::exists(target, ec) && fs::file_size(target, ec) == 0) fs::remove(target, ec);
    return r;
#else
    (void)target;
    r.note = "provisioning only implemented on Windows.";
    return r;
#endif
}

}  // namespace stimlab::docking
