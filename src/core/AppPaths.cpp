#include "core/AppPaths.h"

#include <cstdlib>
#include <system_error>

#if defined(_WIN32)
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#  include <shlobj.h>
#endif

namespace stimlab {

namespace {

std::filesystem::path resolveRoamingAppData() {
    // 1) Explicit override for tests / portable installs.
    if (const char* home = std::getenv("STIMLAB_HOME"); home && *home) {
        return std::filesystem::path(home);
    }

#if defined(_WIN32)
    PWSTR raw = nullptr;
    if (SUCCEEDED(SHGetKnownFolderPath(FOLDERID_RoamingAppData, 0, nullptr, &raw)) && raw) {
        std::filesystem::path p(raw);
        CoTaskMemFree(raw);
        return p / "StimLab";
    }
    if (const char* appdata = std::getenv("APPDATA"); appdata && *appdata) {
        return std::filesystem::path(appdata) / "StimLab";
    }
#endif

    // POSIX / fallback.
    if (const char* xdg = std::getenv("XDG_DATA_HOME"); xdg && *xdg) {
        return std::filesystem::path(xdg) / "StimLab";
    }
    if (const char* hm = std::getenv("HOME"); hm && *hm) {
        return std::filesystem::path(hm) / ".local" / "share" / "StimLab";
    }
    return std::filesystem::temp_directory_path() / "StimLab";
}

}  // namespace

AppPaths::AppPaths() : root_(resolveRoamingAppData()) {}

AppPaths& AppPaths::instance() {
    static AppPaths inst;
    return inst;
}

bool AppPaths::ensureLayout() const {
    std::error_code ec;
    for (const auto& dir : {root_, artifacts(), runtime(), presets(), logs(), cache()}) {
        std::filesystem::create_directories(dir, ec);
        if (ec) return false;
    }
    return true;
}

}  // namespace stimlab
