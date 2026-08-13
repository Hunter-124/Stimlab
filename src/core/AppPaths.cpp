#include "core/AppPaths.h"

#include <cstdlib>
#include <system_error>

#include <spdlog/spdlog.h>

#if defined(_WIN32)
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#  include <shlobj.h>
#endif

namespace biocad {

namespace {

std::filesystem::path resolveRoamingAppData() {
    // 1) Explicit override for tests / portable installs.
    if (const char* home = std::getenv("BIOCAD_HOME"); home && *home) {
        return std::filesystem::path(home);
    }

#if defined(_WIN32)
    PWSTR raw = nullptr;
    if (SUCCEEDED(SHGetKnownFolderPath(FOLDERID_RoamingAppData, 0, nullptr, &raw)) && raw) {
        std::filesystem::path p(raw);
        CoTaskMemFree(raw);
        return p / "BioCAD";
    }
    if (const char* appdata = std::getenv("APPDATA"); appdata && *appdata) {
        return std::filesystem::path(appdata) / "BioCAD";
    }
#endif

    // POSIX / fallback.
    if (const char* xdg = std::getenv("XDG_DATA_HOME"); xdg && *xdg) {
        return std::filesystem::path(xdg) / "BioCAD";
    }
    if (const char* hm = std::getenv("HOME"); hm && *hm) {
        return std::filesystem::path(hm) / ".local" / "share" / "BioCAD";
    }
    return std::filesystem::temp_directory_path() / "BioCAD";
}

}  // namespace

AppPaths::AppPaths() : root_(resolveRoamingAppData()) {}

AppPaths& AppPaths::instance() {
    static AppPaths inst;
    return inst;
}

bool AppPaths::migrateLegacyRoot() const {
    std::error_code ec;
    if (std::filesystem::exists(root_, ec)) return false;

    const auto legacy = root_.parent_path() / "StimLab";
    if (!std::filesystem::exists(legacy, ec)) return false;

    try {
        std::filesystem::rename(legacy, root_);
    } catch (const std::filesystem::filesystem_error& e) {
        // Never abort startup over a migration failure: the app simply starts fresh.
        spdlog::warn("legacy StimLab root migration failed: {}", e.what());
        return false;
    }

    // Cached receptors carry the old REMARK STIMLAB_BOX marker and must be re-prepared.
    std::filesystem::remove_all(root_ / "runtime" / "receptors", ec);

    const auto renameIfPresent = [&](const std::filesystem::path& from,
                                     const std::filesystem::path& to) {
        std::error_code inner;
        if (std::filesystem::exists(from, inner)) {
            std::filesystem::rename(from, to, inner);
        }
    };
    renameIfPresent(root_ / "stimlab.db", root_ / "biocad.db");
    renameIfPresent(root_ / "logs" / "stimlab.log", root_ / "logs" / "biocad.log");
    return true;
}

bool AppPaths::ensureLayout() const {
    migrateLegacyRoot();

    std::error_code ec;
    for (const auto& dir : {root_, artifacts(), runtime(), presets(), logs(), cache()}) {
        std::filesystem::create_directories(dir, ec);
        if (ec) return false;
    }
    return true;
}

}  // namespace biocad
