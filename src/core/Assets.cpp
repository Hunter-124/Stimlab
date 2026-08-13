#include "core/Assets.h"

#include <cstdlib>
#include <system_error>

#if defined(_WIN32)
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#endif

namespace biocad::core {

namespace {

std::filesystem::path exeDirectory() {
    std::filesystem::path p;
#if defined(_WIN32)
    wchar_t buf[MAX_PATH]{};
    if (GetModuleFileNameW(nullptr, buf, MAX_PATH) > 0) p = std::filesystem::path(buf).parent_path();
#else
    std::error_code ec;
    const auto self = std::filesystem::read_symlink("/proc/self/exe", ec);
    if (!ec) p = self.parent_path();
#endif
    return p;
}

}  // namespace

std::filesystem::path assetRoot() {
    std::error_code ec;
    // Explicit override, for tests and for a packaged layout that puts the data
    // somewhere else. Checked first so a broken auto-detection is always fixable
    // from outside the binary.
    if (const char* env = std::getenv("BIOCAD_ASSETS"); env && *env) {
        if (std::filesystem::is_directory(env, ec)) return std::filesystem::path(env);
    }
    const auto exeDir = exeDirectory();
    if (!exeDir.empty()) {
        // Shipped layout: assets/ beside the executable.
        if (std::filesystem::is_directory(exeDir / "assets", ec)) return exeDir / "assets";
        // Dev layout: build/<preset>/bin/BioCAD.exe, with the source tree three up.
        // Walk upwards instead of assuming a fixed depth so a test binary in
        // build/<preset>/tests/ finds the same tree the app does.
        auto up = exeDir;
        for (int i = 0; i < 6 && !up.empty(); ++i) {
            if (std::filesystem::is_directory(up / "assets" / "packs", ec)) return up / "assets";
            const auto parent = up.parent_path();
            if (parent == up) break;
            up = parent;
        }
    }
    if (std::filesystem::is_directory("assets", ec)) return "assets";
    return {};
}

std::filesystem::path assetDir(std::string_view sub) {
    const auto root = assetRoot();
    if (root.empty()) return {};
    return root / std::filesystem::path(sub);
}

}  // namespace biocad::core
