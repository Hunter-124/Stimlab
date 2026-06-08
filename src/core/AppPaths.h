// core/AppPaths.h - resolves and creates the %APPDATA%/StimLab/* layout.
// All persistent application state lives under one root; nothing in Program Files.
#pragma once

#include <filesystem>

namespace stimlab {

class AppPaths {
public:
    // Singleton: resolves the root once (env override STIMLAB_HOME wins, else %APPDATA%/StimLab).
    static AppPaths& instance();

    [[nodiscard]] const std::filesystem::path& root() const { return root_; }

    [[nodiscard]] std::filesystem::path db() const        { return root_ / "stimlab.db"; }
    [[nodiscard]] std::filesystem::path artifacts() const { return root_ / "artifacts"; }
    [[nodiscard]] std::filesystem::path runtime() const   { return root_ / "runtime"; }
    [[nodiscard]] std::filesystem::path presets() const   { return root_ / "presets"; }
    [[nodiscard]] std::filesystem::path logs() const       { return root_ / "logs"; }
    [[nodiscard]] std::filesystem::path cache() const     { return root_ / "cache"; }
    [[nodiscard]] std::filesystem::path config() const    { return root_ / "config.json"; }

    // Create the full directory layout if missing. Returns false on filesystem error.
    bool ensureLayout() const;

private:
    AppPaths();
    std::filesystem::path root_;
};

}  // namespace stimlab
