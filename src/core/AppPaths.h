// core/AppPaths.h - resolves and creates the %APPDATA%/BioCAD/* layout.
// All persistent application state lives under one root; nothing in Program Files.
#pragma once

#include <filesystem>

namespace biocad {

class AppPaths {
public:
    // Singleton: resolves the root once (env override BIOCAD_HOME wins, else %APPDATA%/BioCAD).
    static AppPaths& instance();

    // Explicit root (tests, portable installs); the singleton uses the resolved root.
    explicit AppPaths(std::filesystem::path root) : root_(std::move(root)) {}

    [[nodiscard]] const std::filesystem::path& root() const { return root_; }

    [[nodiscard]] std::filesystem::path db() const        { return root_ / "biocad.db"; }
    [[nodiscard]] std::filesystem::path artifacts() const { return root_ / "artifacts"; }
    [[nodiscard]] std::filesystem::path runtime() const   { return root_ / "runtime"; }
    [[nodiscard]] std::filesystem::path presets() const   { return root_ / "presets"; }
    [[nodiscard]] std::filesystem::path logs() const       { return root_ / "logs"; }
    [[nodiscard]] std::filesystem::path cache() const     { return root_ / "cache"; }
    [[nodiscard]] std::filesystem::path config() const    { return root_ / "config.json"; }
    [[nodiscard]] std::filesystem::path manifest() const  { return root_ / "manifest.json"; }

    // Moves a pre-BioCAD %APPDATA%/StimLab tree to the new root exactly once, so provisioned
    // docking engines (GBs) are not re-downloaded. Returns true if a migration happened.
    bool migrateLegacyRoot() const;

    // Create the full directory layout if missing. Returns false on filesystem error.
    bool ensureLayout() const;

private:
    AppPaths();
    std::filesystem::path root_;
};

}  // namespace biocad
