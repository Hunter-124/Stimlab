// core/Manifest.h - the provisioned-runtime manifest + self-heal (Phase E, WP-B).
//
// A single source of truth for what BioCAD has provisioned under %APPDATA%/BioCAD/
// runtime (engine binaries, prepared receptors, future ML models) and an integrity
// check. After provisioning, the Provisioner records each component's path + size +
// content hash into manifest.json. On launch, verify() re-checks every component:
// a missing or size/hash-mismatched file is reported, and heal() deletes corrupt
// files so the next provisioning run re-fetches them.
//
// The content hash is a fast non-cryptographic FNV-1a (corruption/truncation
// detection) - download AUTHENTICITY is handled separately by the pinned SHA-256 in
// EngineLocator. There is no synthesis content here; this only tracks files on disk.
#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

namespace biocad {

// One provisioned artifact tracked for self-heal.
struct ManifestComponent {
    std::string id;        // stable key, e.g. "engine:vina", "receptor:DAT"
    std::string kind;      // "engine" | "receptor" | "model"
    std::string path;      // absolute path on disk
    long long   size = 0;  // expected size in bytes
    std::string hash;      // FNV-1a content hash (16 hex chars); "" = size-only check
};

// The outcome of verifying a manifest against the filesystem.
struct ManifestStatus {
    int                      total = 0;
    int                      present = 0;     // exist + size (+ hash) match
    std::vector<std::string> missing;         // ids absent on disk
    std::vector<std::string> corrupt;         // ids present but size/hash mismatch
    [[nodiscard]] bool complete() const { return total > 0 && present == total; }
};

// FNV-1a content hash of a file (16 hex chars); empty string if unreadable.
std::string hashFile(const std::filesystem::path& p);

class Manifest {
public:
    void add(ManifestComponent c) { components_.push_back(std::move(c)); }
    [[nodiscard]] const std::vector<ManifestComponent>& components() const { return components_; }
    [[nodiscard]] bool empty() const { return components_.empty(); }

    // Build a component by hashing the file at `path` right now (records its size +
    // content hash). Returns false (and adds nothing) if the file is missing/empty.
    bool addFromFile(const std::string& id, const std::string& kind,
                     const std::filesystem::path& path);

    // Verify every component against the filesystem (existence + size + hash).
    [[nodiscard]] ManifestStatus verify() const;

    // Delete every corrupt component file so a later provision re-fetches it.
    // Returns the number removed. (Missing files are left to provisioning.)
    int heal() const;

    [[nodiscard]] nlohmann::json toJson() const;
    bool save(const std::filesystem::path& p) const;
    static Manifest fromJson(const nlohmann::json& j);
    static Manifest load(const std::filesystem::path& p);  // empty Manifest if absent

private:
    std::vector<ManifestComponent> components_;
};

}  // namespace biocad
