// storage/RunStore.h - SQLite-backed run history (Phase D, WP-4 persistence).
//
// Persists analysis runs to %APPDATA%/StimLab/stimlab.db (WAL mode) so the Runs
// panel survives restarts. Implements the frozen IRunStore contract (recent() +
// the Phase-D record() hook). A path-taking constructor lets tests use a scratch
// database instead of the real profile.
//
// SAFETY SCOPE: stores run metadata (what was analyzed + the verdict/affinity
// summary). No synthesis/route content is stored.
#pragma once

#include <cstddef>
#include <filesystem>
#include <memory>
#include <vector>

#include "contracts/IModules.h"

namespace stimlab {

class SqliteRunStore final : public IRunStore {
public:
    SqliteRunStore();                                       // %APPDATA%/StimLab/stimlab.db
    explicit SqliteRunStore(std::filesystem::path dbPath);  // explicit DB (tests)
    ~SqliteRunStore() override;
    SqliteRunStore(const SqliteRunStore&) = delete;
    SqliteRunStore& operator=(const SqliteRunStore&) = delete;

    std::vector<RunRecord> recent() const override;   // newest first, capped
    void record(const RunRecord& r) override;          // insert (auto id + timestamp)

    [[nodiscard]] bool ok() const;                     // database opened?
    [[nodiscard]] std::size_t count() const;           // total rows

private:
    void openAndInit(const std::filesystem::path& dbPath);

    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace stimlab
