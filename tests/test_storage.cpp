#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <string>

#include "storage/RunStore.h"

using namespace stimlab;

namespace {
std::filesystem::path scratchDb(const char* tag) {
    auto p = std::filesystem::temp_directory_path() / (std::string("stimlab_test_") + tag + ".db");
    std::error_code ec;
    std::filesystem::remove(p, ec);
    std::filesystem::remove(std::filesystem::path(p.string() + "-wal"), ec);
    std::filesystem::remove(std::filesystem::path(p.string() + "-shm"), ec);
    return p;
}
}  // namespace

TEST_CASE("SqliteRunStore round-trips runs and orders newest first", "[storage][runs]") {
    const auto db = scratchDb("roundtrip");
    {
        SqliteRunStore store(db);
        REQUIRE(store.ok());
        REQUIRE(store.count() == 0);

        store.record({"", "Docking", "Amphetamine -> DAT", "complete", "2026-06-08 10:00", "best -7.2"});
        store.record({"", "ADMET", "MDMA", "complete", "2026-06-08 10:05", "WARN"});
        REQUIRE(store.count() == 2);

        const auto rows = store.recent();
        REQUIRE(rows.size() == 2);
        REQUIRE(rows.front().subject == "MDMA");          // newest first
        REQUIRE(rows.front().kind == "ADMET");
        REQUIRE(rows.back().subject == "Amphetamine -> DAT");
        REQUIRE(rows.front().id.rfind("run-", 0) == 0);   // formatted id
    }
    // Reopen: data persisted across store instances (and a default status filled).
    {
        SqliteRunStore reopened(db);
        REQUIRE(reopened.count() == 2);
        reopened.record({"", "Stability", "Cocaine", "", "", "ester hydrolysis"});
        REQUIRE(reopened.count() == 3);
        REQUIRE(reopened.recent().front().status == "complete");  // empty status defaulted
    }
    std::error_code ec;
    std::filesystem::remove(db, ec);
}
