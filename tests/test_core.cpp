#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <string>

#include "core/AppPaths.h"
#include "core/Config.h"
#include "core/Error.h"
#include "core/EventBus.h"
#include "core/Hash.h"
#include "core/Secrets.h"
#include "data/Domain.h"

using namespace biocad;

TEST_CASE("Hash is deterministic and fixed width", "[core][hash]") {
    REQUIRE(hashHex("amphetamine") == hashHex("amphetamine"));
    REQUIRE(hashHex("amphetamine") != hashHex("methamphetamine"));
    REQUIRE(hashHex("anything").size() == 16);
    REQUIRE(hash64("abc") == hash64("abc"));
}

TEST_CASE("Result<T> carries value or error", "[core][result]") {
    Result<int> ok(42);
    REQUIRE(ok.ok());
    REQUIRE(ok.value() == 42);
    REQUIRE(ok.valueOr(-1) == 42);

    Result<int> bad(Error::notFound("missing"));
    REQUIRE_FALSE(bad.ok());
    REQUIRE(bad.error().code == Error::Code::NotFound);
    REQUIRE(bad.valueOr(-1) == -1);

    Status okv;
    REQUIRE(okv.ok());
    Status errv(Error::io("disk"));
    REQUIRE_FALSE(errv.ok());
    REQUIRE(std::string(errv.error().codeName()) == "Io");
}

TEST_CASE("EventBus delivers typed events", "[core][events]") {
    struct Ping { int n; };
    EventBus bus;
    int total = 0;
    auto tok = bus.subscribe<Ping>([&](const Ping& p) { total += p.n; });
    bus.publish(Ping{3});
    bus.publish(Ping{4});
    REQUIRE(total == 7);
    bus.unsubscribe(tok);
    bus.publish(Ping{100});
    REQUIRE(total == 7);
}

TEST_CASE("Config round-trips through disk", "[core][config]") {
    const auto file = std::filesystem::temp_directory_path() / "biocad_test_config.json";
    std::error_code ec;
    std::filesystem::remove(file, ec);
    {
        Config c(file);
        c.load();
        c.set<std::string>("provider", "anthropic");
        c.set<int>("threads", 8);
        REQUIRE(c.save());
    }
    {
        Config c(file);
        c.load();
        REQUIRE(c.has("provider"));
        REQUIRE(c.get<std::string>("provider", "") == "anthropic");
        REQUIRE(c.get<int>("threads", 0) == 8);
        REQUIRE(c.get<int>("missing", -1) == -1);
    }
    std::filesystem::remove(file, ec);
}

TEST_CASE("AppPaths resolves a BioCAD root", "[core][paths]") {
    const auto& root = AppPaths::instance().root();
    REQUIRE_FALSE(root.empty());
    REQUIRE(root.filename().string() == "BioCAD");
    REQUIRE(AppPaths::instance().db().filename().string() == "biocad.db");
}

TEST_CASE("AppPaths migrates a legacy StimLab root exactly once", "[core][paths]") {
    namespace fs = std::filesystem;
    std::error_code ec;
    const auto base = fs::temp_directory_path() / "biocad_migration_test";
    fs::remove_all(base, ec);
    const auto legacy = base / "StimLab";
    const auto target = base / "BioCAD";

    fs::create_directories(legacy / "runtime" / "receptors", ec);
    fs::create_directories(legacy / "logs", ec);
    std::ofstream(legacy / "stimlab.db") << "db";
    std::ofstream(legacy / "logs" / "stimlab.log") << "log";
    std::ofstream(legacy / "runtime" / "receptors" / "x.pdbqt") << "receptor";
    std::ofstream(legacy / "runtime" / "vina.exe") << "engine";

    const AppPaths paths(target);
    REQUIRE(paths.ensureLayout());

    REQUIRE(fs::exists(target));
    REQUIRE_FALSE(fs::exists(legacy));
    REQUIRE(fs::exists(target / "biocad.db"));
    REQUIRE_FALSE(fs::exists(target / "stimlab.db"));
    REQUIRE(fs::exists(target / "logs" / "biocad.log"));
    REQUIRE(fs::exists(target / "runtime" / "vina.exe"));
    REQUIRE_FALSE(fs::exists(target / "runtime" / "receptors"));

    // Second call is a no-op: the root now exists, so nothing is moved or deleted.
    fs::create_directories(base / "StimLab", ec);
    std::ofstream(base / "StimLab" / "stimlab.db") << "second";
    REQUIRE_FALSE(paths.migrateLegacyRoot());
    REQUIRE(fs::exists(base / "StimLab" / "stimlab.db"));
    REQUIRE(fs::exists(target / "runtime" / "vina.exe"));

    fs::remove_all(base, ec);
}

TEST_CASE("Secrets round-trip (DPAPI on Windows, base64 elsewhere)", "[core][secrets]") {
    const std::string secret = "sk-test-12345";
    auto enc = Secrets::protect(secret);
    REQUIRE(enc.ok());
    REQUIRE(enc.value() != secret);
    auto dec = Secrets::unprotect(enc.value());
    REQUIRE(dec.ok());
    REQUIRE(dec.value() == secret);
}

TEST_CASE("Molecule serializes to/from JSON", "[data][json]") {
    Molecule m;
    m.id = "caffeine";
    m.name = "Caffeine";
    m.molWeight = 194.19;
    m.logP = -0.07;
    const nlohmann::json j = m;
    const auto back = j.get<Molecule>();
    REQUIRE(back.id == "caffeine");
    REQUIRE(back.name == "Caffeine");
    REQUIRE(back.molWeight == 194.19);
}
