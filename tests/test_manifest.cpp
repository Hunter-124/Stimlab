#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <string>

#include "core/Manifest.h"

using namespace biocad;
namespace fs = std::filesystem;

namespace {
fs::path scratch(const char* tag) {
    return fs::temp_directory_path() / (std::string("biocad_mf_") + tag);
}
void writeText(const fs::path& p, const std::string& s) {
    std::ofstream o(p, std::ios::binary);
    o << s;
}
}  // namespace

TEST_CASE("Manifest verifies, detects corruption, and self-heals", "[core][manifest]") {
    const fs::path a = scratch("a.bin"), b = scratch("b.bin");
    std::error_code ec;
    fs::remove(a, ec);
    fs::remove(b, ec);
    writeText(a, "engine-binary-bytes");
    writeText(b, "receptor-pdbqt-bytes");

    Manifest m;
    REQUIRE(m.addFromFile("engine:vina", "engine", a));
    REQUIRE(m.addFromFile("receptor:DAT", "receptor", b));
    REQUIRE_FALSE(m.addFromFile("missing", "engine", scratch("nope.bin")));  // absent -> not added
    REQUIRE(m.components().size() == 2);

    // All present + intact.
    auto s = m.verify();
    REQUIRE(s.complete());
    REQUIRE(s.present == 2);

    // Corrupt b (changes size + content hash) -> reported corrupt, not present.
    writeText(b, "receptor-pdbqt-bytes-TAMPERED-LONGER");
    s = m.verify();
    REQUIRE_FALSE(s.complete());
    REQUIRE(s.present == 1);
    REQUIRE(s.corrupt.size() == 1);
    REQUIRE(s.corrupt.front() == "receptor:DAT");

    // heal() deletes the corrupt file so a later provision re-fetches it.
    REQUIRE(m.heal() == 1);
    REQUIRE_FALSE(fs::exists(b, ec));
    s = m.verify();
    REQUIRE(s.missing.size() == 1);
    REQUIRE(s.missing.front() == "receptor:DAT");

    fs::remove(a, ec);
}

TEST_CASE("Manifest round-trips through JSON on disk", "[core][manifest]") {
    const fs::path f = scratch("rt.bin"), mfPath = scratch("manifest.json");
    std::error_code ec;
    writeText(f, "some-bytes");

    Manifest m;
    REQUIRE(m.addFromFile("engine:vina", "engine", f));
    REQUIRE(m.save(mfPath));

    const Manifest loaded = Manifest::load(mfPath);
    REQUIRE(loaded.components().size() == 1);
    REQUIRE(loaded.components().front().id == "engine:vina");
    REQUIRE(loaded.verify().complete());  // file still intact

    fs::remove(f, ec);
    fs::remove(mfPath, ec);
    // Loading an absent manifest yields an empty one (no throw).
    REQUIRE(Manifest::load(scratch("does_not_exist.json")).empty());
}
