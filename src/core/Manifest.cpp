#include "core/Manifest.h"

#include <fstream>
#include <system_error>

#include "core/Hash.h"

namespace biocad {

namespace fs = std::filesystem;

std::string hashFile(const fs::path& p) {
    std::ifstream in(p, std::ios::binary);
    if (!in) return {};
    std::uint64_t h = kFnvOffset;
    char buf[64 * 1024];
    while (in) {
        in.read(buf, sizeof(buf));
        const std::streamsize n = in.gcount();
        if (n > 0) h = fnv1a(buf, static_cast<std::size_t>(n), h);
    }
    static constexpr char kHex[] = "0123456789abcdef";
    std::string out(16, '0');
    for (int i = 15; i >= 0; --i) {
        out[static_cast<std::size_t>(i)] = kHex[h & 0xF];
        h >>= 4;
    }
    return out;
}

bool Manifest::addFromFile(const std::string& id, const std::string& kind, const fs::path& path) {
    std::error_code ec;
    if (!fs::exists(path, ec) || !fs::is_regular_file(path, ec)) return false;
    const auto sz = static_cast<long long>(fs::file_size(path, ec));
    if (ec || sz <= 0) return false;
    components_.push_back({id, kind, path.string(), sz, hashFile(path)});
    return true;
}

ManifestStatus Manifest::verify() const {
    ManifestStatus s;
    s.total = static_cast<int>(components_.size());
    std::error_code ec;
    for (const auto& c : components_) {
        const fs::path p(c.path);
        if (!fs::exists(p, ec) || !fs::is_regular_file(p, ec)) {
            s.missing.push_back(c.id);
            continue;
        }
        const auto sz = static_cast<long long>(fs::file_size(p, ec));
        bool ok = !ec && (c.size == 0 || sz == c.size);
        if (ok && !c.hash.empty()) ok = (hashFile(p) == c.hash);
        if (ok) ++s.present;
        else s.corrupt.push_back(c.id);
    }
    return s;
}

int Manifest::heal() const {
    int removed = 0;
    std::error_code ec;
    for (const auto& c : components_) {
        const fs::path p(c.path);
        if (!fs::exists(p, ec) || !fs::is_regular_file(p, ec)) continue;
        const auto sz = static_cast<long long>(fs::file_size(p, ec));
        bool ok = !ec && (c.size == 0 || sz == c.size);
        if (ok && !c.hash.empty()) ok = (hashFile(p) == c.hash);
        if (!ok) {
            fs::remove(p, ec);
            if (!ec) ++removed;
        }
    }
    return removed;
}

nlohmann::json Manifest::toJson() const {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& c : components_)
        arr.push_back({{"id", c.id}, {"kind", c.kind}, {"path", c.path},
                       {"size", c.size}, {"hash", c.hash}});
    return {{"version", 1}, {"components", arr}};
}

bool Manifest::save(const fs::path& p) const {
    std::ofstream out(p, std::ios::binary);
    if (!out) return false;
    out << toJson().dump(2);
    return static_cast<bool>(out);
}

Manifest Manifest::fromJson(const nlohmann::json& j) {
    Manifest m;
    if (j.contains("components") && j["components"].is_array()) {
        for (const auto& c : j["components"])
            m.components_.push_back({c.value("id", ""), c.value("kind", ""), c.value("path", ""),
                                     c.value("size", 0LL), c.value("hash", "")});
    }
    return m;
}

Manifest Manifest::load(const fs::path& p) {
    std::ifstream in(p, std::ios::binary);
    if (!in) return {};
    try {
        nlohmann::json j;
        in >> j;
        return fromJson(j);
    } catch (...) {
        return {};
    }
}

}  // namespace biocad
