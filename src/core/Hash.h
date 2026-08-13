// core/Hash.h - small, stable, dependency-free hashing (FNV-1a 64).
// Used for content-addressing and DAG cache keys. Deterministic across runs.
#pragma once

#include <cstdint>
#include <string>
#include <string_view>

namespace biocad {

inline constexpr std::uint64_t kFnvOffset = 1469598103934665603ULL;
inline constexpr std::uint64_t kFnvPrime  = 1099511628211ULL;

// FNV-1a over raw bytes.
inline std::uint64_t fnv1a(const void* data, std::size_t len, std::uint64_t seed = kFnvOffset) {
    const auto* p = static_cast<const unsigned char*>(data);
    std::uint64_t h = seed;
    for (std::size_t i = 0; i < len; ++i) {
        h ^= static_cast<std::uint64_t>(p[i]);
        h *= kFnvPrime;
    }
    return h;
}

inline std::uint64_t hash64(std::string_view s, std::uint64_t seed = kFnvOffset) {
    return fnv1a(s.data(), s.size(), seed);
}

// 16-char lowercase hex digest, suitable as an artifact key fragment.
inline std::string hashHex(std::string_view s) {
    static constexpr char kHex[] = "0123456789abcdef";
    std::uint64_t h = hash64(s);
    std::string out(16, '0');
    for (int i = 15; i >= 0; --i) {
        out[static_cast<std::size_t>(i)] = kHex[h & 0xF];
        h >>= 4;
    }
    return out;
}

}  // namespace biocad
