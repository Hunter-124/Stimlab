// core/Uint128.h - minimal portable unsigned 128-bit arithmetic.
#pragma once

#include <cstdint>

namespace biocad::core {

// Two-limb integer with arithmetic modulo 2^128. BioCAD only needs addition and
// multiplication by a 64-bit value for its pinned PCG streams; keeping that narrow
// avoids compiler-specific __int128 while preserving the exact bit sequence.
struct Uint128 {
    std::uint64_t high = 0;
    std::uint64_t low = 0;
};

constexpr Uint128 makeUint128(std::uint64_t high, std::uint64_t low) noexcept {
    return {high, low};
}

constexpr Uint128 add(Uint128 a, Uint128 b) noexcept {
    const std::uint64_t low = a.low + b.low;
    return {a.high + b.high + (low < a.low ? 1ULL : 0ULL), low};
}

constexpr Uint128 multiply(Uint128 value, std::uint64_t factor) noexcept {
    // Split the low-limb product into 32-bit partial products. Only its carry into
    // the high limb matters; value.high * factor naturally wraps modulo 2^64.
    const std::uint64_t a0 = static_cast<std::uint32_t>(value.low);
    const std::uint64_t a1 = value.low >> 32;
    const std::uint64_t b0 = static_cast<std::uint32_t>(factor);
    const std::uint64_t b1 = factor >> 32;
    const std::uint64_t p00 = a0 * b0;
    const std::uint64_t p01 = a0 * b1;
    const std::uint64_t p10 = a1 * b0;
    const std::uint64_t p11 = a1 * b1;
    const std::uint64_t middle = (p00 >> 32) + static_cast<std::uint32_t>(p01) +
                                 static_cast<std::uint32_t>(p10);
    const std::uint64_t low = (middle << 32) | static_cast<std::uint32_t>(p00);
    const std::uint64_t carry = p11 + (p01 >> 32) + (p10 >> 32) + (middle >> 32);
    return {value.high * factor + carry, low};
}

constexpr Uint128 shiftLeftOne(Uint128 value) noexcept {
    return {(value.high << 1) | (value.low >> 63), value.low << 1};
}

}  // namespace biocad::core
