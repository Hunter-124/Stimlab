#include "core/Secrets.h"

#include <string>
#include <vector>

#if defined(_WIN32)
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#  include <wincrypt.h>
#endif

namespace biocad {

namespace {

constexpr char kB64[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

std::string base64Encode(const unsigned char* data, std::size_t len) {
    std::string out;
    out.reserve(((len + 2) / 3) * 4);
    std::size_t i = 0;
    for (; i + 2 < len; i += 3) {
        const unsigned n = (data[i] << 16) | (data[i + 1] << 8) | data[i + 2];
        out.push_back(kB64[(n >> 18) & 63]);
        out.push_back(kB64[(n >> 12) & 63]);
        out.push_back(kB64[(n >> 6) & 63]);
        out.push_back(kB64[n & 63]);
    }
    if (i < len) {
        unsigned n = data[i] << 16;
        const bool two = (i + 1 < len);
        if (two) n |= data[i + 1] << 8;
        out.push_back(kB64[(n >> 18) & 63]);
        out.push_back(kB64[(n >> 12) & 63]);
        out.push_back(two ? kB64[(n >> 6) & 63] : '=');
        out.push_back('=');
    }
    return out;
}

bool base64Decode(const std::string& in, std::vector<unsigned char>& out) {
    auto val = [](char c) -> int {
        if (c >= 'A' && c <= 'Z') return c - 'A';
        if (c >= 'a' && c <= 'z') return c - 'a' + 26;
        if (c >= '0' && c <= '9') return c - '0' + 52;
        if (c == '+') return 62;
        if (c == '/') return 63;
        return -1;
    };
    int buffer = 0;
    int bits = 0;
    for (char c : in) {
        if (c == '=' || c == '\n' || c == '\r') continue;
        const int v = val(c);
        if (v < 0) return false;
        buffer = (buffer << 6) | v;
        bits += 6;
        if (bits >= 8) {
            bits -= 8;
            out.push_back(static_cast<unsigned char>((buffer >> bits) & 0xFF));
        }
    }
    return true;
}

}  // namespace

#if defined(_WIN32)

Result<std::string> Secrets::protect(const std::string& plaintext) {
    DATA_BLOB in{};
    in.pbData = reinterpret_cast<BYTE*>(const_cast<char*>(plaintext.data()));
    in.cbData = static_cast<DWORD>(plaintext.size());

    DATA_BLOB out{};
    if (!CryptProtectData(&in, L"BioCAD", nullptr, nullptr, nullptr,
                          CRYPTPROTECT_UI_FORBIDDEN, &out)) {
        return Error::internal("CryptProtectData failed");
    }
    std::string encoded = base64Encode(out.pbData, out.cbData);
    LocalFree(out.pbData);
    return encoded;
}

Result<std::string> Secrets::unprotect(const std::string& base64Blob) {
    std::vector<unsigned char> raw;
    if (!base64Decode(base64Blob, raw) || raw.empty()) {
        return Error::parse("invalid base64 secret blob");
    }
    DATA_BLOB in{};
    in.pbData = raw.data();
    in.cbData = static_cast<DWORD>(raw.size());

    DATA_BLOB out{};
    if (!CryptUnprotectData(&in, nullptr, nullptr, nullptr, nullptr,
                            CRYPTPROTECT_UI_FORBIDDEN, &out)) {
        return Error::internal("CryptUnprotectData failed");
    }
    std::string plaintext(reinterpret_cast<char*>(out.pbData), out.cbData);
    LocalFree(out.pbData);
    return plaintext;
}

#else  // non-Windows fallback (dev/test only): base64, NOT encrypted.

Result<std::string> Secrets::protect(const std::string& plaintext) {
    return base64Encode(reinterpret_cast<const unsigned char*>(plaintext.data()), plaintext.size());
}

Result<std::string> Secrets::unprotect(const std::string& base64Blob) {
    std::vector<unsigned char> raw;
    if (!base64Decode(base64Blob, raw)) return Error::parse("invalid base64");
    return std::string(raw.begin(), raw.end());
}

#endif

}  // namespace biocad
