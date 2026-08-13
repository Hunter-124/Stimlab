// core/Secrets.h - encrypt small secrets (API keys) at rest via Windows DPAPI.
// Output is base64 of the DPAPI blob; safe to store in config.json / SQLite.
#pragma once

#include <string>

#include "core/Error.h"

namespace biocad {

class Secrets {
public:
    // Encrypt under the current user account. Returns base64 text.
    static Result<std::string> protect(const std::string& plaintext);

    // Reverse of protect(): base64 DPAPI blob -> plaintext.
    static Result<std::string> unprotect(const std::string& base64Blob);
};

}  // namespace biocad
