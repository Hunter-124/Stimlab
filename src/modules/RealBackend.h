// modules/RealBackend.h - the REAL analysis backend.
// Every numeric property is computed from molecular structure by the in-house
// cheminformatics engine (src/chem): SMILES -> graph -> descriptors / TPSA /
// fingerprints / functional groups. No hardcoded property values, no fakes.
#pragma once

#include <memory>

#include "contracts/Services.h"

namespace stimlab {

class RealBackend {
public:
    RealBackend();
    ~RealBackend();
    RealBackend(const RealBackend&) = delete;
    RealBackend& operator=(const RealBackend&) = delete;

    Services services();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace stimlab
