// fakes/FakeBackend.h - "thick fakes": plausible, deterministic data for every
// module so the entire UI renders and is navigable before the real RDKit-backed
// science lands (Phase C). Owns the curated default library.
#pragma once

#include <memory>

#include "contracts/Services.h"

namespace stimlab {

class FakeBackend {
public:
    FakeBackend();
    ~FakeBackend();
    FakeBackend(const FakeBackend&) = delete;
    FakeBackend& operator=(const FakeBackend&) = delete;

    // Non-owning view handed to the UI.
    Services services();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace stimlab
