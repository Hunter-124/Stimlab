// contracts/Services.h - a non-owning bundle of module interfaces.
// The UI receives one Services view and never knows whether it is wired to fakes
// or real RDKit-backed implementations.
#pragma once

#include "contracts/IModules.h"

namespace stimlab {

struct Services {
    ILibrary*           library    = nullptr;
    IStabilityModule*   stability  = nullptr;
    IAdmetModule*       admet      = nullptr;
    IAbsorptionModule*  absorption = nullptr;
    ISimilarityModule*  similarity = nullptr;
    ILegalModule*       legal      = nullptr;
    IDockingModule*     docking    = nullptr;
    IRunStore*          runs       = nullptr;

    [[nodiscard]] bool valid() const {
        return library && stability && admet && absorption &&
               similarity && legal && docking && runs;
    }
};

}  // namespace stimlab
