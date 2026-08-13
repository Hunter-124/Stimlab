// pkpd/Pharmacodynamics.h - the real IPharmacodynamicsModule.
//
// This class is a thin composition seam: the curve fits live in pkpd/Fits.h and the
// exposure simulation in pkpd/PkEngine.h, so the same maths is reachable from the
// tests and the agent tools without going through the Services indirection.
#pragma once

#include "contracts/IModules.h"

namespace biocad::pkpd {

class RealPharmacodynamics final : public IPharmacodynamicsModule {
public:
    CurveFit fitFourParameterLogistic(const std::vector<DoseResponsePoint>&) const override;
    Quantity kiFromIc50(const ChengPrusoffInput&) const override;
    SchildResult schild(const std::vector<SchildPoint>&) const override;
    PkProfile simulate(const PkModelSpec&, const DoseRegimen&) const override;
    OccupancyCurve occupancy(const PkProfile&, const Quantity& kd) const override;
};

}  // namespace biocad::pkpd
