#include "pkpd/Pharmacodynamics.h"

#include "pkpd/Fits.h"
#include "pkpd/PkEngine.h"

namespace biocad::pkpd {

CurveFit RealPharmacodynamics::fitFourParameterLogistic(
    const std::vector<DoseResponsePoint>& points) const {
    // Unweighted by default: 1/y^2 weighting is a deliberate choice about the assay's
    // error model, not something a module should apply behind the caller's back.
    return pkpd::fitFourParameterLogistic(points, false);
}

Quantity RealPharmacodynamics::kiFromIc50(const ChengPrusoffInput& input) const {
    return pkpd::kiFromIc50(input);
}

SchildResult RealPharmacodynamics::schild(const std::vector<SchildPoint>& points) const {
    return pkpd::schild(points);
}

PkProfile RealPharmacodynamics::simulate(const PkModelSpec& spec,
                                         const DoseRegimen& regimen) const {
    return pkpd::simulate(spec, regimen);
}

OccupancyCurve RealPharmacodynamics::occupancy(const PkProfile& profile,
                                               const Quantity& kd) const {
    return pkpd::occupancy(profile, kd);
}

}  // namespace biocad::pkpd
