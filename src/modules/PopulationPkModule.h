// modules/PopulationPkModule.h - the Phase 13 adapter between the sim:: engines and
// the rest of BioCAD.
//
// The engines in src/sim know nothing about services, panels or tools; this is the
// single place that binds them to IPopulationPkModule. It holds no state, because a
// population simulation's only state is its seed and that lives in the request.
#pragma once

#include "contracts/IModules.h"

namespace biocad {

class RealPopulationPk final : public IPopulationPkModule {
public:
    PopulationProfile simulate(const PkModelSpec& model, const DoseRegimen& regimen,
                               const VariabilitySpec& variability) const override;
    NcaResult         nca(const ConcentrationSeries& observed) const override;
    InteractionReport interaction(const PerpetratorSpec& perpetrator,
                                  const VictimSpec& victim) const override;
    EnzymeTimeCourse  enzymeTimeCourse(const PerpetratorSpec& perpetrator,
                                       double horizonH) const override;
    ImpairmentScenario impairment(const VictimSpec& victim, double renalFunctionRatio,
                                  double hepaticClintRatio) const override;
};

}  // namespace biocad
