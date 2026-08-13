#include "modules/PopulationPkModule.h"

#include "sim/Ddi.h"
#include "sim/Nca.h"
#include "sim/Population.h"

namespace biocad {

PopulationProfile RealPopulationPk::simulate(const PkModelSpec& model,
                                             const DoseRegimen& regimen,
                                             const VariabilitySpec& variability) const {
    return sim::simulatePopulation(model, regimen, variability);
}

NcaResult RealPopulationPk::nca(const ConcentrationSeries& observed) const {
    return sim::noncompartmental(observed);
}

InteractionReport RealPopulationPk::interaction(const PerpetratorSpec& perpetrator,
                                                const VictimSpec& victim) const {
    return sim::interaction(perpetrator, victim);
}

EnzymeTimeCourse RealPopulationPk::enzymeTimeCourse(const PerpetratorSpec& perpetrator,
                                                    double horizonH) const {
    return sim::enzymeTimeCourse(perpetrator, horizonH);
}

ImpairmentScenario RealPopulationPk::impairment(const VictimSpec& victim,
                                                double renalFunctionRatio,
                                                double hepaticClintRatio) const {
    return sim::impairment(victim, renalFunctionRatio, hepaticClintRatio);
}

}  // namespace biocad
