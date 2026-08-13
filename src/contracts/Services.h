// contracts/Services.h - a non-owning bundle of module interfaces.
// The UI receives one Services view and codes only against the interfaces, so a
// module can be replaced without touching a panel. There is exactly one
// implementation of each: no fake backend exists, because a test suite that
// validates a double while the product ships the original proves nothing.
#pragma once

#include "contracts/IModules.h"

namespace biocad {

struct Services {
    ILibrary*           library    = nullptr;
    IStabilityModule*   stability  = nullptr;
    IAdmetModule*       admet      = nullptr;
    IAbsorptionModule*  absorption = nullptr;
    ISimilarityModule*  similarity = nullptr;
    ILegalModule*       legal      = nullptr;
    IDockingModule*     docking    = nullptr;
    IRunStore*          runs       = nullptr;
    IPharmacodynamicsModule* pharmacodynamics = nullptr;
    ISequenceModule*    sequence   = nullptr;
    IStructureModule*   structure  = nullptr;
    IMetabolismFactsModule* metabolismFacts = nullptr;
    IAlertsModule*      alerts     = nullptr;
    IIonizationModule*  ionization = nullptr;
    INucleicAcidModule* nucleicAcid = nullptr;
    IAssayModule*       assay      = nullptr;
    IBiologicsModule*   biologics  = nullptr;
    IPopulationPkModule* populationPk = nullptr;
    // Phase 14. `flux` stays a nullptr in a build without BIOCAD_ENABLE_FBA and is
    // deliberately NOT required by valid(): the constraint-based path is an optional
    // feature, so requiring it would make every default build invalid.
    ISimulationModule*  simulation = nullptr;
    IEnrichmentModule*  enrichment = nullptr;
    IFluxModule*        flux       = nullptr;

    [[nodiscard]] bool valid() const {
        return library && stability && admet && absorption &&
               similarity && legal && docking && runs && pharmacodynamics &&
               sequence && structure && metabolismFacts && alerts && ionization &&
               nucleicAcid && assay && biologics && populationPk && simulation && enrichment;
    }
};

}  // namespace biocad
