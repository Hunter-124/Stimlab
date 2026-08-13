// modules/SimulationModule.h - the Phase 14 adapters between the src/sim engines and
// the rest of BioCAD.
//
// Three interfaces, three implementations, one engine set. RealSimulation owns the
// deterministic and stochastic solvers, the chemical-kinetics fits and the SBML
// reader/writer; RealEnrichment owns the over-representation test and the graph
// metrics; RealFlux owns the constraint-based flux path and exists only when
// BIOCAD_ENABLE_FBA is defined, because a flux panel with no reconstruction loaded
// invites a reader to treat a computed growth rate as an organism measurement.
//
// The adapters are thin on purpose: every refusal (an unsupported SBML construct, a
// saturable rate law handed to the SSA, an unbalanced reaction handed to FBA) is
// decided in src/sim and travels out through the DTO, so a panel and an agent tool
// see exactly the same answer.
#pragma once

#include <optional>
#include <string>
#include <vector>

#include "contracts/IModules.h"

namespace biocad {

class RealSimulation final : public ISimulationModule {
public:
    NetworkSpec analyze(const NetworkSpec& network) const override;
    TimeCourse  integrate(const NetworkSpec& network, double horizon, double relTol, double absTol,
                          const std::string& method) const override;
    StochasticEnsemble stochastic(const NetworkSpec& network, double horizon, int replicates,
                                  std::uint64_t seed, bool tauLeap) const override;
    KineticsFit   arrhenius(const std::vector<double>& temperaturesK,
                            const std::vector<double>& rateConstants) const override;
    PhRateProfile phRate(const std::vector<double>& pHValues,
                         const std::vector<double>& rateConstants) const override;
    ControlAnalysis controlAnalysis(const NetworkSpec& network, double horizon) const override;
    std::optional<NetworkSpec> importSbml(const std::string& xml,
                                          std::string* error) const override;
    std::string exportSbml(const NetworkSpec& network) const override;
};

class RealEnrichment final : public IEnrichmentModule {
public:
    EnrichmentReport enrich(const std::vector<std::string>& query,
                            const std::vector<std::string>& background,
                            const std::string& gmtPack) const override;
    GraphMetrics graph(const std::vector<NetworkEdge>& edges) const override;
};

#if BIOCAD_ENABLE_FBA
class RealFlux final : public IFluxModule {
public:
    FluxSolution balance(const NetworkSpec& network) const override;
    FluxSolution fba(const NetworkSpec& network, const std::string& objectiveReactionId,
                     const std::vector<FluxBound>& bounds) const override;
    std::vector<FluxRange> fva(const NetworkSpec& network, const std::string& objectiveReactionId,
                               const std::vector<FluxBound>& bounds,
                               double objectiveFraction) const override;
    FluxSolution parsimonious(const NetworkSpec& network, const std::string& objectiveReactionId,
                              const std::vector<FluxBound>& bounds) const override;
    std::vector<FluxRange> deletions(const NetworkSpec& network,
                                      const std::string& objectiveReactionId,
                                      const std::vector<FluxBound>& bounds,
                                      int order) const override;
};
#endif

}  // namespace biocad
