#include "modules/SimulationModule.h"

#include <algorithm>

#include "sim/Control.h"
#include "sim/Enrichment.h"
#include "sim/Kinetics.h"
#include "sim/Network.h"
#include "sim/Sbml.h"
#include "sim/Solvers.h"

#if BIOCAD_ENABLE_FBA
#include "sim/Flux.h"
#endif

namespace biocad {

NetworkSpec RealSimulation::analyze(const NetworkSpec& network) const {
    return sim::analyze(network);
}

TimeCourse RealSimulation::integrate(const NetworkSpec& network, double horizon, double relTol,
                                     double absTol, const std::string& method) const {
    sim::IntegrationOptions o;
    // A caller that asks for a horizon of zero or a nonsensical tolerance gets the
    // documented default rather than a divide-by-zero deep in the step controller.
    o.horizon = horizon > 0 ? horizon : 10.0;
    o.relativeTolerance = relTol > 0 ? relTol : 1e-8;
    o.absoluteTolerance = absTol > 0 ? absTol : 1e-12;
    o.method = method == "rk4" ? "rk4" : "rosenbrock";
    o.outputPoints = 201;
    // The structural analysis has to have run for the conserved-quantity audit to
    // exist, so it is run here when the caller did not.
    const NetworkSpec analyzed =
        network.conservationLaws.empty() ? sim::analyze(network) : network;
    return sim::integrate(analyzed, o);
}

StochasticEnsemble RealSimulation::stochastic(const NetworkSpec& network, double horizon,
                                              int replicates, std::uint64_t seed,
                                              bool tauLeap) const {
    sim::StochasticOptions o;
    o.horizon = horizon > 0 ? horizon : 10.0;
    o.replicates = std::clamp(replicates, 1, 20000);
    o.seed = seed;
    o.tauLeap = tauLeap;
    o.outputPoints = 101;
    return sim::stochastic(network, o);
}

KineticsFit RealSimulation::arrhenius(const std::vector<double>& temperaturesK,
                                      const std::vector<double>& rateConstants) const {
    return sim::arrhenius(temperaturesK, rateConstants);
}

PhRateProfile RealSimulation::phRate(const std::vector<double>& pHValues,
                                     const std::vector<double>& rateConstants) const {
    return sim::phRate(pHValues, rateConstants);
}

ControlAnalysis RealSimulation::controlAnalysis(const NetworkSpec& network,
                                                double horizon) const {
    return sim::controlAnalysis(network, horizon > 0 ? horizon : 100.0);
}

std::optional<NetworkSpec> RealSimulation::importSbml(const std::string& xml,
                                                      std::string* error) const {
    return sim::importSbml(xml, error);
}

std::string RealSimulation::exportSbml(const NetworkSpec& network) const {
    return sim::exportSbml(network);
}

// ---------------------------------------------------------------------------

EnrichmentReport RealEnrichment::enrich(const std::vector<std::string>& query,
                                        const std::vector<std::string>& background,
                                        const std::string& gmtPack) const {
    const sim::GeneSetPack pack =
        sim::loadGeneSetPack(gmtPack.empty() ? "reactome-human.gmt" : gmtPack);
    EnrichmentReport out = sim::enrich(query, background, pack);
    // A missing pack is a load failure, not an empty result: say so where the panel
    // will show it.
    for (const std::string& w : pack.warnings) out.warnings.push_back(w);
    if (pack.ids.empty())
        out.warnings.push_back("the gene-set pack '" + gmtPack +
                               "' loaded zero pathways, so there was nothing to test");
    return out;
}

GraphMetrics RealEnrichment::graph(const std::vector<NetworkEdge>& edges) const {
    return sim::graph(edges);
}

// ---------------------------------------------------------------------------

#if BIOCAD_ENABLE_FBA
FluxSolution RealFlux::balance(const NetworkSpec& network) const { return sim::balance(network); }

FluxSolution RealFlux::fba(const NetworkSpec& network, const std::string& objectiveReactionId,
                           const std::vector<FluxBound>& bounds) const {
    return sim::fba(network, objectiveReactionId, bounds);
}

std::vector<FluxRange> RealFlux::fva(const NetworkSpec& network,
                                     const std::string& objectiveReactionId,
                                     const std::vector<FluxBound>& bounds,
                                     double objectiveFraction) const {
    return sim::fva(network, objectiveReactionId, bounds,
                    objectiveFraction > 0 && objectiveFraction <= 1.0 ? objectiveFraction : 1.0);
}

FluxSolution RealFlux::parsimonious(const NetworkSpec& network,
                                    const std::string& objectiveReactionId,
                                    const std::vector<FluxBound>& bounds) const {
    return sim::parsimonious(network, objectiveReactionId, bounds);
}

std::vector<FluxRange> RealFlux::deletions(const NetworkSpec& network,
                                           const std::string& objectiveReactionId,
                                           const std::vector<FluxBound>& bounds,
                                           int order) const {
    return sim::deletions(network, objectiveReactionId, bounds, order);
}
#endif

}  // namespace biocad
