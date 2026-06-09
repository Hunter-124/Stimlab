#include "modules/Pipelines.h"

#include <cstdio>
#include <string>

#include "chem/Embed3D.h"
#include "chem/Smiles.h"
#include "data/Domain.h"
#include "modules/docking/EngineLocator.h"
#include "modules/docking/ReceptorPrep.h"

namespace stimlab::workflow {

namespace {
std::string f2(double v) {
    char b[32];
    std::snprintf(b, sizeof(b), "%.2f", v);
    return b;
}
}  // namespace

Dag buildDockingPipeline(const std::string& ligandSmiles, const std::string& targetId, Services svc) {
    // Fingerprint the EXTERNAL inputs (prepared-receptor presence + engine
    // availability) into the node params at build time. The content cache keys on
    // params, not on the live filesystem, so baking these in makes a provisioning
    // change (receptor prepared / vina downloaded) correctly invalidate the cache and
    // re-run the affected nodes - while an unchanged (ligand, target, caps) re-runs as
    // a pure cache hit.
    const bool recReady = docking::locatePreparedReceptor(targetId).ready;
    const bool vinaOk = docking::engineAvailable(docking::Engine::Vina);
    const std::string caps =
        std::string("caps=") + (recReady ? "R" : "-") + (vinaOk ? "V" : "-");

    Dag dag;

    Node ligand;
    ligand.id = "ligand_prep";
    ligand.module = "ligand_prep";
    ligand.params = ligandSmiles;
    ligand.run = [ligandSmiles](const NodeInputs&, const CancelToken&) {
        auto g = chem::parseSmiles(ligandSmiles);
        if (!g) return NodeResult::failure("could not parse SMILES");
        const chem::Conformer conf = chem::embed3D(*g);
        if (conf.empty()) return NodeResult::failure("3D embedding produced no atoms");
        return NodeResult::success("{\"atoms\":" + std::to_string(conf.size()) +
                                   ",\"heavy\":" + std::to_string(conf.heavyCount) + "}");
    };
    dag.add(std::move(ligand));

    Node receptor;
    receptor.id = "receptor_prep";
    receptor.module = "receptor_prep";
    receptor.params = targetId + "|ready=" + (recReady ? "1" : "0");
    receptor.run = [targetId](const NodeInputs&, const CancelToken&) {
        const auto rec = docking::locatePreparedReceptor(targetId);
        // A missing receptor is NOT a failure: docking degrades to the labeled
        // estimate. The node simply reports readiness for the live view.
        return NodeResult::success(std::string("{\"ready\":") + (rec.ready ? "true" : "false") +
                                   ",\"hasBox\":" + (rec.hasBox ? "true" : "false") + "}");
    };
    dag.add(std::move(receptor));

    Node dock;
    dock.id = "dock";
    dock.module = "dock";
    dock.params = ligandSmiles + "|" + targetId + "|" + caps;
    dock.deps = {"ligand_prep", "receptor_prep"};
    dock.run = [ligandSmiles, targetId, svc](const NodeInputs&, const CancelToken&) {
        if (!svc.docking) return NodeResult::failure("no docking service");
        Molecule m;
        m.id = "__wf__";
        m.name = "workflow-ligand";
        m.smiles = ligandSmiles;
        const DockJobResult d = svc.docking->dockDetailed(m, targetId);
        if (d.poses.empty()) return NodeResult::failure("docking produced no poses");
        return NodeResult::success(std::string("{\"engine\":\"") + d.engine + "\",\"real\":" +
                                   (d.real ? "true" : "false") + ",\"best\":" + f2(d.bestAffinity()) +
                                   ",\"poses\":" + std::to_string(d.poses.size()) + "}");
    };
    dag.add(std::move(dock));

    return dag;
}

}  // namespace stimlab::workflow
