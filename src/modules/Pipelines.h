// modules/Pipelines.h - build re-runnable analysis pipelines as workflow DAGs (WP-D).
//
// Rebuilds the legacy prep->dock pipeline as a content-cached, cancellable DAG on the
// WP-D engine (src/workflow). buildDockingPipeline() returns a 3-node graph:
//
//   ligand_prep   : SMILES -> in-house 3D embed (chem::embed3D); emits an atom count
//   receptor_prep : locate the target's prepared receptor PDBQT (cache-only, no net)
//   dock          : dock the ligand into the target via the docking service
//
// Each node's run-fn returns a compact JSON status string that is BOTH the payload
// passed downstream and the cache value, so an unchanged (ligand, target) re-runs as
// a pure cache hit. The nodes call only analysis services; they do not implement
// synthesis, route, or manufacturability workflows.
#pragma once

#include <string>

#include "contracts/Services.h"
#include "workflow/Dag.h"

namespace biocad::workflow {

// Build the prep->dock pipeline DAG for one ligand SMILES against one CNS target id
// (e.g. "DAT"). `svc` is copied (it is a bundle of non-owning service pointers) into
// the node functions, so it must outlive any DagExecutor::run of the returned graph.
Dag buildDockingPipeline(const std::string& ligandSmiles, const std::string& targetId,
                         Services svc);

}  // namespace biocad::workflow
