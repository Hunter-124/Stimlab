#include <catch2/catch_test_macros.hpp>

#include <string>

#include "modules/Pipelines.h"
#include "modules/RealBackend.h"
#include "workflow/Dag.h"
#include "workflow/JobSystem.h"

using namespace stimlab;
using namespace stimlab::workflow;

TEST_CASE("Docking pipeline builds a 3-node prep->dock DAG and runs", "[workflow][pipeline]") {
    RealBackend backend;
    Services svc = backend.services();

    // Dock into a target with no prepared receptor so the dock node deterministically
    // uses the labeled descriptor estimate (no engine subprocess) - hermetic + fast.
    Dag dag = buildDockingPipeline("CC(N)Cc1ccccc1", "__wf_unprepared__", svc);
    REQUIRE(dag.validate().empty());
    REQUIRE(dag.nodes().size() == 3);
    REQUIRE(dag.find("dock") != nullptr);
    REQUIRE(dag.find("dock")->deps.size() == 2);  // ligand_prep + receptor_prep

    JobSystem jobs(4);
    MemoryNodeCache cache;
    DagExecutor exec(jobs, cache);
    const auto r = exec.run(dag, CancelToken{});

    REQUIRE(r.ok);
    REQUIRE(r.ran == 3);
    REQUIRE(r.statusOf("ligand_prep") == NodeStatus::Done);
    REQUIRE(r.statusOf("dock") == NodeStatus::Done);
    // The ligand node reports a non-empty atom count; the dock node reports a labeled
    // (real=false) estimate result with poses.
    REQUIRE(r.output.at("ligand_prep").find("\"atoms\":") != std::string::npos);
    REQUIRE(r.output.at("dock").find("\"real\":false") != std::string::npos);
    REQUIRE(r.output.at("dock").find("\"poses\":") != std::string::npos);
}

TEST_CASE("Re-running the docking pipeline is a full cache hit", "[workflow][pipeline]") {
    RealBackend backend;
    Services svc = backend.services();
    JobSystem jobs(4);
    MemoryNodeCache cache;
    DagExecutor exec(jobs, cache);

    auto build = [&] { return buildDockingPipeline("CNC(C)Cc1ccccc1", "__wf_unprepared__", svc); };
    const auto r1 = exec.run(build(), CancelToken{});
    REQUIRE(r1.ran == 3);

    const auto r2 = exec.run(build(), CancelToken{});  // identical inputs + shared cache
    REQUIRE(r2.ok);
    REQUIRE(r2.cached == 3);
    REQUIRE(r2.ran == 0);
    REQUIRE(r2.output.at("dock") == r1.output.at("dock"));
}
