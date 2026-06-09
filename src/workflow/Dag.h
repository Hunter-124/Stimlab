// workflow/Dag.h - a re-runnable, content-cached, cancellable workflow DAG (WP-D).
//
// A Dag is a set of named Nodes, each with declared upstream deps and a pure-ish run
// function. The DagExecutor schedules independent nodes in parallel on a JobSystem
// and gives three properties the legacy prep->dock pipelines needed:
//
//   * CONTENT-ADDRESSED CACHING / RESUME. Each node has a cache key
//     hash(module, version, params, {dep -> dep cache key}). The key is transitive,
//     so re-running an unchanged graph is a full cache hit (instant resume); changing
//     one node's params/version re-runs that node and everything downstream of it,
//     and nothing else. Cached outputs live behind an INodeCache (in-memory or
//     %APPDATA%/StimLab/cache on disk).
//   * COOPERATIVE CANCEL. A CancelToken stops scheduling new work; in-flight nodes
//     observe it and bail; un-started nodes are marked Cancelled.
//   * PROGRESS. An optional callback fires on every node status transition for a live
//     DAG view.
//
// A node function returns a NodeResult whose `output` string is BOTH the value passed
// to downstream nodes (as NodeInputs) and the cached payload. Keep node functions
// deterministic in their (inputs, params) so caching is sound.
//
// SAFETY SCOPE: the DAG orchestrates analysis modules (ligand prep, receptor prep,
// docking, ADMET screening). It contains no synthesis/route/manufacturability logic;
// node functions are supplied by callers and must honor the same boundary.
#pragma once

#include <functional>
#include <map>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

#include "workflow/JobSystem.h"

namespace stimlab::workflow {

enum class NodeStatus { Pending, Running, Cached, Done, Failed, Cancelled, Skipped };
const char* toString(NodeStatus s);

// The result of running one node. `output` is the payload handed to dependents and
// stored in the cache; on failure it is empty and `error` explains why.
struct NodeResult {
    bool        ok = false;
    std::string output;
    std::string error;

    static NodeResult success(std::string out) { return {true, std::move(out), {}}; }
    static NodeResult failure(std::string err) { return {false, {}, std::move(err)}; }
};

// Resolved upstream outputs handed to a node: dep id -> that dep's output string.
using NodeInputs = std::map<std::string, std::string>;
using NodeFn = std::function<NodeResult(const NodeInputs&, const CancelToken&)>;

// A workflow node. `module` + `version` + `params` are the cache-key inputs that are
// intrinsic to this node; `deps` pull the rest of the key transitively.
struct Node {
    std::string              id;          // unique within the Dag
    std::string              module;      // logical module (cache-key component)
    int                      version = 1; // bump to invalidate this node's cache
    std::string              params;      // node parameters (cache-key component)
    std::vector<std::string> deps;        // upstream node ids
    NodeFn                   run;         // may be empty for a pure data node
};

// A DAG of nodes. Build with add(); validate()/topoOrder() before executing.
class Dag {
public:
    Dag& add(Node n);
    [[nodiscard]] const std::vector<Node>& nodes() const { return nodes_; }
    [[nodiscard]] const Node* find(const std::string& id) const;

    // "" when the graph is well-formed; otherwise a human reason (duplicate id,
    // missing dependency, or cycle).
    [[nodiscard]] std::string validate() const;
    // Node ids in a dependency-respecting order; empty if the graph is cyclic.
    [[nodiscard]] std::vector<std::string> topoOrder() const;

private:
    std::vector<Node> nodes_;
};

// Pluggable content-addressed cache of node outputs keyed by transitive cache key.
class INodeCache {
public:
    virtual ~INodeCache() = default;
    virtual std::optional<std::string> get(const std::string& key) = 0;
    virtual void                       put(const std::string& key, const std::string& output) = 0;
};

// Process-lifetime in-memory cache (the default; demonstrates resume within a run).
class MemoryNodeCache final : public INodeCache {
public:
    std::optional<std::string> get(const std::string& key) override;
    void                       put(const std::string& key, const std::string& output) override;
    [[nodiscard]] std::size_t  size() const;

private:
    std::map<std::string, std::string> store_;
    mutable std::mutex                 mu_;
};

// Disk-backed cache under %APPDATA%/StimLab/cache/<sub> for cross-session resume.
class DiskNodeCache final : public INodeCache {
public:
    explicit DiskNodeCache(std::string subdir = "dag");
    std::optional<std::string> get(const std::string& key) override;
    void                       put(const std::string& key, const std::string& output) override;

private:
    std::string dir_;  // absolute cache directory
};

struct NodeProgress {
    std::string id;
    NodeStatus  status = NodeStatus::Pending;
    std::string detail;
};
using ProgressFn = std::function<void(const NodeProgress&)>;

// The outcome of one DAG run. `output`/`status`/`cacheKey` are keyed by node id.
struct DagRunResult {
    bool                              ok = false;
    std::map<std::string, NodeStatus> status;
    std::map<std::string, std::string> output;
    std::map<std::string, std::string> cacheKey;
    int                               ran = 0;        // nodes actually executed
    int                               cached = 0;      // nodes served from cache
    int                               failed = 0;
    int                               cancelled = 0;   // incl. skipped-after-failure
    std::string                       error;           // graph-level error (invalid DAG)

    [[nodiscard]] NodeStatus statusOf(const std::string& id) const;
};

// Schedules a Dag on a JobSystem with caching/cancel/progress. One executor can run
// many Dags; results accumulate in the shared cache so re-runs resume.
class DagExecutor {
public:
    DagExecutor(JobSystem& jobs, INodeCache& cache) : jobs_(jobs), cache_(cache) {}

    DagRunResult run(const Dag& dag, const CancelToken& cancel, ProgressFn onProgress = {});

private:
    JobSystem&  jobs_;
    INodeCache& cache_;
};

}  // namespace stimlab::workflow
