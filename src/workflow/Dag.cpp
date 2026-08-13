#include "workflow/Dag.h"

#include <algorithm>
#include <condition_variable>
#include <filesystem>
#include <fstream>
#include <future>
#include <queue>
#include <sstream>
#include <system_error>
#include <unordered_map>
#include <unordered_set>

#include "core/AppPaths.h"
#include "core/Hash.h"

namespace biocad::workflow {

namespace fs = std::filesystem;

const char* toString(NodeStatus s) {
    switch (s) {
        case NodeStatus::Pending:   return "pending";
        case NodeStatus::Running:   return "running";
        case NodeStatus::Cached:    return "cached";
        case NodeStatus::Done:      return "done";
        case NodeStatus::Failed:    return "failed";
        case NodeStatus::Cancelled: return "cancelled";
        case NodeStatus::Skipped:   return "skipped";
    }
    return "?";
}

// --------------------------------------------------------------------- Dag
Dag& Dag::add(Node n) {
    nodes_.push_back(std::move(n));
    return *this;
}

const Node* Dag::find(const std::string& id) const {
    for (const auto& n : nodes_)
        if (n.id == id) return &n;
    return nullptr;
}

std::string Dag::validate() const {
    std::unordered_set<std::string> ids;
    for (const auto& n : nodes_) {
        if (n.id.empty()) return "a node has an empty id";
        if (!ids.insert(n.id).second) return "duplicate node id: " + n.id;
    }
    for (const auto& n : nodes_)
        for (const auto& d : n.deps) {
            if (d == n.id) return "node depends on itself: " + n.id;
            if (!ids.count(d)) return "node " + n.id + " depends on missing node: " + d;
        }
    if (!nodes_.empty() && topoOrder().empty()) return "the graph has a cycle";
    return {};
}

std::vector<std::string> Dag::topoOrder() const {
    std::unordered_map<std::string, int> indeg;
    std::unordered_map<std::string, std::vector<std::string>> deps;  // dep -> dependents
    for (const auto& n : nodes_) indeg.emplace(n.id, 0);
    for (const auto& n : nodes_) {
        for (const auto& d : n.deps) {
            if (!indeg.count(d)) return {};  // missing dep -> not orderable
            ++indeg[n.id];
            deps[d].push_back(n.id);
        }
    }
    // Deterministic order: seed ready set in declaration order.
    std::vector<std::string> order;
    std::queue<std::string> ready;
    for (const auto& n : nodes_)
        if (indeg[n.id] == 0) ready.push(n.id);
    while (!ready.empty()) {
        const std::string id = ready.front();
        ready.pop();
        order.push_back(id);
        for (const auto& dep : deps[id])
            if (--indeg[dep] == 0) ready.push(dep);
    }
    if (order.size() != nodes_.size()) return {};  // cycle
    return order;
}

// --------------------------------------------------------------- caches
std::optional<std::string> MemoryNodeCache::get(const std::string& key) {
    std::lock_guard<std::mutex> lk(mu_);
    auto it = store_.find(key);
    if (it == store_.end()) return std::nullopt;
    return it->second;
}

void MemoryNodeCache::put(const std::string& key, const std::string& output) {
    std::lock_guard<std::mutex> lk(mu_);
    store_[key] = output;
}

std::size_t MemoryNodeCache::size() const {
    std::lock_guard<std::mutex> lk(mu_);
    return store_.size();
}

DiskNodeCache::DiskNodeCache(std::string subdir) {
    std::error_code ec;
    const fs::path dir = AppPaths::instance().cache() / subdir;
    fs::create_directories(dir, ec);
    dir_ = dir.string();
}

std::optional<std::string> DiskNodeCache::get(const std::string& key) {
    std::error_code ec;
    const fs::path p = fs::path(dir_) / (key + ".txt");
    if (!fs::exists(p, ec)) return std::nullopt;
    std::ifstream in(p, std::ios::binary);
    if (!in) return std::nullopt;
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

void DiskNodeCache::put(const std::string& key, const std::string& output) {
    const fs::path p = fs::path(dir_) / (key + ".txt");
    std::ofstream out(p, std::ios::binary);
    out << output;
}

// --------------------------------------------------------- DagRunResult
NodeStatus DagRunResult::statusOf(const std::string& id) const {
    auto it = status.find(id);
    return it == status.end() ? NodeStatus::Pending : it->second;
}

// ----------------------------------------------------------- DagExecutor
DagRunResult DagExecutor::run(const Dag& dag, const CancelToken& cancel, ProgressFn onProgress) {
    DagRunResult res;
    if (const std::string verr = dag.validate(); !verr.empty()) {
        res.error = verr;
        return res;
    }

    const auto& nodes = dag.nodes();
    const int N = static_cast<int>(nodes.size());
    if (N == 0) {
        res.ok = true;
        return res;
    }

    std::unordered_map<std::string, int> idx;
    for (int i = 0; i < N; ++i) idx[nodes[i].id] = i;

    // Transitive cache keys, computed in topological order so dep keys are ready.
    std::unordered_map<std::string, std::string> key;
    for (const std::string& id : dag.topoOrder()) {
        const Node& n = nodes[idx[id]];
        std::string mat = n.module;
        mat += '\x1f';
        mat += std::to_string(n.version);
        mat += '\x1f';
        mat += n.params;
        mat += '\x1f';
        std::vector<std::string> ds = n.deps;
        std::sort(ds.begin(), ds.end());
        for (const auto& d : ds) {
            mat += d;
            mat += '=';
            mat += key[d];
            mat += '\x1f';
        }
        key[id] = hashHex(mat);
    }
    for (const auto& [id, k] : key) res.cacheKey[id] = k;

    std::vector<std::vector<int>> dependents(N);
    std::vector<int> indeg(N, 0);
    for (int i = 0; i < N; ++i) {
        indeg[i] = static_cast<int>(nodes[i].deps.size());
        for (const auto& d : nodes[i].deps) dependents[idx[d]].push_back(i);
    }

    std::vector<NodeStatus> status(N, NodeStatus::Pending);
    std::vector<std::string> output(N);

    std::mutex m;
    std::condition_variable cv;
    int finished = 0;
    std::queue<int> ready;
    std::vector<NodeProgress> pending;        // progress events, drained by the scheduler
    std::vector<std::future<void>> futs;      // one per scheduled worker (joined before return)
    for (int i = 0; i < N; ++i)
        if (indeg[i] == 0) ready.push(i);

    // MUST hold m. Marks node i finished, releases dependents, records a progress event.
    auto finalize = [&](int i, NodeStatus s, std::string detail) {
        status[i] = s;
        ++finished;
        pending.push_back({nodes[i].id, s, std::move(detail)});
        for (int j : dependents[i])
            if (--indeg[j] == 0) ready.push(j);
    };

    std::unique_lock<std::mutex> lk(m);
    while (finished < N) {
        while (!ready.empty()) {
            const int i = ready.front();
            ready.pop();

            bool depBad = false;
            for (const auto& d : nodes[i].deps) {
                const NodeStatus ds = status[idx[d]];
                if (ds == NodeStatus::Failed || ds == NodeStatus::Skipped ||
                    ds == NodeStatus::Cancelled)
                    depBad = true;
            }
            if (cancel.cancelled()) {
                finalize(i, NodeStatus::Cancelled, "run cancelled");
                continue;
            }
            if (depBad) {
                finalize(i, NodeStatus::Skipped, "upstream failed/cancelled");
                continue;
            }
            if (auto hit = cache_.get(key[nodes[i].id])) {
                output[i] = *hit;
                finalize(i, NodeStatus::Cached, "cache hit");
                continue;
            }
            if (!nodes[i].run) {  // pure data node: empty success, still cached
                output[i].clear();
                cache_.put(key[nodes[i].id], output[i]);
                finalize(i, NodeStatus::Done, "no-op");
                continue;
            }

            // Schedule the node on the pool. Status flips to Running now; completion
            // is finalized by the worker under the same lock.
            status[i] = NodeStatus::Running;
            pending.push_back({nodes[i].id, NodeStatus::Running, "running"});
            NodeInputs inputs;
            for (const auto& d : nodes[i].deps) inputs[d] = output[idx[d]];
            const std::string k = key[nodes[i].id];
            futs.push_back(jobs_.submit([&, i, inputs, k] {
                // A throwing node fn must still finalize the node (else `finished`
                // never reaches N and run() hangs) - treat it as a failure.
                NodeResult r;
                try {
                    r = nodes[i].run(inputs, cancel);
                } catch (const std::exception& e) {
                    r = NodeResult::failure(std::string("exception: ") + e.what());
                } catch (...) {
                    r = NodeResult::failure("unknown exception");
                }
                std::lock_guard<std::mutex> lg(m);
                if (r.ok) {
                    output[i] = r.output;
                    cache_.put(k, r.output);
                    finalize(i, NodeStatus::Done, "ok");
                } else {
                    finalize(i, NodeStatus::Failed, r.error.empty() ? "failed" : r.error);
                }
                cv.notify_all();
            }));
        }

        // Fire progress callbacks outside the lock (single-threaded here, so the
        // callback never races and may safely touch UI state via the caller).
        if (onProgress && !pending.empty()) {
            std::vector<NodeProgress> batch;
            batch.swap(pending);
            lk.unlock();
            for (const auto& p : batch) onProgress(p);
            lk.lock();
            continue;  // re-check ready (a worker may have enqueued more)
        }

        if (finished < N) cv.wait(lk);
    }
    // Drain any final progress events.
    if (onProgress && !pending.empty()) {
        std::vector<NodeProgress> batch;
        batch.swap(pending);
        lk.unlock();
        for (const auto& p : batch) onProgress(p);
        lk.lock();
    }

    // CRITICAL: every scheduled worker must FULLY return before this frame's locals
    // (m, cv, status, ...) are destroyed. The last worker calls cv.notify_all() after
    // finalize, which can wake the scheduler to exit here; without this join the
    // scheduler could destroy cv/m while that worker is still unwinding the condvar
    // machinery (a notify-then-destroy use-after-free). A future completes only after
    // its packaged_task - and thus the whole worker lambda - has returned.
    if (lk.owns_lock()) lk.unlock();
    for (auto& f : futs) f.get();

    for (int i = 0; i < N; ++i) {
        res.status[nodes[i].id] = status[i];
        res.output[nodes[i].id] = output[i];
        switch (status[i]) {
            case NodeStatus::Done:      ++res.ran; break;
            case NodeStatus::Cached:    ++res.cached; break;
            case NodeStatus::Failed:    ++res.failed; break;
            case NodeStatus::Cancelled:
            case NodeStatus::Skipped:   ++res.cancelled; break;
            default: break;
        }
    }
    res.ok = (res.failed == 0 && res.cancelled == 0);
    return res;
}

}  // namespace biocad::workflow
