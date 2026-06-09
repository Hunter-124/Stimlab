#include <catch2/catch_test_macros.hpp>

#include <atomic>
#include <chrono>
#include <map>
#include <mutex>
#include <string>
#include <thread>

#include "workflow/Dag.h"
#include "workflow/JobSystem.h"

using namespace stimlab::workflow;

namespace {

// A shared run-counter so tests can assert which nodes actually executed (vs were
// served from cache). Keyed by node id.
struct RunLog {
    std::mutex mu;
    std::map<std::string, int> counts;
    int count(const std::string& id) {
        std::lock_guard<std::mutex> lk(mu);
        return counts[id];
    }
    void bump(const std::string& id) {
        std::lock_guard<std::mutex> lk(mu);
        ++counts[id];
    }
};

// Build a node whose output is "<concat of dep outputs>+<id>" and which records that
// it ran. `log` and `module`/`version`/`params` let tests probe caching behaviour.
Node makeNode(const std::string& id, std::vector<std::string> deps, RunLog& log,
              std::string module, int version, std::string params) {
    Node n;
    n.id = id;
    n.module = std::move(module);
    n.version = version;
    n.params = std::move(params);
    n.deps = std::move(deps);
    n.run = [id, &log](const NodeInputs& in, const CancelToken&) {
        log.bump(id);
        std::string out;
        for (const auto& [k, v] : in) out += v + ",";
        out += id;
        return NodeResult::success(out);
    };
    return n;
}

}  // namespace

TEST_CASE("JobSystem runs submitted tasks and satisfies their futures", "[workflow][jobs]") {
    JobSystem jobs(4);
    REQUIRE(jobs.threadCount() == 4);
    std::atomic<int> n{0};
    std::vector<std::future<void>> fs;
    for (int i = 0; i < 200; ++i) fs.push_back(jobs.submit([&n] { n.fetch_add(1); }));
    for (auto& f : fs) f.get();
    REQUIRE(n.load() == 200);
}

TEST_CASE("Invalid DAGs are rejected", "[workflow][dag]") {
    RunLog log;
    {
        Dag d;
        d.add(makeNode("a", {}, log, "A", 1, ""));
        d.add(makeNode("a", {}, log, "A", 1, ""));  // duplicate id
        REQUIRE(d.validate().find("duplicate") != std::string::npos);
    }
    {
        Dag d;
        d.add(makeNode("a", {"ghost"}, log, "A", 1, ""));  // missing dep
        REQUIRE(d.validate().find("missing") != std::string::npos);
    }
    {
        Dag d;  // cycle a->b->a
        Node a = makeNode("a", {"b"}, log, "A", 1, "");
        Node b = makeNode("b", {"a"}, log, "B", 1, "");
        d.add(a).add(b);
        REQUIRE(d.validate().find("cycle") != std::string::npos);
        JobSystem jobs(2);
        MemoryNodeCache cache;
        DagExecutor exec(jobs, cache);
        const auto r = exec.run(d, CancelToken{});
        REQUIRE_FALSE(r.ok);
        REQUIRE_FALSE(r.error.empty());
    }
}

TEST_CASE("A 3-node DAG runs in dependency order and produces outputs", "[workflow][dag]") {
    RunLog log;
    Dag d;
    d.add(makeNode("a", {}, log, "A", 1, "x"));
    d.add(makeNode("b", {"a"}, log, "B", 1, "x"));
    d.add(makeNode("c", {"b"}, log, "C", 1, "x"));

    JobSystem jobs(4);
    MemoryNodeCache cache;
    DagExecutor exec(jobs, cache);
    const auto r = exec.run(d, CancelToken{});

    REQUIRE(r.ok);
    REQUIRE(r.ran == 3);
    REQUIRE(r.cached == 0);
    REQUIRE(r.statusOf("c") == NodeStatus::Done);
    REQUIRE(r.output.at("c") == "a,b,c");  // c saw b's output, which saw a's
    REQUIRE(log.count("a") == 1);
}

TEST_CASE("Re-running an unchanged DAG is a full cache hit (resume)", "[workflow][cache]") {
    RunLog log;
    auto build = [&log] {
        Dag d;
        d.add(makeNode("a", {}, log, "A", 1, "x"));
        d.add(makeNode("b", {"a"}, log, "B", 1, "x"));
        d.add(makeNode("c", {"b"}, log, "C", 1, "x"));
        return d;
    };
    JobSystem jobs(4);
    MemoryNodeCache cache;
    DagExecutor exec(jobs, cache);

    const auto r1 = exec.run(build(), CancelToken{});
    REQUIRE(r1.ran == 3);

    const auto r2 = exec.run(build(), CancelToken{});  // same cache
    REQUIRE(r2.ok);
    REQUIRE(r2.cached == 3);
    REQUIRE(r2.ran == 0);
    REQUIRE(r2.output.at("c") == "a,b,c");
    // No node function ran a second time.
    REQUIRE(log.count("a") == 1);
    REQUIRE(log.count("b") == 1);
    REQUIRE(log.count("c") == 1);
}

TEST_CASE("Changing one node re-runs it and only its downstream", "[workflow][cache]") {
    RunLog log;
    JobSystem jobs(4);
    MemoryNodeCache cache;
    DagExecutor exec(jobs, cache);

    Dag d1;
    d1.add(makeNode("a", {}, log, "A", 1, "x"));
    d1.add(makeNode("b", {"a"}, log, "B", 1, "pb"));
    d1.add(makeNode("c", {"b"}, log, "C", 1, "x"));
    REQUIRE(exec.run(d1, CancelToken{}).ran == 3);

    // Same graph but b's params change -> b's key changes -> b and (transitively) c
    // miss; a is unchanged and is served from cache.
    Dag d2;
    d2.add(makeNode("a", {}, log, "A", 1, "x"));
    d2.add(makeNode("b", {"a"}, log, "B", 1, "pb-CHANGED"));
    d2.add(makeNode("c", {"b"}, log, "C", 1, "x"));
    const auto r = exec.run(d2, CancelToken{});

    REQUIRE(r.statusOf("a") == NodeStatus::Cached);
    REQUIRE(r.statusOf("b") == NodeStatus::Done);
    REQUIRE(r.statusOf("c") == NodeStatus::Done);
    REQUIRE(log.count("a") == 1);  // never re-ran
    REQUIRE(log.count("b") == 2);
    REQUIRE(log.count("c") == 2);
}

TEST_CASE("A pre-cancelled run marks all nodes cancelled", "[workflow][cancel]") {
    RunLog log;
    Dag d;
    d.add(makeNode("a", {}, log, "A", 1, "x"));
    d.add(makeNode("b", {"a"}, log, "B", 1, "x"));

    JobSystem jobs(2);
    MemoryNodeCache cache;
    DagExecutor exec(jobs, cache);
    CancelToken tok;
    tok.cancel();

    const auto r = exec.run(d, tok);
    REQUIRE_FALSE(r.ok);
    REQUIRE(r.cancelled == 2);
    REQUIRE(r.ran == 0);
    REQUIRE(log.count("a") == 0);  // nothing executed
}

TEST_CASE("A node that cancels mid-run skips its downstream", "[workflow][cancel]") {
    RunLog log;
    CancelToken tok;
    Dag d;
    // 'a' cancels the run from inside; 'b'/'c' downstream must not execute.
    Node a;
    a.id = "a";
    a.module = "A";
    a.run = [&tok, &log](const NodeInputs&, const CancelToken&) {
        log.bump("a");
        tok.cancel();
        return NodeResult::success("a");
    };
    d.add(a);
    d.add(makeNode("b", {"a"}, log, "B", 1, "x"));
    d.add(makeNode("c", {"b"}, log, "C", 1, "x"));

    JobSystem jobs(2);
    MemoryNodeCache cache;
    DagExecutor exec(jobs, cache);
    const auto r = exec.run(d, tok);

    REQUIRE_FALSE(r.ok);
    REQUIRE(r.statusOf("a") == NodeStatus::Done);
    REQUIRE(r.statusOf("b") == NodeStatus::Cancelled);  // un-started after cancel
    REQUIRE(log.count("c") == 0);
}

TEST_CASE("A failed node skips its downstream but the run reports it", "[workflow][dag]") {
    RunLog log;
    Dag d;
    Node a;
    a.id = "a";
    a.module = "A";
    a.run = [](const NodeInputs&, const CancelToken&) { return NodeResult::failure("boom"); };
    d.add(a);
    d.add(makeNode("b", {"a"}, log, "B", 1, "x"));

    JobSystem jobs(2);
    MemoryNodeCache cache;
    DagExecutor exec(jobs, cache);
    const auto r = exec.run(d, CancelToken{});

    REQUIRE_FALSE(r.ok);
    REQUIRE(r.failed == 1);
    REQUIRE(r.statusOf("a") == NodeStatus::Failed);
    REQUIRE(r.statusOf("b") == NodeStatus::Skipped);
    REQUIRE(log.count("b") == 0);
}

TEST_CASE("Wide DAG runs nodes concurrently and the sink sees all of them",
          "[workflow][concurrency]") {
    RunLog log;
    constexpr int M = 16;
    std::atomic<int> live{0}, maxLive{0}, ran{0};

    Dag d;
    d.add(makeNode("root", {}, log, "ROOT", 1, "x"));
    for (int i = 0; i < M; ++i) {
        Node n;
        n.id = "m" + std::to_string(i);
        n.module = "M";
        n.params = std::to_string(i);
        n.deps = {"root"};
        n.run = [&live, &maxLive, &ran](const NodeInputs&, const CancelToken&) {
            const int now = live.fetch_add(1) + 1;
            int prev = maxLive.load();
            while (now > prev && !maxLive.compare_exchange_weak(prev, now)) {}
            std::this_thread::sleep_for(std::chrono::milliseconds(15));
            ran.fetch_add(1);
            live.fetch_sub(1);
            return NodeResult::success("m");
        };
        d.add(n);
    }
    std::vector<std::string> sinkDeps;
    for (int i = 0; i < M; ++i) sinkDeps.push_back("m" + std::to_string(i));
    d.add(makeNode("sink", sinkDeps, log, "SINK", 1, "x"));

    JobSystem jobs(4);
    MemoryNodeCache cache;
    DagExecutor exec(jobs, cache);
    const auto r = exec.run(d, CancelToken{});

    REQUIRE(r.ok);
    REQUIRE(ran.load() == M);                  // every middle node ran exactly once
    REQUIRE(r.statusOf("sink") == NodeStatus::Done);
    if (jobs.threadCount() >= 2) REQUIRE(maxLive.load() >= 2);  // genuine parallelism
}

TEST_CASE("DiskNodeCache round-trips an output", "[workflow][cache]") {
    DiskNodeCache cache("dag-unit-test");
    const std::string key = "stimlab_unit_key_123";
    REQUIRE_FALSE(cache.get(key).has_value());
    cache.put(key, "hello-payload");
    auto got = cache.get(key);
    REQUIRE(got.has_value());
    REQUIRE(*got == "hello-payload");
}
