// workflow/JobSystem.h - a small fixed-size worker pool + cooperative cancel token.
//
// The DAG executor (Dag.h) submits node work here so independent nodes run in
// parallel without the UI thread ever blocking. CancelToken is a shared, thread-safe
// flag: a workflow run holds one, the UI flips it, and long-running node functions
// poll cancelled() to stop cooperatively. Deliberately tiny and dependency-free
// (std::thread + a queue) - the project hand-rolls this rather than pulling Taskflow
// because the value-add (content-addressed caching / resume / cancel) lives in the
// DAG layer, not the scheduler.
//
// SAFETY SCOPE: this only schedules analysis work (prep / dock / screen). It carries
// no synthesis/route content and is domain-agnostic.
#pragma once

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

namespace biocad::workflow {

// A shared cooperative cancellation flag. Copies share one atomic, so the UI can
// cancel a run that node functions (which captured a copy) observe via cancelled().
class CancelToken {
public:
    CancelToken() : flag_(std::make_shared<std::atomic<bool>>(false)) {}

    void cancel() const noexcept { flag_->store(true, std::memory_order_relaxed); }
    [[nodiscard]] bool cancelled() const noexcept {
        return flag_->load(std::memory_order_relaxed);
    }

private:
    std::shared_ptr<std::atomic<bool>> flag_;
};

// Fixed-size thread pool. submit() enqueues a task and returns a future the caller
// can wait on. The destructor drains all queued tasks, then joins - so every future
// handed out is satisfied (no broken promises) before teardown.
class JobSystem {
public:
    // threads == 0 -> hardware_concurrency() - 1, clamped to >= 1.
    explicit JobSystem(unsigned threads = 0);
    ~JobSystem();
    JobSystem(const JobSystem&) = delete;
    JobSystem& operator=(const JobSystem&) = delete;

    std::future<void> submit(std::function<void()> fn);
    [[nodiscard]] unsigned threadCount() const noexcept {
        return static_cast<unsigned>(workers_.size());
    }

private:
    void worker();

    std::vector<std::thread>               workers_;
    std::queue<std::packaged_task<void()>> tasks_;
    std::mutex                             mu_;
    std::condition_variable                cv_;
    bool                                   stop_ = false;
};

}  // namespace biocad::workflow
