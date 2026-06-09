#include "workflow/JobSystem.h"

#include <algorithm>

namespace stimlab::workflow {

JobSystem::JobSystem(unsigned threads) {
    if (threads == 0) {
        const unsigned hw = std::thread::hardware_concurrency();
        threads = hw > 1 ? hw - 1 : 1;
    }
    workers_.reserve(threads);
    for (unsigned i = 0; i < threads; ++i) workers_.emplace_back([this] { worker(); });
}

JobSystem::~JobSystem() {
    {
        std::lock_guard<std::mutex> lk(mu_);
        stop_ = true;
    }
    cv_.notify_all();
    for (auto& t : workers_)
        if (t.joinable()) t.join();
}

std::future<void> JobSystem::submit(std::function<void()> fn) {
    std::packaged_task<void()> task(std::move(fn));
    std::future<void> fut = task.get_future();
    {
        std::lock_guard<std::mutex> lk(mu_);
        tasks_.push(std::move(task));
    }
    cv_.notify_one();
    return fut;
}

void JobSystem::worker() {
    for (;;) {
        std::packaged_task<void()> task;
        {
            std::unique_lock<std::mutex> lk(mu_);
            cv_.wait(lk, [this] { return stop_ || !tasks_.empty(); });
            // Drain remaining tasks even after stop so handed-out futures complete.
            if (tasks_.empty()) return;
            task = std::move(tasks_.front());
            tasks_.pop();
        }
        task();
    }
}

}  // namespace stimlab::workflow
