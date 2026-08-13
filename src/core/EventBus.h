// core/EventBus.h - thread-safe, typed publish/subscribe.
// Handlers are keyed by event type; publish copies the handler list under lock,
// then invokes outside the lock so subscribers may publish/unsubscribe re-entrantly.
#pragma once

#include <cstddef>
#include <functional>
#include <mutex>
#include <typeindex>
#include <unordered_map>
#include <vector>

namespace biocad {

class EventBus {
public:
    using Token = std::size_t;

    template <typename E>
    Token subscribe(std::function<void(const E&)> fn) {
        std::lock_guard<std::mutex> lock(mutex_);
        const Token token = ++counter_;
        auto& entries = handlers_[std::type_index(typeid(E))];
        entries.push_back(Entry{
            token,
            [callback = std::move(fn)](const void* payload) {
                callback(*static_cast<const E*>(payload));
            }});
        return token;
    }

    template <typename E>
    void publish(const E& event) {
        std::vector<RawHandler> snapshot;
        {
            std::lock_guard<std::mutex> lock(mutex_);
            auto it = handlers_.find(std::type_index(typeid(E)));
            if (it == handlers_.end()) return;
            snapshot.reserve(it->second.size());
            for (const auto& entry : it->second) snapshot.push_back(entry.callback);
        }
        for (const auto& handler : snapshot) handler(&event);
    }

    void unsubscribe(Token token) {
        std::lock_guard<std::mutex> lock(mutex_);
        for (auto& [type, entries] : handlers_) {
            for (auto it = entries.begin(); it != entries.end(); ++it) {
                if (it->token == token) {
                    entries.erase(it);
                    return;
                }
            }
        }
    }

private:
    using RawHandler = std::function<void(const void*)>;
    struct Entry {
        Token token;
        RawHandler callback;
    };

    std::mutex mutex_;
    std::unordered_map<std::type_index, std::vector<Entry>> handlers_;
    Token counter_{0};
};

}  // namespace biocad
