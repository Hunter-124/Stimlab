// core/Config.h - simple JSON-backed key/value settings persisted to disk.
#pragma once

#include <filesystem>
#include <mutex>
#include <string>

#include <nlohmann/json.hpp>

namespace stimlab {

class Config {
public:
    explicit Config(std::filesystem::path file);

    void load();              // reads file if present; tolerant of missing/corrupt file
    bool save() const;        // writes pretty JSON; returns false on IO error

    template <typename T>
    T get(const std::string& key, T fallback) const {
        std::lock_guard<std::mutex> lock(mutex_);
        if (auto it = data_.find(key); it != data_.end()) {
            try {
                return it->template get<T>();
            } catch (const std::exception&) {
                return fallback;
            }
        }
        return fallback;
    }

    template <typename T>
    void set(const std::string& key, T value) {
        std::lock_guard<std::mutex> lock(mutex_);
        data_[key] = std::move(value);
    }

    [[nodiscard]] bool has(const std::string& key) const;

private:
    std::filesystem::path file_;
    nlohmann::json data_ = nlohmann::json::object();
    mutable std::mutex mutex_;
};

}  // namespace stimlab
