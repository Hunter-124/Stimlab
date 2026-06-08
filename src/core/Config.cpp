#include "core/Config.h"

#include <fstream>

namespace stimlab {

Config::Config(std::filesystem::path file) : file_(std::move(file)) {}

void Config::load() {
    std::lock_guard<std::mutex> lock(mutex_);
    std::ifstream in(file_);
    if (!in) {
        data_ = nlohmann::json::object();
        return;
    }
    try {
        in >> data_;
        if (!data_.is_object()) data_ = nlohmann::json::object();
    } catch (const std::exception&) {
        data_ = nlohmann::json::object();  // corrupt file -> start fresh
    }
}

bool Config::save() const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::error_code ec;
    std::filesystem::create_directories(file_.parent_path(), ec);
    std::ofstream out(file_, std::ios::trunc);
    if (!out) return false;
    out << data_.dump(2);
    return static_cast<bool>(out);
}

bool Config::has(const std::string& key) const {
    std::lock_guard<std::mutex> lock(mutex_);
    return data_.contains(key);
}

}  // namespace stimlab
