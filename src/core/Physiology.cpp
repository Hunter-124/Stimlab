#include "core/Physiology.h"

#include <fstream>

#include <nlohmann/json.hpp>

#include "core/Assets.h"

namespace biocad::core {
namespace {

Physiology load() {
    Physiology p;
    const auto file = assetDir("packs") / "physiology.json";
    p.path = file.string();
    std::ifstream in(file);
    if (!in) {
        p.error = "assets/packs/physiology.json could not be opened (looked at " + p.path + ")";
        return p;
    }
    nlohmann::json j;
    try {
        in >> j;
    } catch (const std::exception& e) {
        p.error = std::string("assets/packs/physiology.json is not valid JSON: ") + e.what();
        return p;
    }
    if (j.value("schemaVersion", 0) != 1) {
        p.error = "assets/packs/physiology.json has an unsupported schemaVersion";
        return p;
    }
    const auto flows = j.value("flows", nlohmann::json::object());
    auto flow = [&](const char* key) {
        return flows.contains(key) ? flows[key].value("value", 0.0) : 0.0;
    };
    p.hepaticBloodFlowLPerH = flow("hepaticBloodFlow");
    p.enterocyteBloodFlowLPerH = flow("enterocyteBloodFlow");
    p.glomerularFiltrationLPerH = flow("glomerularFiltrationRate");
    const auto deg = j.value("enzymeDegradation", nlohmann::json::object());
    if (deg.contains("values"))
        for (const auto& [k, v] : deg["values"].items())
            if (v.is_number()) p.enzymeDegradation[k] = v.get<double>();
    p.source = j.value("source", std::string{});
    p.loaded = p.hepaticBloodFlowLPerH > 0.0;
    if (!p.loaded) p.error = "assets/packs/physiology.json carries no positive hepatic blood flow";
    return p;
}

}  // namespace

const Physiology& physiology() {
    static const Physiology cached = load();
    return cached;
}

double enzymeDegradationRate(const std::string& enzymeKey) {
    const auto& p = physiology();
    const auto it = p.enzymeDegradation.find(enzymeKey);
    return it == p.enzymeDegradation.end() ? -1.0 : it->second;
}

}  // namespace biocad::core
