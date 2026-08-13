#include "modules/docking/Presets.h"

#include <algorithm>
#include <cctype>
#include <mutex>

namespace biocad::docking {
namespace {

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

// The cached view of the loaded packs. Built once (or after reloadTargetPacks())
// and never mutated afterwards, so the pointers findPreset() hands out stay valid.
struct Cache {
    packs::LoadReport report;
    std::vector<ReceptorTarget> dockable;   // targets with a real pdb + box
    std::vector<ReceptorTarget> headline;   // subset flagged headline in a pack
};

void rebuild(Cache& c) {
    c.dockable.clear();
    c.headline.clear();
    for (const auto& t : c.report.targets()) {
        if (!t.hasBox) continue;  // coverage gap, not a dockable preset
        c.dockable.push_back(t.target);
        if (t.headline) c.headline.push_back(t.target);
    }
}

Cache& built() {
    static Cache c;
    static std::once_flag once;
    std::call_once(once, [] {
        c.report = packs::loadBuiltin();
        rebuild(c);
    });
    return c;
}

}  // namespace

const std::vector<ReceptorTarget>& targetPresets() { return built().dockable; }

const packs::LoadReport& targetPackReport() { return built().report; }

void reloadTargetPacks() {
    Cache& c = built();
    c.report = packs::loadBuiltin();
    rebuild(c);
}

const ReceptorTarget* findPreset(const std::string& nameOrId) {
    const std::string key = lower(nameOrId);
    for (const auto& t : targetPresets()) {
        if (lower(t.id) == key || lower(t.name) == key) return &t;
    }
    // Fall back to a prefix match on the display name so "DAT" resolves even when
    // the caller passes "DAT (dopamine transporter)" and vice versa.
    for (const auto& t : targetPresets()) {
        if (lower(t.name).rfind(key, 0) == 0) return &t;
    }
    return nullptr;
}

std::vector<std::string> presetNames() {
    std::vector<std::string> names;
    names.reserve(targetPresets().size());
    for (const auto& t : targetPresets()) names.push_back(t.name);
    return names;
}

std::vector<ReceptorTarget> headlinePresets() { return built().headline; }

}  // namespace biocad::docking
