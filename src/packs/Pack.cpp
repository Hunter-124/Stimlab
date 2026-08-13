#include "packs/Pack.h"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <system_error>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "core/AppPaths.h"
#include "core/Error.h"

#if defined(_WIN32)
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#endif

namespace biocad::packs {
namespace {

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

// A required string field. Missing or empty is a structural error, because a
// compound without an id cannot be referenced and a target without a name cannot
// be shown.
std::string requireString(const nlohmann::json& j, const char* key, const std::string& where) {
    const auto it = j.find(key);
    if (it == j.end() || !it->is_string() || it->get<std::string>().empty()) {
        throw Error::parse(where + ": missing required string field \"" + key + "\"");
    }
    return it->get<std::string>();
}

std::string optString(const nlohmann::json& j, const char* key) {
    const auto it = j.find(key);
    return (it != j.end() && it->is_string()) ? it->get<std::string>() : std::string{};
}

Xrefs parseXrefs(const nlohmann::json& j) {
    Xrefs x;
    const auto it = j.find("xrefs");
    if (it == j.end() || !it->is_object()) return x;
    x.chembl = optString(*it, "chembl");
    x.inchikey = optString(*it, "inchikey");
    x.uniprot = optString(*it, "uniprot");
    if (const auto cid = it->find("pubchemCid");
        cid != it->end() && cid->is_number_integer()) {
        x.pubchemCid = cid->get<long long>();
    }
    return x;
}

PackCompound parseCompound(const nlohmann::json& j, const std::string& where) {
    PackCompound c;
    c.id = requireString(j, "id", where);
    c.name = requireString(j, "name", where);
    c.smiles = requireString(j, "smiles", where + " (" + c.id + ")");
    c.drugClass = optString(j, "drugClass");
    c.legalStatus = optString(j, "legalUs");
    if (c.legalStatus.empty()) c.legalStatus = optString(j, "legalStatus");
    c.notes = optString(j, "notes");
    c.xrefs = parseXrefs(j);

    const auto props = j.find("properties");
    if (props != j.end() && props->is_object()) {
        c.formula = optString(*props, "formula");
        c.molWeight = props->value("molWeight", 0.0);
        c.logP = props->value("logP", 0.0);
        c.tpsa = props->value("tpsa", 0.0);
        c.hbd = props->value("hbd", 0);
        c.hba = props->value("hba", 0);
        c.rotatableBonds = props->value("rotatableBonds", 0);
        c.hasProperties = c.molWeight > 0.0;
    }
    return c;
}

PackTarget parseTarget(const nlohmann::json& j, const std::string& where) {
    PackTarget t;
    t.target.id = requireString(j, "id", where);
    t.target.name = requireString(j, "name", where);
    if (const auto pdb = j.find("pdb"); pdb != j.end() && pdb->is_string()) {
        t.target.pdb = pdb->get<std::string>();
    }
    t.uniprot = optString(j, "uniprot");
    t.chembl = optString(j, "chembl");
    t.headline = j.value("headline", false);

    if (const auto box = j.find("box"); box != j.end() && box->is_object()) {
        DockBox b;
        b.cx = box->value("cx", 0.0);
        b.cy = box->value("cy", 0.0);
        b.cz = box->value("cz", 0.0);
        b.sx = box->value("sx", 22.0);
        b.sy = box->value("sy", 22.0);
        b.sz = box->value("sz", 22.0);
        if (b.sx <= 0.0 || b.sy <= 0.0 || b.sz <= 0.0) {
            throw Error::parse(where + " (" + t.target.id + "): box edges must be positive");
        }
        t.target.box = b;
        t.hasBox = true;
    }
    if (t.hasBox && t.target.pdb.empty()) {
        throw Error::parse(where + " (" + t.target.id +
                           "): a box without a PDB id is an unverifiable binding site");
    }

    if (const auto panels = j.find("panels"); panels != j.end() && panels->is_array()) {
        for (const auto& p : *panels) {
            if (p.is_string()) t.panels.push_back(p.get<std::string>());
        }
    }
    return t;
}

}  // namespace

Molecule PackCompound::molecule() const {
    Molecule m;
    m.id = id;
    m.name = name;
    m.smiles = smiles;
    m.drugClass = drugClass;
    m.legalStatus = legalStatus;
    m.notes = notes;
    if (hasProperties) {
        m.formula = formula;
        m.molWeight = molWeight;
        m.logP = logP;
        m.tpsa = tpsa;
        m.hbd = hbd;
        m.hba = hba;
        m.rotatableBonds = rotatableBonds;
    }
    return m;
}

Pack parse(const nlohmann::json& j, std::string sourcePath) {
    if (!j.is_object()) throw Error::parse(sourcePath + ": pack must be a JSON object");

    const auto sv = j.find("schemaVersion");
    if (sv == j.end() || !sv->is_number_integer()) {
        throw Error::parse(sourcePath + ": missing integer \"schemaVersion\"");
    }
    Pack p;
    p.schemaVersion = sv->get<int>();
    if (p.schemaVersion != kSchemaVersion) {
        throw Error::unsupported(sourcePath + ": schemaVersion " +
                                 std::to_string(p.schemaVersion) + " is not supported (this build reads " +
                                 std::to_string(kSchemaVersion) + ")");
    }

    p.id = requireString(j, "id", sourcePath);
    p.title = requireString(j, "title", sourcePath);
    p.description = optString(j, "description");
    p.sourcePath = std::move(sourcePath);

    std::unordered_set<std::string> seenCompounds;
    if (const auto cs = j.find("compounds"); cs != j.end()) {
        if (!cs->is_array()) throw Error::parse(p.sourcePath + ": \"compounds\" must be an array");
        p.compounds.reserve(cs->size());
        for (const auto& c : *cs) {
            auto compound = parseCompound(c, p.sourcePath);
            if (!seenCompounds.insert(compound.id).second) {
                throw Error::parse(p.sourcePath + ": duplicate compound id \"" + compound.id +
                                   "\" inside one pack");
            }
            p.compounds.push_back(std::move(compound));
        }
    }

    std::unordered_set<std::string> seenTargets;
    if (const auto ts = j.find("targets"); ts != j.end()) {
        if (!ts->is_array()) throw Error::parse(p.sourcePath + ": \"targets\" must be an array");
        p.targets.reserve(ts->size());
        for (const auto& t : *ts) {
            auto target = parseTarget(t, p.sourcePath);
            if (!seenTargets.insert(target.target.id).second) {
                throw Error::parse(p.sourcePath + ": duplicate target id \"" + target.target.id +
                                   "\" inside one pack");
            }
            p.targets.push_back(std::move(target));
        }
    }
    return p;
}

Pack parseString(const std::string& text, std::string sourcePath) {
    nlohmann::json j = nlohmann::json::parse(text, nullptr, /*allow_exceptions=*/false);
    if (j.is_discarded()) throw Error::parse(sourcePath + ": not valid JSON");
    return parse(j, std::move(sourcePath));
}

LoadReport loadFrom(const std::filesystem::path& dir, bool builtin) {
    LoadReport report;
    std::error_code ec;
    if (!std::filesystem::is_directory(dir, ec)) return report;

    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::directory_iterator(dir, ec)) {
        if (entry.is_regular_file(ec) && entry.path().extension() == ".json") {
            files.push_back(entry.path());
        }
    }
    std::sort(files.begin(), files.end());  // deterministic load order

    for (const auto& file : files) {
        std::ifstream in(file, std::ios::binary);
        if (!in) {
            report.errors.push_back(file.string() + ": cannot open");
            continue;
        }
        std::ostringstream buf;
        buf << in.rdbuf();
        try {
            Pack p = parseString(buf.str(), file.string());
            p.builtin = builtin;
            report.packs.push_back(std::move(p));
        } catch (const Error& e) {
            report.errors.push_back(e.message);
        } catch (const std::exception& e) {
            report.errors.push_back(file.string() + ": " + e.what());
        }
    }
    return report;
}

std::filesystem::path builtinPackDir() {
    std::error_code ec;
#if defined(_WIN32)
    wchar_t buf[MAX_PATH]{};
    if (GetModuleFileNameW(nullptr, buf, MAX_PATH) > 0) {
        const auto exeDir = std::filesystem::path(buf).parent_path();
        // Shipped layout: assets/packs beside the exe.
        if (std::filesystem::is_directory(exeDir / "assets" / "packs", ec)) {
            return exeDir / "assets" / "packs";
        }
        // Dev layout: build/<preset>/bin/BioCAD.exe with the source tree three up.
        const auto devDir = exeDir.parent_path().parent_path().parent_path() / "assets" / "packs";
        if (std::filesystem::is_directory(devDir, ec)) return devDir;
    }
#endif
    if (std::filesystem::is_directory("assets/packs", ec)) return "assets/packs";
    return {};
}

LoadReport loadBuiltin() {
    LoadReport report = loadFrom(builtinPackDir(), /*builtin=*/true);
    if (report.packs.empty() && report.errors.empty()) {
        report.errors.push_back(
            "no built-in packs found - the application has no compound or target catalog. "
            "Expected assets/packs/*.json beside the executable.");
    }

    LoadReport user = loadFrom(AppPaths::instance().root() / "packs", /*builtin=*/false);
    for (auto& p : user.packs) report.packs.push_back(std::move(p));
    for (auto& e : user.errors) report.errors.push_back(std::move(e));
    return report;
}

namespace {

// Later packs win by pack id: a user pack named "cns-monoamine" replaces the
// built-in one wholesale rather than merging into it, so the user can always
// tell exactly which document produced a row.
std::vector<const Pack*> effectivePacks(const std::vector<Pack>& packs) {
    std::unordered_map<std::string, std::size_t> byId;
    std::vector<const Pack*> out;
    for (const auto& p : packs) {
        const auto it = byId.find(p.id);
        if (it == byId.end()) {
            byId.emplace(p.id, out.size());
            out.push_back(&p);
        } else {
            out[it->second] = &p;  // later wins
        }
    }
    return out;
}

}  // namespace

std::vector<PackCompound> LoadReport::compounds() {
    std::vector<PackCompound> out;
    std::unordered_map<std::string, std::string> seen;  // compound id -> pack id
    for (const Pack* p : effectivePacks(packs)) {
        for (const auto& c : p->compounds) {
            const auto it = seen.find(c.id);
            if (it != seen.end()) {
                errors.push_back("duplicate compound id \"" + c.id + "\" in packs \"" +
                                 it->second + "\" and \"" + p->id +
                                 "\" - the later definition was dropped");
                continue;
            }
            seen.emplace(c.id, p->id);
            out.push_back(c);
        }
    }
    return out;
}

std::vector<PackTarget> LoadReport::targets() {
    std::vector<PackTarget> out;
    std::unordered_map<std::string, std::string> seen;  // target id -> pack id
    for (const Pack* p : effectivePacks(packs)) {
        for (const auto& t : p->targets) {
            const auto it = seen.find(t.target.id);
            if (it != seen.end()) {
                errors.push_back("duplicate target id \"" + t.target.id + "\" in packs \"" +
                                 it->second + "\" and \"" + p->id +
                                 "\" - the later definition was dropped");
                continue;
            }
            seen.emplace(t.target.id, p->id);
            out.push_back(t);
        }
    }
    return out;
}

std::optional<PackTarget> LoadReport::findTarget(const std::string& nameOrId) const {
    const std::string key = lower(nameOrId);
    const auto effective = effectivePacks(packs);
    for (const Pack* p : effective) {
        for (const auto& t : p->targets) {
            if (lower(t.target.id) == key || lower(t.target.name) == key) return t;
        }
    }
    // Prefix match so "DAT" resolves "DAT (dopamine transporter)" and vice versa.
    for (const Pack* p : effective) {
        for (const auto& t : p->targets) {
            if (lower(t.target.name).rfind(key, 0) == 0) return t;
        }
    }
    return std::nullopt;
}

}  // namespace biocad::packs
