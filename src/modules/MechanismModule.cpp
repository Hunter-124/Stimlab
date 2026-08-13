#include "modules/MechanismModule.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <future>
#include <sstream>

#include <nlohmann/json.hpp>

#include "agent/WebTools.h"
#include "modules/docking/ReceptorPrep.h"
#include "packs/Pack.h"
#include "workflow/JobSystem.h"

namespace biocad {
namespace {

namespace fs = std::filesystem;

std::string lower(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::vector<std::string> stringArray(const nlohmann::json& j, const char* key) {
    std::vector<std::string> out;
    if (!j.contains(key) || !j.at(key).is_array()) return out;
    for (const auto& v : j.at(key))
        if (v.is_string()) out.push_back(v.get<std::string>());
    return out;
}

std::string readFile(const fs::path& p) {
    std::ifstream in(p, std::ios::binary);
    if (!in) return {};
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

// The single sentence that keeps an empty result honest. An empty list is a
// statement about the query and the source, never about the compound.
std::string coverageNoteFor(const std::string& compoundId, std::size_t entries, bool attempted,
                            bool online) {
    std::string s = "Retrieved records only: BioCAD does not infer a mechanism from structure, "
                    "fingerprint similarity or a docking pose. ";
    if (!online) {
        s += "No retrieval was possible in this build, so the empty list for '" + compoundId +
             "' says nothing at all about the compound.";
        return s;
    }
    if (!attempted) {
        s += "No retrieval was attempted for '" + compoundId +
             "', so the empty list says nothing at all about the compound.";
        return s;
    }
    if (entries == 0) {
        s += "The source returned no mechanism record for '" + compoundId +
             "'. That is a statement about this query and this database's curation, not about the "
             "compound: an uncurated compound and a mechanism-free compound look identical here.";
        return s;
    }
    s += "The " + std::to_string(entries) + " record(s) below are what the source curates for '" +
         compoundId + "'; anything it has not curated is absent, not excluded.";
    return s;
}

constexpr const char* kChemblSource =
    "ChEMBL web services (CC BY-SA 3.0; attribution: ChEMBL, EMBL-EBI)";
constexpr const char* kReactomeSource =
    "Reactome ContentService (CC0; attribution: Reactome, reactome.org)";

}  // namespace

// ---------------------------------------------------------------- pack lookups
bool ActionTypePack::recognises(const std::string& actionType) const {
    return std::any_of(values.begin(), values.end(),
                       [&](const ActionTypeValue& v) { return v.value == actionType; });
}

const StackMember* InteractionPack::find(const std::string& idOrName) const {
    const std::string key = lower(idOrName);
    for (const auto& m : members)
        if (lower(m.id) == key || lower(m.name) == key) return &m;
    return nullptr;
}

const MechanismPanelPack* MechanismPacks::panel(const std::string& panelId) const {
    const std::string key = lower(panelId);
    for (const auto& p : panels)
        if (lower(p.panelId) == key) return &p;
    return nullptr;
}

// ------------------------------------------------------------------ pack loading
fs::path defaultMechanismPackDir() {
    const auto root = packs::builtinPackDir();
    if (root.empty()) return {};
    return root / "mechanism";
}

void parseMechanismDocument(const std::string& text, const std::string& sourcePath,
                            MechanismPacks& into) {
    nlohmann::json j;
    try {
        j = nlohmann::json::parse(text);
    } catch (const std::exception& e) {
        into.errors.push_back(sourcePath + ": not valid JSON (" + e.what() + ")");
        return;
    }
    if (!j.is_object()) {
        into.errors.push_back(sourcePath + ": top level is not an object");
        return;
    }
    const int version = j.value("schemaVersion", 0);
    if (version != kMechanismPackSchemaVersion) {
        into.errors.push_back(sourcePath + ": unsupported schemaVersion " + std::to_string(version) +
                              " (this build understands " +
                              std::to_string(kMechanismPackSchemaVersion) + ")");
        return;
    }
    const std::string kind = j.value("kind", "");

    if (kind == "panel") {
        MechanismPanelPack p;
        p.panelId = j.value("panelId", "");
        p.title = j.value("title", "");
        p.declaredSize = j.value("declaredSize", 0);
        p.licence = j.value("licence", "");
        p.sourceNote = j.value("sourceNote", "");
        p.coverageNote = j.value("coverageNote", "");
        if (p.panelId.empty()) {
            into.errors.push_back(sourcePath + ": panel pack has no panelId");
            return;
        }
        for (const auto& t : j.value("targets", nlohmann::json::array())) {
            MechanismPanelTarget row;
            row.id = t.value("id", "");
            row.name = t.value("name", "");
            row.targetClass = t.value("targetClass", "");
            row.targetRef = t.value("targetRef", "");
            if (row.id.empty()) {
                into.errors.push_back(sourcePath + ": a target row has no id");
                continue;
            }
            p.targets.push_back(std::move(row));
        }
        // A declared size below the enumerated roster would understate the coverage
        // gap, which is the one number this module exists to protect.
        if (p.declaredSize < static_cast<int>(p.targets.size()))
            p.declaredSize = static_cast<int>(p.targets.size());
        into.panels.push_back(std::move(p));
        return;
    }

    if (kind == "action-types") {
        ActionTypePack v;
        v.id = j.value("id", "");
        v.title = j.value("title", "");
        v.declaredSize = j.value("declaredSize", 0);
        v.licence = j.value("licence", "");
        v.sourceNote = j.value("sourceNote", "");
        for (const auto& e : j.value("values", nlohmann::json::array())) {
            ActionTypeValue av;
            av.value = e.value("value", "");
            av.gloss = e.value("gloss", "");
            if (!av.value.empty()) v.values.push_back(std::move(av));
        }
        if (v.values.empty()) {
            into.errors.push_back(sourcePath + ": action-type pack has no values");
            return;
        }
        into.actionTypes = std::move(v);
        return;
    }

    if (kind == "interactions") {
        InteractionPack p;
        p.id = j.value("id", "");
        p.title = j.value("title", "");
        p.licence = j.value("licence", "");
        p.boundaryNote = j.value("boundaryNote", "");
        p.coverageNote = j.value("coverageNote", "");
        if (p.boundaryNote.empty()) {
            into.errors.push_back(sourcePath +
                                  ": interaction pack has no boundaryNote; every flag must carry "
                                  "the not-a-recommendation sentence, so the pack is rejected");
            return;
        }
        for (const auto& e : j.value("members", nlohmann::json::array())) {
            StackMember m;
            m.id = e.value("id", "");
            m.name = e.value("name", "");
            m.kind = e.value("kind", "");
            m.substrateOf = stringArray(e, "substrateOf");
            m.inhibitorOf = stringArray(e, "inhibitorOf");
            m.inducerOf = stringArray(e, "inducerOf");
            m.mechanisms = stringArray(e, "mechanisms");
            m.citation = e.value("citation", "");
            m.note = e.value("note", "");
            if (m.id.empty()) continue;
            if (m.citation.empty()) {
                into.errors.push_back(sourcePath + ": member '" + m.id +
                                      "' has no citation and was dropped");
                continue;
            }
            p.members.push_back(std::move(m));
        }
        for (const auto& e : j.value("classRules", nlohmann::json::array())) {
            StackClassRule r;
            r.left = e.value("left", "");
            r.right = e.value("right", "");
            r.mechanism = e.value("mechanism", "");
            r.evidence = e.value("evidence", "");
            r.citation = e.value("citation", "");
            if (r.left.empty() || r.right.empty() || r.citation.empty()) {
                into.errors.push_back(sourcePath +
                                      ": a classRule is missing left/right/citation and was dropped");
                continue;
            }
            p.classRules.push_back(std::move(r));
        }
        into.interactions = std::move(p);
        return;
    }

    if (kind == "pharmacogenomics") {
        PharmacogenomicsPack p;
        p.id = j.value("id", "");
        p.title = j.value("title", "");
        p.licence = j.value("licence", "");
        p.boundaryStatement = j.value("boundaryStatement", "");
        for (const auto& e : j.value("phenotypes", nlohmann::json::array())) {
            PgxPhenotype ph;
            ph.code = e.value("code", "");
            ph.term = e.value("term", "");
            // "extensive metabolizer" is deprecated nomenclature. A pack that ships it
            // is rejected at the boundary rather than rendered.
            if (lower(ph.term).find("extensive metabolizer") != std::string::npos) {
                into.errors.push_back(sourcePath + ": phenotype '" + ph.code +
                                      "' uses the deprecated term 'extensive metabolizer'");
                continue;
            }
            if (!ph.code.empty()) p.phenotypes.push_back(std::move(ph));
        }
        for (const auto& e : j.value("activityScoreBands", nlohmann::json::array())) {
            PgxActivityBand b;
            b.gene = e.value("gene", "");
            b.phenotype = e.value("phenotype", "");
            b.band = e.value("band", "");
            if (!b.gene.empty() && !b.band.empty()) p.bands.push_back(std::move(b));
        }
        for (const auto& e : j.value("notes", nlohmann::json::array())) {
            PgxNote n;
            n.compound = e.value("compound", "");
            n.gene = e.value("gene", "");
            n.phenotype = e.value("phenotype", "");
            n.implication = e.value("implication", "");
            if (!n.compound.empty() && !n.implication.empty()) p.notes.push_back(std::move(n));
        }
        into.pharmacogenomics = std::move(p);
        return;
    }

    into.errors.push_back(sourcePath + ": unknown pack kind '" + kind + "'");
}

MechanismPacks loadMechanismPacks(const fs::path& dir) {
    MechanismPacks out;
    out.sourceDir = dir.string();
    std::error_code ec;
    if (dir.empty() || !fs::is_directory(dir, ec)) {
        out.errors.push_back("mechanism pack directory not found: " +
                             (dir.empty() ? std::string("<unresolved>") : dir.string()));
        return out;
    }
    std::vector<fs::path> files;
    for (const auto& e : fs::directory_iterator(dir, ec))
        if (e.is_regular_file(ec) && e.path().extension() == ".json") files.push_back(e.path());
    std::sort(files.begin(), files.end());  // deterministic load order
    for (const auto& f : files) {
        const std::string text = readFile(f);
        if (text.empty()) {
            out.errors.push_back(f.string() + ": empty or unreadable");
            continue;
        }
        parseMechanismDocument(text, f.string(), out);
    }
    return out;
}

bool mechanismNetworkAvailable() {
#ifdef BIOCAD_HAVE_SCIENCE
    return true;
#else
    return false;
#endif
}

const char* mechanismOfflineReason() {
    return "not retrieved: this build has no network transport compiled in (the science feature is "
           "off), so nothing was queried. This is an absence of retrieval, not an absence of "
           "mechanism.";
}

std::string keggPathwayUrl(const std::string& mapId) {
    // A page URL, deliberately not an API URL: there is no code path in BioCAD that
    // fetches kegg.jp, and adding one would need a KEGG licence.
    return "https://www.kegg.jp/pathway/" + mapId;
}

const char* keggDeepLinkNote() {
    return "KEGG is not queried and none of its content is bundled: its REST API is licensed for "
           "academic use only. This is a deep link to a public KEGG pathway page for you to open "
           "yourself - linking is not redistribution.";
}

// ------------------------------------------------------------------ the module
RealMechanism::RealMechanism(const IDockingModule* docking,
                             std::vector<packs::PackCompound> catalog)
    : docking_(docking), catalog_(std::move(catalog)),
      packs_(loadMechanismPacks(defaultMechanismPackDir())) {}

RealMechanism::RealMechanism(const IDockingModule* docking,
                             std::vector<packs::PackCompound> catalog, MechanismPacks packs)
    : docking_(docking), catalog_(std::move(catalog)), packs_(std::move(packs)) {}

MechanismReport RealMechanism::mechanisms(const std::string& compoundId) const {
    MechanismReport r;
    r.compoundId = compoundId;
    r.networkAvailable = mechanismNetworkAvailable();

    // Resolve the ChEMBL cross-reference from the catalog. A raw CHEMBL id is
    // accepted directly so the panel can query something not in the packs.
    std::string chembl;
    if (lower(compoundId).rfind("chembl", 0) == 0) {
        chembl = compoundId;
    } else {
        const std::string key = lower(compoundId);
        for (const auto& c : catalog_)
            if (lower(c.id) == key || lower(c.name) == key) { chembl = c.xrefs.chembl; break; }
    }

    if (!r.networkAvailable) {
        r.warnings.emplace_back(mechanismOfflineReason());
        r.warnings.emplace_back(
            "Mechanism records are retrieved from ChEMBL (CC BY-SA 3.0); there is no bundled "
            "mechanism table to fall back on, because a stale local copy would be presented as a "
            "current retrieval.");
        r.coverageNote = coverageNoteFor(compoundId, 0, false, false);
        return r;
    }
    if (chembl.empty()) {
        r.warnings.push_back("no ChEMBL cross-reference for '" + compoundId +
                             "' in the loaded compound packs, so no query could be formed; add "
                             "xrefs.chembl to the pack entry or pass a CHEMBL id directly");
        r.coverageNote = coverageNoteFor(compoundId, 0, false, true);
        return r;
    }

    r.retrievalAttempted = true;
    const std::string url = "https://www.ebi.ac.uk/chembl/api/data/mechanism.json"
                            "?molecule_chembl_id=" + chembl + "&limit=50";
    const agent::ApiGetResult got = agent::apiGet(url);
    if (!got.ok) {
        r.warnings.push_back("ChEMBL mechanism query failed (" + got.error + "); nothing retrieved");
        r.coverageNote = coverageNoteFor(compoundId, 0, true, true);
        return r;
    }

    nlohmann::json j;
    try {
        j = nlohmann::json::parse(got.body);
    } catch (const std::exception& e) {
        r.warnings.push_back(std::string("ChEMBL response was not valid JSON (") + e.what() + ")");
        r.coverageNote = coverageNoteFor(compoundId, 0, true, true);
        return r;
    }

    for (const auto& m : j.value("mechanisms", nlohmann::json::array())) {
        MechanismEntry e;
        e.actionType = m.value("action_type", "");
        // Free text is carried VERBATIM. It is never parsed into a vocabulary,
        // because that would invent structure the source does not have.
        if (m.contains("mechanism_of_action") && m.at("mechanism_of_action").is_string())
            e.freeTextMechanism = m.at("mechanism_of_action").get<std::string>();
        e.source = kChemblSource;
        e.provenance = Provenance::Measured;
        // BioCAD's own orthogonal axis. ChEMBL cannot express inhibition modality, so
        // it stays Unknown until an assay fit supplies one; it is never guessed from
        // the action type.
        e.modality = InhibitionModality::Unknown;
        for (const auto& ref : m.value("mechanism_refs", nlohmann::json::array())) {
            std::string s = ref.value("ref_type", "");
            const std::string id = ref.value("ref_id", "");
            const std::string ru = ref.contains("ref_url") && ref.at("ref_url").is_string()
                                       ? ref.at("ref_url").get<std::string>() : std::string();
            if (!id.empty()) s += " " + id;
            if (!ru.empty()) s += " <" + ru + ">";
            if (!s.empty()) e.references.push_back(std::move(s));
        }
        if (!e.actionType.empty() && !packs_.actionTypes.recognises(e.actionType))
            r.warnings.push_back("action_type '" + e.actionType +
                                 "' is not in the bundled controlled vocabulary; it is reported "
                                 "verbatim rather than mapped onto a neighbouring value");

        // Target identity, from ChEMBL's own target record.
        const std::string tid = m.contains("target_chembl_id") && m.at("target_chembl_id").is_string()
                                    ? m.at("target_chembl_id").get<std::string>() : std::string();
        e.targetName = tid;
        if (!tid.empty()) {
            const agent::ApiGetResult t =
                agent::apiGet("https://www.ebi.ac.uk/chembl/api/data/target/" + tid + ".json");
            if (t.ok) {
                try {
                    const auto tj = nlohmann::json::parse(t.body);
                    if (tj.contains("pref_name") && tj.at("pref_name").is_string())
                        e.targetName = tj.at("pref_name").get<std::string>();
                    if (tj.contains("organism") && tj.at("organism").is_string())
                        e.organism = tj.at("organism").get<std::string>();
                    for (const auto& comp : tj.value("target_components", nlohmann::json::array())) {
                        if (comp.contains("accession") && comp.at("accession").is_string()) {
                            e.targetAccession = comp.at("accession").get<std::string>();
                            break;
                        }
                    }
                } catch (const std::exception&) {
                    r.warnings.push_back("ChEMBL target record " + tid + " was not valid JSON");
                }
            } else {
                r.warnings.push_back("ChEMBL target record " + tid + " not retrieved (" + t.error +
                                     "); the target id is shown instead of its name");
            }
        }
        r.entries.push_back(std::move(e));
    }

    if (got.fromCache)
        r.warnings.emplace_back("served from the local ChEMBL response cache under cache/api");
    r.warnings.emplace_back(
        "ChEMBL content is CC BY-SA 3.0: attribution is required and anything derived from these "
        "records inherits the share-alike obligation.");
    r.coverageNote = coverageNoteFor(compoundId, r.entries.size(), true, true);
    return r;
}

PanelScreenReport RealMechanism::screenPanel(const Molecule& m, const std::string& panelId) const {
    PanelScreenReport r;
    r.panelId = panelId;

    const MechanismPanelPack* pack = packs_.panel(panelId);
    if (!pack) {
        r.warnings.push_back("no panel pack with id '" + panelId +
                             "' is loaded; nothing was screened");
        r.coverageStatement = "0 of an unknown number of targets screened: the panel roster itself "
                              "could not be loaded, so the coverage gap cannot even be counted.";
        r.hergSafetyMargin = notComputed("panel roster");
        return r;
    }
    for (const std::string& e : packs_.errors) r.warnings.push_back("pack: " + e);

    r.panelSize = pack->declaredSize;

    // Resolve every row to a dockable receptor preset FIRST, so the unscreened count
    // is known before any docking runs.
    std::vector<ReceptorTarget> presets;
    if (docking_) presets = docking_->presets();
    auto findPreset = [&presets](const std::string& ref) -> const ReceptorTarget* {
        if (ref.empty()) return nullptr;
        const std::string key = lower(ref);
        for (const auto& t : presets)
            if (lower(t.id) == key) return &t;
        return nullptr;
    };

    r.results.reserve(pack->targets.size());
    std::vector<const ReceptorTarget*> jobs;   // parallel to r.results, nullptr = skip
    for (const MechanismPanelTarget& t : pack->targets) {
        PanelTargetResult row;
        row.targetId = t.id;
        row.targetName = t.name.empty() ? t.id : t.name;
        row.affinity = notComputed("a docked pose for " + t.id);
        const ReceptorTarget* p = findPreset(t.targetRef);
        if (!docking_) {
            row.skipReason = "no docking module is wired in this build";
        } else if (!p) {
            row.skipReason = "no receptor preset with a binding box in the loaded target packs"
                             " (declared coverage gap, not a docking failure)";
        }
        jobs.push_back(row.skipReason.empty() ? p : nullptr);
        r.results.push_back(std::move(row));
    }

    // Dock the resolvable rows in parallel on the EXISTING job system. Each row is
    // independent; nothing is aggregated across them, by design.
    if (std::any_of(jobs.begin(), jobs.end(), [](const ReceptorTarget* p) { return p != nullptr; })) {
        workflow::JobSystem pool;
        std::vector<std::future<void>> futures;
        futures.reserve(jobs.size());
        for (std::size_t i = 0; i < jobs.size(); ++i) {
            if (!jobs[i]) continue;
            const ReceptorTarget* target = jobs[i];
            PanelTargetResult* row = &r.results[i];
            futures.push_back(pool.submit([this, &m, target, row] {
                const docking::ReceptorPrepResult prep =
                    docking::locatePreparedReceptor(target->id);
                // Every row states what was actually docked into and in which box.
                // Two rows from different preparations are not comparable, which is
                // why the report offers no field in which to compare them.
                row->receptorPreparation =
                    prep.ready ? ("prepared receptor " + prep.path + " (PDB " + target->pdb + ")")
                               : ("no prepared receptor cached for PDB " + target->pdb);
                std::ostringstream box;
                box.setf(std::ios::fixed);
                box.precision(1);
                const bool own = prep.ready && prep.hasBox;
                box << "center (" << (own ? prep.cx : target->box.cx) << ", "
                    << (own ? prep.cy : target->box.cy) << ", " << (own ? prep.cz : target->box.cz)
                    << ") size (" << target->box.sx << ", " << target->box.sy << ", "
                    << target->box.sz << ") A, "
                    << (own ? "center from the prepared receptor's own frame"
                            : "center from the target pack");
                row->boxDefinition = box.str();

                const DockJobResult d = docking_->dockDetailed(m, target->id);
                if (d.provenance == Provenance::Model && !d.poses.empty()) {
                    row->screened = true;
                    row->affinity = makeQuantity(d.bestAffinity(), "kcal/mol", 0.0,
                                                 Provenance::Model, d.engine);
                } else {
                    // The descriptor fallback is not a screen. Recording it as one
                    // would move a target out of the unknown column on the strength
                    // of a rank-ordering heuristic.
                    row->screened = false;
                    row->skipReason = "no docking engine produced a pose (" +
                                      (d.log.empty() ? std::string("engine or prepared receptor "
                                                                   "unavailable")
                                                     : d.log) +
                                      "); the descriptor estimate is not a screen";
                    row->affinity = notComputed("a docked pose from a real engine");
                }
            }));
        }
        for (auto& f : futures) f.get();
    }

    r.screened = static_cast<int>(std::count_if(r.results.begin(), r.results.end(),
                                                [](const PanelTargetResult& x) { return x.screened; }));
    r.unscreened = r.panelSize - r.screened;

    const int enumerated = static_cast<int>(pack->targets.size());
    std::ostringstream cs;
    cs << r.unscreened << " of " << r.panelSize << " targets in '" << pack->panelId
       << "' were NOT screened; " << r.screened << " were. The unknown fraction therefore "
       << (r.unscreened >= r.screened ? "dominates" : "remains material")
       << " this result: a panel says nothing about the " << r.unscreened
       << " targets it did not interrogate.";
    if (enumerated < r.panelSize)
        cs << " " << (r.panelSize - enumerated) << " of those are not even enumerated in the "
           << "roster, so BioCAD cannot name them.";
    cs << " There is no composite safety score and no cross-target ranking here: receptor "
       << "preparations, box volumes and rotatable-bond penalties differ per row, so the scores "
       << "are not on a common scale.";
    if (!pack->coverageNote.empty()) cs << " " << pack->coverageNote;
    r.coverageStatement = cs.str();

    // hERG margin, only from user-supplied MEASURED inputs.
    if (herg_.measuredIc50Molar <= 0.0 && herg_.freeCmaxMolar <= 0.0) {
        r.hergSafetyMargin = notComputed(
            "a measured hERG IC50 and a measured free Cmax (both user-supplied; BioCAD does not "
            "predict a hERG IC50 and derives no QT or TdP risk)");
    } else if (herg_.measuredIc50Molar <= 0.0) {
        r.hergSafetyMargin = notComputed("a measured hERG IC50 (user-supplied only)");
    } else if (herg_.freeCmaxMolar <= 0.0) {
        r.hergSafetyMargin = notComputed("a measured free Cmax");
    } else {
        std::string src = "measured hERG IC50 / measured free Cmax";
        if (!herg_.citation.empty()) src += "; IC50 from " + herg_.citation;
        r.hergSafetyMargin = makeQuantity(herg_.measuredIc50Molar / herg_.freeCmaxMolar, "", 0.0,
                                          Provenance::Measured, src);
        if (r.hergSafetyMargin.value < r.hergMarginFlagBelow)
            r.warnings.push_back("hERG margin is below the flag threshold of " +
                                 std::to_string(static_cast<int>(r.hergMarginFlagBelow)) +
                                 "-fold; that is a flag on the ratio of two measurements, not a "
                                 "proarrhythmia prediction");
    }
    return r;
}

PathwayContext RealMechanism::pathways(const std::string& uniprotAccession) const {
    PathwayContext r;
    r.accession = uniprotAccession;
    r.source = kReactomeSource;
    r.networkAvailable = mechanismNetworkAvailable();
    r.warnings.emplace_back(
        "Pathway membership only. There is no pathway impact score: no database supports "
        "propagating a docking score through a pathway graph, so such a number would be "
        "fabrication with a scientific veneer.");
    if (!r.networkAvailable) {
        r.warnings.emplace_back(mechanismOfflineReason());
        return r;
    }
    if (uniprotAccession.empty()) {
        r.warnings.emplace_back("no UniProt accession given, so no query could be formed");
        return r;
    }

    const agent::ApiGetResult got = agent::apiGet(
        "https://reactome.org/ContentService/data/mapping/UniProt/" + uniprotAccession + "/pathways");
    if (!got.ok) {
        r.warnings.push_back("Reactome mapping query failed (" + got.error + "); nothing retrieved");
        return r;
    }
    nlohmann::json j;
    try {
        j = nlohmann::json::parse(got.body);
    } catch (const std::exception& e) {
        r.warnings.push_back(std::string("Reactome response was not valid JSON (") + e.what() + ")");
        return r;
    }
    if (!j.is_array()) {
        r.warnings.emplace_back("Reactome returned no pathway list for this accession");
        return r;
    }
    for (const auto& p : j) {
        PathwayNode n;
        n.stableId = p.value("stId", "");
        n.name = p.value("displayName", "");
        n.species = p.value("speciesName", "");
        if (n.stableId.empty()) continue;
        n.url = "https://reactome.org/content/detail/" + n.stableId;
        // Ancestors turn the flat mapping into the hierarchy the panel renders.
        const agent::ApiGetResult anc = agent::apiGet(
            "https://reactome.org/ContentService/data/event/" + n.stableId + "/ancestors");
        if (anc.ok) {
            try {
                const auto aj = nlohmann::json::parse(anc.body);
                // The endpoint returns an array of ancestor CHAINS, outermost first.
                for (const auto& chain : aj) {
                    if (!chain.is_array()) continue;
                    for (const auto& a : chain) {
                        const std::string id = a.value("stId", "");
                        if (!id.empty() && id != n.stableId &&
                            std::find(n.ancestorIds.begin(), n.ancestorIds.end(), id) ==
                                n.ancestorIds.end())
                            n.ancestorIds.push_back(id);
                    }
                }
            } catch (const std::exception&) {
                r.warnings.push_back("ancestors for " + n.stableId + " were not valid JSON");
            }
        } else {
            r.warnings.push_back("ancestors for " + n.stableId + " not retrieved (" + anc.error +
                                 "); it is rendered as a root");
        }
        r.pathways.push_back(std::move(n));
    }
    if (r.pathways.empty())
        r.warnings.emplace_back("Reactome curates no pathway for this accession, which is a "
                                "statement about Reactome's coverage and not about the protein");
    return r;
}

StackReport RealMechanism::checkStack(const std::vector<std::string>& memberIds) const {
    StackReport r;
    const InteractionPack& pack = packs_.interactions;
    r.coverageNote = pack.coverageNote.empty()
                         ? std::string("no interaction pack is loaded, so no flag could be raised")
                         : pack.coverageNote;
    for (const std::string& e : packs_.errors) r.warnings.push_back("pack: " + e);

    std::vector<const StackMember*> resolved;
    for (const std::string& id : memberIds) {
        if (id.empty()) continue;
        const StackMember* m = pack.find(id);
        if (!m) {
            r.members.push_back(id);
            r.unknownMembers.push_back(id);
            continue;
        }
        r.members.push_back(m->id);
        resolved.push_back(m);
    }
    if (!r.unknownMembers.empty())
        r.warnings.push_back(std::to_string(r.unknownMembers.size()) +
                             " member(s) are not in the interaction pack and were NOT screened; "
                             "they are listed in unknownMembers rather than silently ignored");

    auto has = [](const std::vector<std::string>& v, const std::string& x) {
        return std::find(v.begin(), v.end(), x) != v.end();
    };

    for (std::size_t a = 0; a < resolved.size(); ++a) {
        for (std::size_t b = 0; b < resolved.size(); ++b) {
            if (a == b) continue;
            const StackMember& p = *resolved[a];   // candidate perpetrator
            const StackMember& v = *resolved[b];   // candidate victim

            for (const std::string& cyp : p.inhibitorOf) {
                if (!has(v.substrateOf, cyp)) continue;
                InteractionFlag f;
                f.leftId = p.id;
                f.rightId = v.id;
                f.mechanism = cyp + " inhibition by " + p.name + "; " + v.name + " is a " + cyp +
                              " substrate";
                f.direction = p.name + " is the perpetrator; " + v.name + " is the affected member";
                f.evidence = "enzyme role assignment for both members";
                f.citation = p.citation;
                f.boundaryNote = pack.boundaryNote;
                r.flags.push_back(std::move(f));
            }
            for (const std::string& cyp : p.inducerOf) {
                if (!has(v.substrateOf, cyp)) continue;
                InteractionFlag f;
                f.leftId = p.id;
                f.rightId = v.id;
                f.mechanism = cyp + " induction by " + p.name + "; " + v.name + " is a " + cyp +
                              " substrate";
                f.direction = p.name + " is the perpetrator; " + v.name + " is the affected member";
                f.evidence = "enzyme role assignment for both members";
                f.citation = p.citation;
                f.boundaryNote = pack.boundaryNote;
                r.flags.push_back(std::move(f));
            }
            for (const StackClassRule& rule : pack.classRules) {
                if (!has(p.mechanisms, rule.left) || !has(v.mechanisms, rule.right)) continue;
                if (a > b && rule.left == rule.right) continue;  // symmetric rule, once per pair
                InteractionFlag f;
                f.leftId = p.id;
                f.rightId = v.id;
                f.mechanism = rule.mechanism + " (" + p.name + " + " + v.name + ")";
                f.direction = rule.left == rule.right
                                  ? std::string("both members contribute the same mechanism")
                                  : (p.name + " contributes " + rule.left + "; " + v.name +
                                     " contributes " + rule.right);
                f.evidence = rule.evidence;
                f.citation = rule.citation;
                f.boundaryNote = pack.boundaryNote;
                r.flags.push_back(std::move(f));
            }
        }
    }

    if (r.flags.empty() && resolved.size() > 1)
        r.warnings.emplace_back("no mechanism in this pack links the known members. That is the "
                                "pack's coverage speaking, not a safety conclusion.");
    return r;
}

PharmacogenomicReport RealMechanism::pharmacogenomics(const std::string& compoundId) const {
    PharmacogenomicReport r;
    r.compoundId = compoundId;
    // The CPIC content is a bundled CC0 pack; no live query is made, so this is
    // false in EVERY build and the warning says why rather than leaving the reader
    // to guess whether a lookup failed.
    r.networkAvailable = false;
    r.boundaryStatement = packs_.pharmacogenomics.boundaryStatement;
    if (r.boundaryStatement.empty())
        r.boundaryStatement = "Conditional notes only: no genotype is interpreted, no phenotype is "
                              "assigned, and no dose or dose adjustment is emitted.";
    if (!mechanismNetworkAvailable())
        r.warnings.emplace_back(mechanismOfflineReason());
    r.warnings.emplace_back(
        "notes come from the bundled CPIC pack (CC0), not from a live query, so they are as current "
        "as this build's assets and no more");
    for (const std::string& e : packs_.errors) r.warnings.push_back("pack: " + e);

    const std::string key = lower(compoundId);
    for (const PgxNote& n : packs_.pharmacogenomics.notes) {
        if (lower(n.compound) != key) continue;
        PharmacogenomicNote out;
        out.gene = n.gene;
        out.phenotype = n.phenotype;
        out.implication = n.implication;
        out.source = packs_.pharmacogenomics.licence;
        out.citation = "CPIC guideline for " + n.gene +
                       "; CYP2D6 activity-score bands from Caudle and co-workers, Clinical and "
                       "Translational Science (2020)";
        for (const PgxActivityBand& b : packs_.pharmacogenomics.bands)
            if (b.gene == n.gene && b.phenotype == n.phenotype) { out.activityScoreBand = b.band; break; }
        r.notes.push_back(std::move(out));
    }
    if (r.notes.empty())
        r.warnings.push_back("no CPIC note in the bundled pack mentions '" + compoundId +
                             "'; absence here means the pack has no note, not that the compound has "
                             "no pharmacogenomic relevance");
    return r;
}

}  // namespace biocad
