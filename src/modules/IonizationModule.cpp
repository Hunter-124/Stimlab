#include "modules/IonizationModule.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>

#include <nlohmann/json.hpp>

#include "chem/Crippen.h"
#include "chem/Descriptors.h"
#include "chem/Formula.h"
#include "chem/Perceive.h"
#include "chem/Solubility.h"
#include "chem/Speciation.h"
#include "core/Assets.h"

namespace biocad {
namespace {

// The pH window and grid every curve in a report shares, so the microspecies
// plot, the net-charge trace and the solubility profile can be read against one
// another without the reader rescaling anything in their head.
constexpr double kPhMin  = 0.0;
constexpr double kPhMax  = 14.0;
constexpr double kPhStep = 0.05;

// The solubility profile is drawn over the physiologically reachable window only:
// a pH-1 stomach to a pH-10 upper limit. Outside it the monoprotic expression is
// arithmetically fine and physically meaningless.
constexpr double kSolubilityPhMin = 1.0;
constexpr double kSolubilityPhMax = 10.0;

// Total molarity assumed when a compound is treated as its own buffer. It is
// stated in BufferReport::assumptions rather than hidden, because buffer capacity
// is linear in this number and a reader who does not know it cannot use the plot.
constexpr double kSelfBufferMolar = 0.1;

// Envelope pruning threshold. 1e-4 keeps every peak a bench mass spectrometer
// would resolve above baseline and discards the combinatorial tail.
constexpr double kEnvelopePrune = 1.0e-4;

std::string requiredString(const nlohmann::json& j, const char* key) {
    if (!j.contains(key) || !j.at(key).is_string()) return {};
    return j.at(key).get<std::string>();
}

std::filesystem::path findPackFile() {
    std::error_code ec;
    if (const char* env = std::getenv("BIOCAD_DESCRIPTOR_DIR")) {
        const auto p = std::filesystem::path(env) / "ionization.json";
        if (std::filesystem::is_regular_file(p, ec)) return p;
    }
    const auto dir = core::assetDir("packs/descriptors");
    if (!dir.empty()) {
        const auto p = dir / "ionization.json";
        if (std::filesystem::is_regular_file(p, ec)) return p;
    }
    return {};
}

// A formula string for a molecule that always exists: the compound's own recorded
// formula when it has one, else the formula perceived from its SMILES. Mass and
// isotope envelope are the two things this phase can ALWAYS compute, so they must
// never depend on the pack. chem::molecularFormula takes the PERCEIVED graph, not
// the catalog DTO, so the SMILES goes through chem::parsePerceived - the raw
// unperceived SMILES parser is never used here, because a number follows.
std::string formulaTextFor(const Molecule& m) {
    if (!m.formula.empty()) return m.formula;
    const auto perceived = chem::parsePerceived(m.smiles);
    return perceived ? chem::molecularFormula(*perceived) : std::string{};
}

// Groups whose pKa sits inside the water window, as buffer components. A group at
// pKa 13.4 contributes nothing measurable to buffer value at any pH a solution is
// actually held at, and including it would put a phantom shoulder on the plot.
std::vector<BufferComponent> selfBufferComponents(const std::vector<IonizableGroup>& groups) {
    std::vector<BufferComponent> out;
    for (const auto& g : groups) {
        if (g.pKa.provenance == Provenance::NotComputed) continue;
        if (g.pKa.value < 1.0 || g.pKa.value > 13.0) continue;
        out.push_back(BufferComponent{g.label, g.pKa.value, kSelfBufferMolar});
    }
    return out;
}

// chem::phSolubility models the monoprotic acid and monoprotic base cases. For a
// polyprotic compound one group has to be chosen, and the choice is recorded in
// the report's warnings rather than made silently: the group nearest the middle of
// the plotted window is the one whose kink is actually visible there.
const IonizableGroup* dominantGroup(const std::vector<IonizableGroup>& groups) {
    const double mid = 0.5 * (kSolubilityPhMin + kSolubilityPhMax);
    const IonizableGroup* best = nullptr;
    double bestDist = 0.0;
    for (const auto& g : groups) {
        if (g.pKa.provenance == Provenance::NotComputed) continue;
        const double d = std::fabs(g.pKa.value - mid);
        if (!best || d < bestDist) {
            best = &g;
            bestDist = d;
        }
    }
    return best;
}

}  // namespace

const IonizationEntry* IonizationPack::find(const std::string& moleculeId) const {
    const auto it = entries.find(moleculeId);
    return it == entries.end() ? nullptr : &it->second;
}

std::filesystem::path defaultIonizationPackPath() { return findPackFile(); }

const char* ionizationInputNote() {
    return "pKa and melting point are INPUTS on this surface, read from the cited ionization "
           "pack or typed by you. BioCAD contains no pKa predictor and no melting-point "
           "predictor and never guesses either one: a compound that is not in the pack gets a "
           "'not computed' naming the input it lacked, not a curve drawn from a plausible "
           "default. Formula, exact mass and the isotope envelope are different in kind - they "
           "are arithmetic on measured isotope masses, so they are always available.";
}

IonizationPack parseIonizationPack(const std::string& text, std::string sourcePath) {
    IonizationPack pack;
    pack.sourcePath = std::move(sourcePath);

    nlohmann::json j;
    try {
        j = nlohmann::json::parse(text);
    } catch (const std::exception& e) {
        pack.errors.push_back(pack.sourcePath + ": not valid JSON - " + e.what());
        return pack;
    }
    if (!j.is_object()) {
        pack.errors.push_back(pack.sourcePath + ": top level is not a JSON object.");
        return pack;
    }

    pack.schemaVersion = j.value("schemaVersion", 0);
    if (pack.schemaVersion != kIonizationSchemaVersion) {
        pack.errors.push_back(pack.sourcePath + ": schemaVersion " +
                              std::to_string(pack.schemaVersion) +
                              " is not understood by this build (expected " +
                              std::to_string(kIonizationSchemaVersion) + ").");
        return pack;
    }
    pack.id          = requiredString(j, "id");
    pack.title       = requiredString(j, "title");
    pack.description = requiredString(j, "description");
    pack.note        = requiredString(j, "note");
    if (pack.note.empty()) {
        pack.errors.push_back(pack.sourcePath +
                              ": no top-level \"note\" - the statement that these values are "
                              "inputs rather than predictions is not optional.");
        return pack;
    }

    if (j.contains("sources")) {
        if (!j.at("sources").is_object()) {
            pack.errors.push_back(pack.sourcePath + ": \"sources\" is not an object.");
            return pack;
        }
        for (const auto& [key, val] : j.at("sources").items()) {
            if (!val.is_string() || val.get<std::string>().empty()) {
                pack.errors.push_back(pack.sourcePath + ": source key \"" + key +
                                      "\" has no prose description.");
                continue;
            }
            pack.sources[key] = val.get<std::string>();
        }
    }

    if (!j.contains("compounds") || !j.at("compounds").is_array()) {
        pack.errors.push_back(pack.sourcePath + ": missing the \"compounds\" array.");
        return pack;
    }

    // A source key resolves to its prose description; anything else is taken as
    // inline prose. Either way the Quantity ends up carrying a describable source,
    // which is the whole requirement - an unattributed pKa is not an input.
    const auto resolveSource = [&pack](const std::string& key) {
        const auto it = pack.sources.find(key);
        return it == pack.sources.end() ? key : it->second;
    };

    std::size_t index = 0;
    for (const auto& e : j.at("compounds")) {
        const std::string where = pack.sourcePath + " compound #" + std::to_string(index++);
        if (!e.is_object()) {
            pack.errors.push_back(where + ": not an object.");
            continue;
        }
        IonizationEntry entry;
        entry.id = requiredString(e, "id");
        if (entry.id.empty()) {
            pack.errors.push_back(where + ": no \"id\".");
            continue;
        }
        if (pack.entries.count(entry.id)) {
            pack.errors.push_back(where + " (" + entry.id +
                                  "): duplicate id - which entry wins would be arbitrary.");
            continue;
        }
        entry.groupsNote = requiredString(e, "groupsNote");

        bool bad = false;
        if (e.contains("groups")) {
            if (!e.at("groups").is_array()) {
                pack.errors.push_back(where + " (" + entry.id + "): \"groups\" is not an array.");
                continue;
            }
            entry.hasGroups = true;
            for (const auto& g : e.at("groups")) {
                if (!g.is_object() || !g.contains("pKa") || !g.at("pKa").is_number()) {
                    pack.errors.push_back(where + " (" + entry.id +
                                          "): a group has no numeric \"pKa\".");
                    bad = true;
                    break;
                }
                const std::string label = requiredString(g, "label");
                const std::string src   = requiredString(g, "source");
                const double      pKa   = g.at("pKa").get<double>();
                if (label.empty()) {
                    pack.errors.push_back(where + " (" + entry.id +
                                          "): a group has no \"label\" - an unnamed site cannot "
                                          "be read off a microspecies plot.");
                    bad = true;
                    break;
                }
                if (src.empty()) {
                    pack.errors.push_back(where + " (" + entry.id + ", " + label +
                                          "): no \"source\" - an uncited pKa is not an input.");
                    bad = true;
                    break;
                }
                if (!(pKa > -5.0 && pKa < 20.0)) {
                    pack.errors.push_back(where + " (" + entry.id + ", " + label + "): pKa " +
                                          std::to_string(pKa) +
                                          " is outside the range a dissociation constant is "
                                          "reported in.");
                    bad = true;
                    break;
                }
                if (!g.contains("acidic") || !g.at("acidic").is_boolean()) {
                    pack.errors.push_back(where + " (" + entry.id + ", " + label +
                                          "): \"acidic\" is required - whether the site loses or "
                                          "gains a proton sets the sign of the charge.");
                    bad = true;
                    break;
                }
                IonizableGroup ig;
                ig.label  = label;
                ig.acidic = g.at("acidic").get<bool>();
                ig.pKa    = makeQuantity(pKa, "", 0.0, Provenance::Measured, resolveSource(src));
                entry.groups.push_back(std::move(ig));
            }
        }
        if (bad) continue;

        if (e.contains("meltingPointC")) {
            if (!e.at("meltingPointC").is_number()) {
                pack.errors.push_back(where + " (" + entry.id +
                                      "): \"meltingPointC\" is not a number.");
                continue;
            }
            entry.meltingPointSource = requiredString(e, "meltingPointSource");
            if (entry.meltingPointSource.empty()) {
                pack.errors.push_back(where + " (" + entry.id +
                                      "): a melting point with no \"meltingPointSource\" - the "
                                      "General Solubility Equation would carry an uncited input "
                                      "into every solubility number it produces.");
                continue;
            }
            entry.meltingPointC   = e.at("meltingPointC").get<double>();
            entry.hasMeltingPoint = true;
            entry.meltingPointSource = resolveSource(entry.meltingPointSource);
        }

        if (e.contains("logPMeasured")) {
            if (!e.at("logPMeasured").is_number()) {
                pack.errors.push_back(where + " (" + entry.id +
                                      "): \"logPMeasured\" is not a number.");
                continue;
            }
            entry.logPSource = requiredString(e, "logPSource");
            if (entry.logPSource.empty()) {
                pack.errors.push_back(where + " (" + entry.id +
                                      "): a measured logP with no \"logPSource\".");
                continue;
            }
            entry.logPMeasured    = e.at("logPMeasured").get<double>();
            entry.hasLogPMeasured = true;
            entry.logPSource      = resolveSource(entry.logPSource);
        }

        if (!entry.hasGroups && !entry.hasMeltingPoint && !entry.hasLogPMeasured) {
            pack.errors.push_back(where + " (" + entry.id +
                                  "): the entry supplies no input at all, so its presence would "
                                  "only disguise a missing one.");
            continue;
        }
        pack.entries.emplace(entry.id, std::move(entry));
    }
    return pack;
}

IonizationPack loadIonizationPack(const std::filesystem::path& file) {
    IonizationPack pack;
    pack.sourcePath = file.string();
    if (file.empty()) {
        pack.errors.push_back("assets/packs/descriptors/ionization.json not found (set "
                              "BIOCAD_DESCRIPTOR_DIR to its directory).");
        return pack;
    }
    std::error_code ec;
    if (!std::filesystem::is_regular_file(file, ec)) {
        pack.errors.push_back(pack.sourcePath + ": file not found.");
        return pack;
    }
    std::ifstream in(file, std::ios::binary);
    if (!in) {
        pack.errors.push_back(pack.sourcePath + ": could not be opened for reading.");
        return pack;
    }
    std::ostringstream buf;
    buf << in.rdbuf();
    return parseIonizationPack(buf.str(), file.string());
}

RealIonization::RealIonization() : pack_(loadIonizationPack(defaultIonizationPackPath())) {}

RealIonization::RealIonization(IonizationPack pack) : pack_(std::move(pack)) {}

std::optional<FormulaMass> RealIonization::formula(const std::string& text) const {
    const auto parsed = chem::parseFormula(text);
    if (!parsed) return std::nullopt;
    return chem::toFormulaMass(*parsed);
}

IsotopeEnvelope RealIonization::envelope(const std::string& formula, double minIntensity) const {
    const auto parsed = chem::parseFormula(formula);
    if (!parsed) {
        // No peaks and a source that says why, rather than an empty plot the reader
        // would read as "this molecule has one isotopologue".
        IsotopeEnvelope env;
        env.formula = formula;
        env.source  = std::string("not computed - \"") + formula +
                     "\" could not be parsed as a molecular formula";
        return env;
    }
    return chem::isotopeEnvelope(*parsed, minIntensity);
}

BalancedEquation RealIonization::balance(const std::vector<std::string>& reactants,
                                        const std::vector<std::string>& products,
                                        const std::vector<double>&      reactantGrams) const {
    return chem::balanceEquation(reactants, products, reactantGrams);
}

SpeciationResult RealIonization::solve(const SpeciationProblem& p) const {
    return chem::solveSpeciation(p);
}

SpeciationResult RealIonization::solvePh(const SpeciationProblem& p) const {
    return chem::solveSpeciationPh(p);
}

SpeciationCurve RealIonization::titrate(const Molecule& m,
                                       const std::vector<IonizableGroup>& groups,
                                       const Quantity&                    logP) const {
    SpeciationCurve c = chem::titrateGroups(groups, logP, kPhMin, kPhMax, kPhStep);
    c.moleculeId = m.id.empty() ? m.name : m.id;
    return c;
}

BufferReport RealIonization::buffer(const std::vector<BufferComponent>& components) const {
    chem::BufferSpec spec;
    spec.components = components;
    return chem::bufferCapacity(spec);
}

SolubilityReport RealIonization::solubility(const Molecule&                    m,
                                           const std::vector<IonizableGroup>& groups,
                                           const Quantity&                    logP,
                                           double meltingPointC) const {
    chem::SolubilityInput in;
    in.moleculeId = m.id.empty() ? m.name : m.id;
    in.pHMin = kSolubilityPhMin;
    in.pHMax = kSolubilityPhMax;
    in.pHStep = kPhStep;

    if (logP.provenance != Provenance::NotComputed) {
        in.logP    = logP.value;
        in.hasLogP = true;
    }
    // meltingPointC <= 0 is the contract's "not supplied". A melting point below
    // 0 C would be a liquid at room temperature, for which the GSE's crystal-
    // lattice term has no meaning anyway.
    if (meltingPointC > 0.0) {
        in.meltingPointC   = meltingPointC;
        in.hasMeltingPoint = true;
    }
    // The catalog's recorded weight, else the average mass of the same formula the
    // mass panel shows - so the two surfaces cannot disagree about one compound.
    double mw = m.molWeight;
    if (mw <= 0.0) {
        if (const auto parsed = chem::parseFormula(formulaTextFor(m))) mw = chem::averageMass(*parsed);
    }
    if (mw > 0.0) {
        in.molWeight    = mw;
        in.hasMolWeight = true;
    }

    std::vector<std::string> polyproticWarning;
    if (const IonizableGroup* g = dominantGroup(groups)) {
        in.kind    = g->acidic ? chem::IonizationKind::MonoproticAcid
                               : chem::IonizationKind::MonoproticBase;
        in.pKa     = g->pKa.value;
        in.hasPKa  = true;
        std::size_t usable = 0;
        for (const auto& other : groups)
            if (other.pKa.provenance != Provenance::NotComputed) ++usable;
        if (usable > 1) {
            polyproticWarning.push_back(
                "This compound has " + std::to_string(usable) +
                " cited ionizable groups, but a pH-solubility profile is a monoprotic model: the "
                "curve below uses \"" + g->label + "\" (pKa " +
                std::to_string(g->pKa.value) +
                "), the group nearest the middle of the plotted window, and ignores the others. "
                "Read the microspecies plot for the full ladder.");
        }
    } else {
        in.kind = chem::IonizationKind::Neutral;
    }

    SolubilityReport r = chem::phSolubility(in);
    r.warnings.insert(r.warnings.end(), polyproticWarning.begin(), polyproticWarning.end());
    return r;
}

Quantity RealIonization::logPFor(const Molecule& m) const {
    // A measured logP wins outright: Wildman-Crippen's own published RMS error is
    // about 0.67 log units, which propagates straight into logD and from there
    // into every partitioning statement, so a shake-flask number is not merely
    // preferable but categorically different in tier.
    const std::string id = m.id.empty() ? m.name : m.id;
    if (const IonizationEntry* e = pack_.find(id); e && e->hasLogPMeasured) {
        return makeQuantity(e->logPMeasured, "", 0.0, Provenance::Measured, e->logPSource);
    }
    const auto perceived = chem::parsePerceived(m.smiles);
    if (!perceived) {
        return notComputed("a parsable SMILES (needed for the Wildman-Crippen logP)");
    }
    const auto cr = chem::crippen(*perceived);
    if (!cr.ok) return notComputed("the Wildman-Crippen descriptor pack (" + cr.note + ")");
    return makeQuantity(cr.logP, "", 0.67, Provenance::Predicted, chem::crippenCitation());
}

IonizationReport RealIonization::analyze(const Molecule& m) const {
    IonizationReport r;
    r.moleculeId = m.id.empty() ? m.name : m.id;

    // ---- Always computable: composition arithmetic on measured isotope masses.
    const std::string formulaText = formulaTextFor(m);
    if (const auto parsed = chem::parseFormula(formulaText)) {
        r.mass     = chem::toFormulaMass(*parsed);
        r.envelope = chem::isotopeEnvelope(*parsed, kEnvelopePrune);
    } else {
        r.mass.formula      = formulaText;
        r.mass.monoisotopic = notComputed("a parsable molecular formula");
        r.mass.average      = notComputed("a parsable molecular formula");
        r.mass.mz           = notComputed("a parsable molecular formula");
        r.mass.warnings.push_back("\"" + formulaText +
                                  "\" could not be parsed as a molecular formula.");
        r.envelope = envelope(formulaText, kEnvelopePrune);
    }

    const Quantity logP = logPFor(m);

    // ---- Pack-dependent from here down.
    const IonizationEntry* e = pack_.find(r.moleculeId);
    if (!e || !e->hasGroups) {
        const std::string missing =
            e ? "pKa (the cited ionization pack has an entry for " + r.moleculeId +
                    " but no group list)"
              : "pKa (not in the cited ionization pack)";
        r.speciation.moleculeId       = r.moleculeId;
        r.speciation.logP             = logP;
        r.speciation.isoelectricPoint = notComputed(missing);
        r.speciation.logDAtPh74       = notComputed(missing);
        r.speciation.assumptions.emplace_back(ionizationInputNote());

        r.solubility.moleculeId       = r.moleculeId;
        r.solubility.intrinsic        = notComputed(
            e && e->hasMeltingPoint
                ? "a measured intrinsic solubility, or a logP to pair with the melting point"
                : "a melting point (the General Solubility Equation needs one)");
        r.solubility.pHmax            = notComputed(missing);
        r.solubility.solubilityAtPh74 = notComputed(missing);
        r.solubility.doseNumber       = notComputed("a dose (BCS Do is dose over solubility)");
        r.solubility.dissolutionNumber = notComputed("a particle radius, density and diffusivity");
        r.solubility.absorptionNumber  = notComputed("an effective permeability Peff");

        r.buffer.betaAtPh74   = notComputed(missing);
        r.buffer.maxCapacity  = notComputed(missing);
        r.buffer.maxCapacityPh = notComputed(missing);
        r.buffer.assumptions.emplace_back(
            "Buffer value is the sum over conjugate pairs, so with no cited pKa there is no pair "
            "to sum: this is an absent input, not a buffer capacity of zero.");

        r.dissolution.timeTo85Pct = notComputed(missing);
        r.dissolution.assumptions.emplace_back(
            "A dissolution time course needs a solubility to dissolve toward, which needs a pKa "
            "and a melting point or a measured S0.");
        return r;
    }

    r.speciation = titrate(m, e->groups, logP);
    if (!e->groupsNote.empty()) r.speciation.assumptions.push_back(e->groupsNote);
    r.speciation.assumptions.emplace_back(ionizationInputNote());

    r.solubility = solubility(m, e->groups, logP,
                              e->hasMeltingPoint ? e->meltingPointC : 0.0);
    if (e->hasMeltingPoint) {
        r.solubility.assumptions.push_back("Melting point " +
                                           std::to_string(e->meltingPointC) +
                                           " C, an input from: " + e->meltingPointSource);
    }

    r.buffer = buffer(selfBufferComponents(e->groups));
    r.buffer.assumptions.push_back(
        "The compound is treated as its own buffer at " + std::to_string(kSelfBufferMolar) +
        " M total per ionizable group; buffer value is linear in that concentration, so scale "
        "accordingly. Groups with pKa outside 1-13 are excluded because their contribution is "
        "below the water terms everywhere.");

    // Dissolution needs FORMULATION inputs - a dose, a particle radius, a density,
    // a diffusivity - none of which are properties of the compound and none of
    // which this pack carries. The panel collects them; analyze() names them.
    r.dissolution.timeTo85Pct = notComputed(
        "a dose, particle radius, density and diffusivity (formulation inputs, entered in the "
        "panel)");
    r.dissolution.assumptions.emplace_back(
        "Dissolution is a property of a formulation, not of a molecule: BioCAD will not assume a "
        "particle size distribution on your behalf.");
    return r;
}

}  // namespace biocad
