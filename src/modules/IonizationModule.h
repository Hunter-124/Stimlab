// modules/IonizationModule.h - the Phase 11 exact-chemistry adapter, plus the
// cited input pack that feeds it.
//
// WHAT THIS FILE IS FOR. chem/Formula.*, chem/Speciation.* and chem/Solubility.*
// know physics but know nothing about compounds. This adapter is the single place
// where a Molecule becomes a set of equilibrium inputs, and therefore the single
// place where the question "do we actually have a pKa for this thing?" is asked
// and answered. There is exactly one answer path: the cited pack, or NotComputed.
//
// WHY THE PACK IS DATA AND NOT A PREDICTOR. A pKa determines the sign of a
// molecule's charge at every pH, and through logD it determines its partitioning.
// It is also, unlike a mass, not derivable - the published structure-based
// estimators disagree with measurement by a pKa unit or more on exactly the
// polyprotic zwitterions where the answer matters most. So BioCAD has no pKa
// predictor and no entry point that could acquire one by accident: pKa, melting
// point, Ksp and precipitation rate constants arrive as inputs or the dependent
// curve is absent. An absent compound yields notComputed() naming the input it
// lacked; it never yields a curve computed from a plausible default, because a
// plausible default is indistinguishable from a measurement once it is plotted.
//
// A malformed pack fails loudly. Every parse failure lands in
// IonizationPack::errors and is rendered by the panel, because a pack that
// silently failed to load is indistinguishable from a compound that genuinely has
// no cited pKa - and those two situations call for opposite responses.
#pragma once

#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include "contracts/IModules.h"
#include "data/Ionization.h"

namespace biocad {

// One compound's measured inputs. `groups` being present but EMPTY is a cited
// statement ("the source reports no dissociation in 0-14"), which is why
// `hasGroups` exists separately from `groups.empty()`.
struct IonizationEntry {
    std::string                 id;
    std::vector<IonizableGroup> groups;
    bool                        hasGroups = false;
    std::string                 groupsNote;

    double      meltingPointC = 0.0;
    bool        hasMeltingPoint = false;
    std::string meltingPointSource;

    double      logPMeasured = 0.0;
    bool        hasLogPMeasured = false;
    std::string logPSource;
};

// One parsed pack document. `errors` is surfaced, never swallowed.
struct IonizationPack {
    int         schemaVersion = 0;
    std::string id;
    std::string title;
    std::string description;
    std::string note;        // the INPUTS-not-predictions statement, rendered verbatim
    std::string sourcePath;
    std::map<std::string, std::string>     sources;   // key -> prose description
    std::map<std::string, IonizationEntry> entries;   // compound id -> inputs
    std::vector<std::string>               errors;

    // nullptr when the compound has no entry, which is what makes every
    // pack-dependent Quantity NotComputed rather than defaulted.
    [[nodiscard]] const IonizationEntry* find(const std::string& moleculeId) const;
};

// The only schema version this build understands. An unknown version is an error
// in `errors`, never a silent skip.
inline constexpr int kIonizationSchemaVersion = 1;

// assets/packs/descriptors/ionization.json, resolved the way every other
// descriptor pack is: BIOCAD_DESCRIPTOR_DIR first, then core::assetDir.
std::filesystem::path defaultIonizationPackPath();

// Parse from disk. A missing or malformed file yields an empty pack with a
// populated `errors`; it never throws.
IonizationPack loadIonizationPack(const std::filesystem::path& file);

// Parse from a JSON string, for tests and harnesses.
IonizationPack parseIonizationPack(const std::string& text,
                                   std::string        sourcePath = "<embedded>");

// The fixed statement the panel and the agent tools both render, so the framing
// cannot drift between surfaces.
const char* ionizationInputNote();

class RealIonization final : public IIonizationModule {
public:
    RealIonization();
    explicit RealIonization(IonizationPack pack);

    std::optional<FormulaMass> formula(const std::string& text) const override;
    IsotopeEnvelope            envelope(const std::string& formula,
                                       double minIntensity) const override;
    BalancedEquation           balance(const std::vector<std::string>& reactants,
                                       const std::vector<std::string>& products,
                                       const std::vector<double>& reactantGrams) const override;
    SpeciationResult           solve(const SpeciationProblem& p) const override;
    SpeciationResult           solvePh(const SpeciationProblem& p) const override;
    SpeciationCurve            titrate(const Molecule& m,
                                       const std::vector<IonizableGroup>& groups,
                                       const Quantity& logP) const override;
    BufferReport               buffer(const std::vector<BufferComponent>& components) const override;
    SolubilityReport           solubility(const Molecule& m,
                                          const std::vector<IonizableGroup>& groups,
                                          const Quantity& logP,
                                          double meltingPointC) const override;
    IonizationReport           analyze(const Molecule& m) const override;

    [[nodiscard]] const IonizationPack& pack() const { return pack_; }

    // The logP one call to analyze() would use for this molecule: the pack's
    // measured value when it has one, else Wildman-Crippen at Provenance::
    // Predicted with crippenCitation() as its source. Exposed because the panel
    // shows which of the two it got, and a second derivation would be a second
    // answer.
    [[nodiscard]] Quantity logPFor(const Molecule& m) const;

private:
    IonizationPack pack_;
};

}  // namespace biocad
