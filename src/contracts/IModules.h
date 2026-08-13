// contracts/IModules.h - frozen analysis-module interfaces.
// Every analysis capability is a pure-virtual interface; real impls (RDKit-backed)
// implement these. The UI codes ONLY against these contracts, never against an
// implementation - which is what lets a module be replaced without touching a panel.
//
// SAFETY SCOPE: interfaces describe identity, pharmacology (binding affinity),
// stability, absorption/PK, similarity and legal status. There is intentionally
// no synthesis/route/manufacturability interface.
#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "bio/Structure.h"
#include "contracts/IDockingBackend.h"
#include "data/Assay.h"
#include "data/Biologics.h"
#include "data/Domain.h"
#include "data/Ionization.h"
#include "data/Nucleic.h"
#include "data/Population.h"
#include "data/Systems.h"

namespace biocad {

// Curated compound library (defaults + user imports).
class ILibrary {
public:
    virtual ~ILibrary() = default;
    virtual std::vector<Molecule> all() const = 0;
    virtual std::optional<Molecule> byId(const std::string& id) const = 0;
    virtual std::size_t count() const = 0;
};

// Molecular stability (replaces the out-of-scope manufacturability score).
class IStabilityModule {
public:
    virtual ~IStabilityModule() = default;
    virtual StabilityReport analyze(const Molecule& m) const = 0;
};

// ADMET / metabolism (D, M, E, T) - harmful metabolites, DDIs, safety flags.
class IAdmetModule {
public:
    virtual ~IAdmetModule() = default;
    virtual AdmetReport screen(const Molecule& m) const = 0;
};

// Structural alerts for metabolic bioactivation.
//
// SAFETY SCOPE: this interface emits LIABILITY FLAGS. A flag says a substructure
// has been associated with reactive-metabolite formation in the literature; it
// does not say the compound is toxic, and there is deliberately no toxicity
// verdict, no composite risk score and no Verdict::Danger anywhere in
// AlertReport. An empty flag list means "no alert in the pack matched", which is
// not a safety claim about the compound.
class IAlertsModule {
public:
    virtual ~IAlertsModule() = default;
    virtual AlertReport screen(const Molecule& m) const = 0;
};

// Curated, cited biotransformations for a library compound.
//
// SAFETY SCOPE: this interface RETRIEVES facts; it never enumerates hypotheses.
// Rule-based metabolite prediction was measured at 1.1-29% precision and
// 14.7-28.3% sensitivity across every major tool (Boyce et al. 2022, Comput
// Toxicol 21:100208), so hypothesis generation is a separate surface and there is
// no entry point here that would produce one. An empty MetabolismReport::known is
// a statement about BioCAD's curation, not about the compound, which is why
// coverageNote is a required field of the result.
class IMetabolismFactsModule {
public:
    virtual ~IMetabolismFactsModule() = default;
    virtual MetabolismReport known(const Molecule& m) const = 0;
};

// Absorption / pharmacokinetics (the "A" of ADMET) - permeability, F%, BBB.
class IAbsorptionModule {
public:
    virtual ~IAbsorptionModule() = default;
    virtual AbsorptionReport predict(const Molecule& m) const = 0;
};

// Structural + pharmacophore similarity vs the known-substance reference set.
class ISimilarityModule {
public:
    virtual ~ISimilarityModule() = default;
    virtual SimilarityReport search(const Molecule& m) const = 0;
};

// Legal-analog "substantially similar" scorecard.
class ILegalModule {
public:
    virtual ~ILegalModule() = default;
    virtual LegalScorecard score(const Molecule& m) const = 0;
};

// Docking = ligand->protein binding affinity (pharmacology/activity).
class IDockingModule {
public:
    virtual ~IDockingModule() = default;
    virtual std::vector<std::string> targets() const = 0;
    virtual DockingResult dock(const Molecule& m, const std::string& target) const = 0;

    // Phase D: richer result carrying 3D poses (for the molecular viewer) and the
    // engine/fallback provenance. Default = no 3D detail, so a legacy impl keeps
    // compiling; the real backend-backed impl overrides both.
    virtual DockJobResult dockDetailed(const Molecule& m, const std::string& target) const {
        (void)m; (void)target; return {};
    }
    // CNS receptor presets (PDB ref + binding-site box). Default = none.
    virtual std::vector<ReceptorTarget> presets() const { return {}; }
};

// Run history.
class IRunStore {
public:
    virtual ~IRunStore() = default;
    virtual std::vector<RunRecord> recent() const = 0;
    // Phase D: persist a run. Default no-op so a legacy store keeps compiling;
    // the SQLite-backed store overrides it.
    virtual void record(const RunRecord&) {}
};

// Pharmacodynamics and pharmacokinetics.
//
// SAFETY SCOPE: this interface emits EXPOSURE SCENARIOS under stated assumptions.
// It has no dose-recommendation entry point, and adding one is out of scope by
// design. Every returned number carries a Provenance; a missing prerequisite is a
// NotComputed Quantity naming what is missing, never a silently assumed default.
class IPharmacodynamicsModule {
public:
    virtual ~IPharmacodynamicsModule() = default;

    // Four-parameter logistic fit in log-concentration space.
    virtual CurveFit fitFourParameterLogistic(const std::vector<DoseResponsePoint>&) const = 0;

    // Cheng-Prusoff Ki from an IC50. Returns NotComputed naming the missing field
    // when the modality's required inputs are absent.
    virtual Quantity kiFromIc50(const ChengPrusoffInput&) const = 0;

    // Schild regression: pA2, slope with a confidence interval, and KB only when
    // the slope CI includes 1.
    virtual SchildResult schild(const std::vector<SchildPoint>&) const = 0;

    // Integrate an exposure profile for a dosing regimen.
    virtual PkProfile simulate(const PkModelSpec&, const DoseRegimen&) const = 0;

    // Fractional target occupancy from a free-concentration profile and a Kd.
    virtual OccupancyCurve occupancy(const PkProfile&, const Quantity& kd) const = 0;
};

// Pairwise sequence alignment (Gotoh affine gaps, BLOSUM62 from a pack).
//
// SAFETY SCOPE: alignGlobal() returns a NotComputed E-value. Karlin-Altschul
// statistics are defined for local alignments only; offering an E-value on a
// global alignment would attach a significance claim to a number that has none.
class ISequenceModule {
public:
    virtual ~ISequenceModule() = default;
    virtual SequenceAlignment alignGlobal(const std::string& a, const std::string& b) const = 0;
    virtual SequenceAlignment alignLocal(const std::string& a, const std::string& b) const = 0;
};

// Protein structure I/O and comparison.
//
// There is deliberately no entry point taking a chem::Conformer: that type is a
// small-molecule distance-geometry embedding, and comparing it as a protein
// would yield an lDDT/TM-score that looks meaningful and is not. bio::Structure
// has no such constructor, so the misuse is a compile error.
class IStructureModule {
public:
    virtual ~IStructureModule() = default;

    // Reads a local .pdb/.cif file (including one previously downloaded into the
    // cache). std::nullopt when the file is unreadable or the format is unknown;
    // recoverable problems arrive in Structure::warnings instead.
    virtual std::optional<bio::Structure> load(const std::filesystem::path& file) const = 0;

    virtual StructureComparison compare(const bio::Structure& ref,
                                        const bio::Structure& model) const = 0;

    // Solvent-accessible surface area. `source` carries the full parameter
    // string (algorithm, probe radius, point count, radii set, hydrogen policy)
    // because a SASA without them is not reproducible.
    virtual Quantity sasa(const bio::Structure& s) const = 0;
};

// Exact chemistry: formula/mass arithmetic, acid-base speciation, buffers and
// pH-dependent solubility/dissolution.
//
// SAFETY SCOPE: pKa, melting point, Ksp and precipitation rate constants are
// INPUTS. There is no entry point that predicts a pKa, and no curve is returned
// from an assumed default: analyze() fills the dependent field with a
// NotComputed Quantity naming the prerequisite it lacked. balance() is
// stoichiometry on an equation the user wrote and deliberately has no route,
// condition, precursor or scale-up output.
class IIonizationModule {
public:
    virtual ~IIonizationModule() = default;

    // Parses a formula string (nesting, isotope labels, trailing charge) into
    // both masses. std::nullopt only when the string is not a formula.
    virtual std::optional<FormulaMass> formula(const std::string& text) const = 0;

    // Theoretical isotope envelope by sparse convolution, pruned below
    // `minIntensity` relative to the base peak.
    virtual IsotopeEnvelope envelope(const std::string& formula,
                                     double minIntensity) const = 0;

    // Integer-null-space balancing plus limiting reagent, theoretical yield and
    // atom economy when the caller supplied reactant amounts (grams, parallel to
    // `reactants`; empty for balancing only).
    virtual BalancedEquation balance(const std::vector<std::string>& reactants,
                                     const std::vector<std::string>& products,
                                     const std::vector<double>& reactantGrams) const = 0;

    // General polyprotic equilibrium by damped Newton on the component tableau.
    virtual SpeciationResult solve(const SpeciationProblem& p) const = 0;

    // pH of a solution whose totals are given, by closing on charge balance.
    virtual SpeciationResult solvePh(const SpeciationProblem& p) const = 0;

    // Microspecies fractions, net charge and logD over a pH range for one
    // molecule's ionizable groups.
    virtual SpeciationCurve titrate(const Molecule& m,
                                    const std::vector<IonizableGroup>& groups,
                                    const Quantity& logP) const = 0;

    virtual BufferReport buffer(const std::vector<BufferComponent>& components) const = 0;

    // pH-solubility. `meltingPointC` <= 0 means "not supplied", which makes
    // SolubilityReport::intrinsic NotComputed rather than defaulted.
    virtual SolubilityReport solubility(const Molecule& m,
                                        const std::vector<IonizableGroup>& groups,
                                        const Quantity& logP,
                                        double meltingPointC) const = 0;

    // Everything the panel renders for one compound, with pKa and melting point
    // taken from the compound's pack entry when present.
    virtual IonizationReport analyze(const Molecule& m) const = 0;
};

// Real experimental data: plate import, QC, curve/kinetics/thermodynamics fits,
// and prospective assay-design simulation.
//
// SAFETY SCOPE: this interface consumes measurements and returns fits. A raw well
// value is Measured and a fitted parameter is Model, and nothing here converts a
// potency into a dose. The design entry point simulates plates, which is
// experimental design; it is not a procedure for making anything.
class IAssayModule {
public:
    virtual ~IAssayModule() = default;

    // Parses long CSV/TSV or a 96/384/1536 grid export. Unknown columns survive as
    // metadata; recoverable problems arrive in AssayDataset::warnings rather than
    // as a silent drop. std::nullopt only when the text is not tabular at all.
    virtual std::optional<AssayDataset> import(const std::string& text) const = 0;

    // Per-plate QC. Z-prime is NotComputed without both a positive and a negative
    // control, and edge/row/column effects are reported, never corrected.
    virtual QcReport qc(const Plate& p) const = 0;

    // Fits one concentration series (or one sensorgram/melt/isotherm) with the
    // requested model. `robust` selects Tukey-biweight IRLS over plain weighting.
    virtual FitResult fit(const std::vector<Well>& series, AssayModel model,
                          bool robust) const = 0;

    // Fits every candidate model and ranks by AICc. `decisive` is false when the
    // top two differ by less than 2, which is the only honest answer there.
    virtual ModelComparison compare(const std::vector<Well>& series,
                                    const std::vector<AssayModel>& candidates) const = 0;

    // Global fit over the full [S] x [I] matrix. This is the ONE producer of
    // InhibitionModality; it returns Unknown rather than guessing when AICc cannot
    // separate the modalities.
    virtual ModelComparison inhibitionModality(const std::vector<Well>& matrix) const = 0;

    // Forward-simulates plates from a stated truth model and error structure, then
    // sends them through the same import/QC/fit path, and reports what the design
    // would actually recover.
    virtual AssayDesignReport simulate(const AssayDesignSpec& spec) const = 0;
};

// DNA/RNA workbench: sequence and feature I/O, restriction mapping, translation,
// oligo thermodynamics, codon metrics and guide search.
//
// BIOSECURITY BOUNDARY, permanent: there is no synthesis-vendor entry point, no
// order export, no batch pathogen design and no therapeutic/germline framing.
// findGuides() returns the exact reference it searched and how many bases it
// examined, because an off-target count without its scope is not a specificity
// claim. optimizeCodons() performs constraint satisfaction; it does not predict
// expression, and the DTO has no field in which to pretend otherwise.
class INucleicAcidModule {
public:
    virtual ~INucleicAcidModule() = default;

    // FASTA or GenBank. std::nullopt when the text is neither; recoverable
    // problems (unknown IUPAC symbol, unparseable location) arrive in
    // NucRecord::warnings.
    virtual std::optional<NucRecord> parse(const std::string& text) const = 0;

    // FASTA/GenBank round-trip. These are the ONLY export formats.
    virtual std::string toFasta(const NucRecord& r) const = 0;
    virtual std::string toGenBank(const NucRecord& r) const = 0;

    virtual std::string reverseComplement(const std::string& seq) const = 0;

    virtual TranslationResult translate(const NucRecord& r, int geneticCodeId,
                                        int minOrfAminoAcids) const = 0;

    virtual RestrictionDigest digest(const NucRecord& r,
                                     const std::vector<std::string>& enzymes) const = 0;

    // Nearest-neighbour thermodynamics. Salt and oligo concentrations are inputs
    // with no defaults worth hiding: a Tm without them is not reproducible.
    virtual OligoThermo oligo(const std::string& seq, double naMolar, double mgMolar,
                              double oligoMolar, double dntpMolar) const = 0;

    virtual std::vector<SecondaryStructure> selfStructures(const std::string& seq,
                                                           double naMolar) const = 0;

    virtual std::vector<PrimerPair> designPrimers(const NucRecord& r, int begin, int end,
                                                  double targetTmC) const = 0;

    virtual CodonMetrics codonMetrics(const std::string& cds,
                                      const std::string& usageTable) const = 0;

    virtual CodonOptimizationResult optimizeCodons(
        const std::string& cds, const std::string& usageTable,
        const std::vector<std::string>& forbiddenSites) const = 0;

    // `reference` is whatever the user actually supplied - a plasmid, a contig, a
    // genome. The result states which, and how many bases were searched.
    virtual GuideSearchResult findGuides(const NucRecord& target, const NucRecord& reference,
                                        const std::string& pam) const = 0;
};

// Population PK, noncompartmental analysis and mechanistic drug interactions.
//
// SAFETY SCOPE, permanent: this interface emits EXPOSURE SCENARIOS. There is no
// entry point that returns a dose, a dose change or a regimen, and adding one is
// out of scope by design. An AUC ratio is a statement about exposure under stated
// in vitro inputs; converting it into a personal dose adjustment is a clinical
// decision that this software does not make.
//
// Reproducibility is part of the contract: simulate() must produce byte-identical
// output for the same VariabilitySpec including its seed, because a band that
// cannot be reproduced cannot be checked.
class IPopulationPkModule {
public:
    virtual ~IPopulationPkModule() = default;

    virtual PopulationProfile simulate(const PkModelSpec& model, const DoseRegimen& regimen,
                                       const VariabilitySpec& variability) const = 0;

    // Noncompartmental analysis of observed data. lambda_z is selected strictly
    // after Tmax by adjusted R-squared over at least three terminal points, and
    // everything extrapolated is flagged unreliable above 20% extrapolation.
    virtual NcaResult nca(const ConcentrationSeries& observed) const = 0;

    // FDA screening R-values plus the mechanistic static AUCR. Returns
    // NotComputed for the AUCR when fm is absent, and reports hepatic-only when
    // Fg is absent, rather than assuming either.
    virtual InteractionReport interaction(const PerpetratorSpec& perpetrator,
                                          const VictimSpec& victim) const = 0;

    // Dynamic enzyme activity. At constant inhibitor concentration its steady
    // state must equal the static model exactly; the result carries both so the
    // agreement is visible rather than asserted.
    virtual EnzymeTimeCourse enzymeTimeCourse(const PerpetratorSpec& perpetrator,
                                              double horizonH) const = 0;

    // Renal and hepatic impairment as explicit, editable exposure scenarios.
    virtual ImpairmentScenario impairment(const VictimSpec& victim, double renalFunctionRatio,
                                          double hepaticClintRatio) const = 0;
};

// Antibody and protein-biologics analysis.
//
// HONESTY SCOPE, permanent: number() returns NotComputed numbering when an IMGT
// anchor fails, because plausible wrong numbering is worse than none.
// "Closest germline set" is never species identification. alanineScan() is
// GEOMETRIC and unit-free Heuristic - it is not a delta-delta-G, and there is no
// entry point here for humanness, immunogenicity, affinity maturation,
// aggregation, viscosity, titre or shelf life.
class IBiologicsModule {
public:
    virtual ~IBiologicsModule() = default;

    virtual AbDomain number(const std::string& sequence, NumberingScheme scheme) const = 0;

    // The same domain re-rendered in another scheme. IMGT is canonical
    // internally; the other schemes are table-driven views of it.
    virtual AbDomain convertScheme(const AbDomain& domain, NumberingScheme to) const = 0;

    // Pack-defined, cited motif flags. `structure` may be null; when it is,
    // exposure is reported as unknown rather than assumed.
    virtual std::vector<SequenceLiability> liabilities(const AbDomain& domain,
                                                       const bio::Structure* structure) const = 0;

    virtual DevelopabilityReport developability(const std::vector<std::string>& chains,
                                                const bio::Structure* structure) const = 0;

    virtual MassLadder massLadder(const std::vector<std::string>& chains,
                                  int disulfideCount) const = 0;

    virtual PeptideMap digest(const std::string& chain, const std::string& protease,
                              int maxMissedCleavages) const = 0;

    virtual InterfaceReport interfaceOf(const bio::Structure& complex,
                                        const std::string& chainsA,
                                        const std::string& chainsB) const = 0;

    virtual AlanineScanReport alanineScan(const bio::Structure& complex,
                                          const std::string& chainsA,
                                          const std::string& chainsB) const = 0;
};

// Reaction-network simulation, chemical kinetics, metabolic flux, pathway
// enrichment and graph metrics.
//
// SAFETY SCOPE: a simulated growth rate is a property of the model, never a
// measurement of an organism. No docking score may propagate into a flux, a
// pathway or a network number, and there is no entry point that would let it.
// Every returned TimeCourse embeds the solver settings that produced it.
class ISimulationModule {
public:
    virtual ~ISimulationModule() = default;

    // Fills the structural analysis (conservation laws, Wegscheider cycles) of a
    // network the caller assembled.
    virtual NetworkSpec analyze(const NetworkSpec& network) const = 0;

    virtual TimeCourse integrate(const NetworkSpec& network, double horizon, double relTol,
                                 double absTol, const std::string& method) const = 0;

    virtual StochasticEnsemble stochastic(const NetworkSpec& network, double horizon,
                                          int replicates, std::uint64_t seed,
                                          bool tauLeap) const = 0;

    virtual KineticsFit arrhenius(const std::vector<double>& temperaturesK,
                                  const std::vector<double>& rateConstants) const = 0;

    virtual PhRateProfile phRate(const std::vector<double>& pHValues,
                                 const std::vector<double>& rateConstants) const = 0;

    virtual ControlAnalysis controlAnalysis(const NetworkSpec& network,
                                            double horizon) const = 0;

    // SBML Level 3 Version 2 Core subset. std::nullopt when the document uses a
    // construct BioCAD does not implement - which is reported by name, never
    // silently ignored.
    virtual std::optional<NetworkSpec> importSbml(const std::string& xml,
                                                  std::string* error) const = 0;
    virtual std::string exportSbml(const NetworkSpec& network) const = 0;
};

// Constraint-based metabolic flux. Gated at build time by BIOCAD_ENABLE_FBA.
//
// balance() MUST pass before fba() may run: a flux distribution over a reaction
// that does not conserve mass and charge is arithmetic about nothing.
class IFluxModule {
public:
    virtual ~IFluxModule() = default;

    virtual FluxSolution balance(const NetworkSpec& network) const = 0;
    virtual FluxSolution fba(const NetworkSpec& network, const std::string& objectiveReactionId,
                             const std::vector<FluxBound>& bounds) const = 0;
    virtual std::vector<FluxRange> fva(const NetworkSpec& network,
                                       const std::string& objectiveReactionId,
                                       const std::vector<FluxBound>& bounds,
                                       double objectiveFraction) const = 0;
    virtual FluxSolution parsimonious(const NetworkSpec& network,
                                       const std::string& objectiveReactionId,
                                       const std::vector<FluxBound>& bounds) const = 0;
    virtual std::vector<FluxRange> deletions(const NetworkSpec& network,
                                              const std::string& objectiveReactionId,
                                              const std::vector<FluxBound>& bounds,
                                              int order) const = 0;
};

// Pathway over-representation and network topology.
//
// enrich() REQUIRES a background set: the hypergeometric answer is a function of
// it, so there is no defaulted background to get wrong.
class IEnrichmentModule {
public:
    virtual ~IEnrichmentModule() = default;

    virtual EnrichmentReport enrich(const std::vector<std::string>& query,
                                    const std::vector<std::string>& background,
                                    const std::string& gmtPack) const = 0;
    virtual GraphMetrics graph(const std::vector<NetworkEdge>& edges) const = 0;
};

}  // namespace biocad
