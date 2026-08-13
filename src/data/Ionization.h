// data/Ionization.h - JSON-serializable DTOs for exact chemistry: formulas and
// masses, acid/base speciation, buffers, and solubility/dissolution.
//
// SAFETY SCOPE: everything here is stoichiometry, equilibrium and dissolution
// physics for a stated composition. There is deliberately no reaction route, no
// precursor, no reagent-sourcing and no scale-up type: balancing an equation the
// user wrote is arithmetic, whereas proposing how to make something is out of
// scope by design. Nothing here emits a dose.
//
// PROVENANCE RULE: pKa, melting point, Ksp and precipitation kinetics are INPUTS
// - from the user or from a cited pack - and are never guessed. A curve whose
// prerequisite is absent is a NotComputed Quantity naming that prerequisite,
// never a curve computed from an assumed default.
#pragma once

#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "data/Domain.h"

namespace biocad {

// ---------------------------------------------------------------------------
// 11.1 - Formula and mass.
// ---------------------------------------------------------------------------

// One element of a parsed molecular formula. `count` is the total atom count
// after expanding any nesting, so C6H5(CH3) arrives as C7 H8.
struct FormulaElement {
    std::string symbol;
    int         count = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(FormulaElement, symbol, count)

// A parsed formula with the two masses that are NOT interchangeable: the
// monoisotopic mass (sum of the most abundant isotope of each element) and the
// average mass (sum of standard atomic weights). Above ~10 kDa only the average
// mass is meaningful for an observed envelope centroid, so both are always
// carried rather than one being labelled "the" mass.
struct FormulaMass {
    std::string                 formula;         // canonical Hill-order string
    std::vector<FormulaElement> elements;
    int                         charge = 0;      // net charge of the written species
    Quantity                    monoisotopic;    // Da
    Quantity                    average;         // Da
    Quantity                    mz;              // Da/charge, NotComputed when charge == 0
    int                         electrons = 0;   // total electron count
    double                      unsaturation = 0;  // rings + double bonds equivalent
    std::vector<std::string>    warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(FormulaMass, formula, elements, charge, monoisotopic,
                                   average, mz, electrons, unsaturation, warnings)

// One peak of a theoretical isotope envelope, relative intensity normalized so
// the most abundant peak is 1.0.
struct IsotopePeak {
    double mass = 0;        // Da
    double intensity = 0;   // 0..1, relative to the base peak
    int    nominalShift = 0;  // integer mass offset from the monoisotopic peak
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(IsotopePeak, mass, intensity, nominalShift)

struct IsotopeEnvelope {
    std::string              formula;
    std::vector<IsotopePeak> peaks;         // ascending mass
    double                   prunedBelow = 0;  // intensity threshold actually applied
    std::string              source;        // isotope table citation
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(IsotopeEnvelope, formula, peaks, prunedBelow, source)

// A balanced chemical equation. Coefficients come from the integer null space of
// the element-conservation matrix, so an unbalanceable equation is reported as
// such rather than approximated. This is stoichiometry only: no conditions, no
// reagent quantities beyond the user's own, no route.
struct BalancedEquation {
    std::vector<std::string> reactants;      // formula strings as entered
    std::vector<std::string> products;
    std::vector<int>         reactantCoefficients;
    std::vector<int>         productCoefficients;
    bool                     balanced = false;
    std::string              limitingReagent;   // empty when no amounts were given
    Quantity                 theoreticalYield;  // g, NotComputed without amounts
    Quantity                 atomEconomy;       // percent, NotComputed when unbalanced
    std::vector<std::string> warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(BalancedEquation, reactants, products,
                                   reactantCoefficients, productCoefficients, balanced,
                                   limitingReagent, theoreticalYield, atomEconomy, warnings)

// ---------------------------------------------------------------------------
// 11.2 - Ionization and speciation.
// ---------------------------------------------------------------------------

// One ionizable group. `pKa` is an input: `provenance` must be Measured (a cited
// value) or Model (a value the user typed for exploration); a Predicted pKa is
// only legitimate with a named model and its benchmark error.
struct IonizableGroup {
    std::string label;         // e.g. "carboxyl", "amine"
    Quantity    pKa;
    bool        acidic = true; // true: HA -> A- ; false: BH+ -> B
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(IonizableGroup, label, pKa, acidic)

// A component-tableau equilibrium problem: species i is formed from components j
// with stoichiometry `stoichiometry[i][j]` and formation constant log10 `logK[i]`.
// Components carry total concentrations `totals` (mol/L). H+ is a component like
// any other; fixing its activity is what "at pH x" means.
struct SpeciationProblem {
    std::vector<std::string>         components;    // component names, e.g. {"H", "A"}
    std::vector<double>              totals;        // mol/L per component
    std::vector<std::string>         species;       // species names
    std::vector<std::vector<double>> stoichiometry; // species x components
    std::vector<double>              logK;          // per species, log10 beta
    std::vector<double>              charges;       // per species, for charge balance
    int                              fixedComponent = -1;  // index held at fixedLog10Activity
    double                           fixedLog10Activity = 0;  // e.g. -pH for H
    double                           ionicStrength = 0;   // mol/L; 0 selects ideal activities
    bool                             daviesActivities = false;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SpeciationProblem, components, totals, species,
                                   stoichiometry, logK, charges, fixedComponent,
                                   fixedLog10Activity, ionicStrength, daviesActivities)

// The solved distribution. `massBalanceResidual` and `chargeBalanceResidual` are
// reported, not hidden: a solver that does not show its residual is asking to be
// trusted on faith. Convergence to <1e-10 relative is the acceptance criterion.
struct SpeciationResult {
    std::vector<double>      concentrations;   // mol/L per species
    std::vector<double>      fractions;        // 0..1 per species, of its own component
    std::vector<double>      componentFree;    // mol/L free per component
    double                   pH = 0;
    double                   netCharge = 0;    // equivalents/L
    double                   massBalanceResidual = 0;   // max relative residual
    double                   chargeBalanceResidual = 0;
    int                      iterations = 0;
    bool                     converged = false;
    double                   ionicStrength = 0;
    std::vector<double>      activityCoefficients;  // per species; all 1 when ideal
    std::vector<std::string> warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SpeciationResult, concentrations, fractions,
                                   componentFree, pH, netCharge, massBalanceResidual,
                                   chargeBalanceResidual, iterations, converged,
                                   ionicStrength, activityCoefficients, warnings)

// One pH point of a titration/distribution curve.
struct SpeciationPoint {
    double              pH = 0;
    std::vector<double> microspeciesFractions;  // parallel to SpeciationCurve::labels
    double              netCharge = 0;
    double              logD = 0;               // logP + log10(f_neutral)
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SpeciationPoint, pH, microspeciesFractions, netCharge, logD)

struct SpeciationCurve {
    std::string                  moleculeId;
    std::vector<std::string>     labels;      // microspecies labels
    std::vector<SpeciationPoint> points;      // ascending pH
    Quantity                     isoelectricPoint;   // pH of zero net charge
    Quantity                     logDAtPh74;
    Quantity                     logP;               // the input logP that logD used
    std::vector<std::string>     assumptions;
    std::vector<std::string>     warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SpeciationCurve, moleculeId, labels, points,
                                   isoelectricPoint, logDAtPh74, logP, assumptions, warnings)

// ---------------------------------------------------------------------------
// 11.3 - Buffers, solubility, dissolution.
// ---------------------------------------------------------------------------

// One buffer component: a conjugate pair at total concentration `totalMolar`.
struct BufferComponent {
    std::string label;
    double      pKa = 0;
    double      totalMolar = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(BufferComponent, label, pKa, totalMolar)

struct BufferCapacityPoint {
    double pH = 0;
    double beta = 0;   // mol/L per pH unit, Van Slyke
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(BufferCapacityPoint, pH, beta)

struct BufferReport {
    std::vector<BufferComponent>     components;
    std::vector<BufferCapacityPoint> curve;
    Quantity                         betaAtPh74;
    Quantity                         maxCapacity;
    Quantity                         maxCapacityPh;
    std::vector<std::string>         assumptions;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(BufferReport, components, curve, betaAtPh74,
                                   maxCapacity, maxCapacityPh, assumptions)

struct SolubilityPoint {
    double pH = 0;
    double logS = 0;      // log10 mol/L of TOTAL dissolved species
    bool   saltLimited = false;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SolubilityPoint, pH, logS, saltLimited)

// pH-solubility with the salt-limited branch. `intrinsic` is S0 - the neutral
// species' solubility - and is NotComputed unless a melting point was supplied
// (the General Solubility Equation needs one) or a measured S0 was entered.
struct SolubilityReport {
    std::string                  moleculeId;
    Quantity                     intrinsic;      // S0, mol/L
    Quantity                     pHmax;          // pH of the salt-solubility kink
    Quantity                     solubilityAtPh74;
    std::vector<SolubilityPoint> curve;
    Quantity                     doseNumber;         // Do, dimensionless (BCS)
    Quantity                     dissolutionNumber;  // Dn
    Quantity                     absorptionNumber;   // An
    std::vector<std::string>     assumptions;
    std::vector<std::string>     warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SolubilityReport, moleculeId, intrinsic, pHmax,
                                   solubilityAtPh74, curve, doseNumber, dissolutionNumber,
                                   absorptionNumber, assumptions, warnings)

// One step of a dissolution/precipitation time course. The three masses are
// carried separately so their sum can be asserted constant to 1e-9 every step:
// a dissolution model that loses mass is wrong in a way a single curve hides.
struct DissolutionPoint {
    double timeS = 0;
    double dissolvedMolar = 0;
    double solidMg = 0;
    double precipitatedMg = 0;
    double particleRadiusUm = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(DissolutionPoint, timeS, dissolvedMolar, solidMg,
                                   precipitatedMg, particleRadiusUm)

struct DissolutionReport {
    std::vector<DissolutionPoint> points;
    Quantity                      timeTo85Pct;   // s
    double                        maxMassImbalance = 0;  // mg, must stay < 1e-9
    std::vector<std::string>      assumptions;
    std::vector<std::string>      warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(DissolutionReport, points, timeTo85Pct,
                                   maxMassImbalance, assumptions, warnings)

// Everything the Ionization panel needs for one compound, so the panel makes one
// call and renders whatever came back - including the NotComputed branches.
struct IonizationReport {
    std::string       moleculeId;
    FormulaMass       mass;
    IsotopeEnvelope   envelope;
    SpeciationCurve   speciation;
    SolubilityReport  solubility;
    BufferReport      buffer;
    DissolutionReport dissolution;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(IonizationReport, moleculeId, mass, envelope,
                                   speciation, solubility, buffer, dissolution)

}  // namespace biocad
