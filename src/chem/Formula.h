// chem/Formula.h - exact formula arithmetic: parsing, Hill notation, exact and
// average mass, isotope envelopes, formula finding, and equation balancing.
//
// DATA. Every mass and abundance is read from
// assets/packs/descriptors/isotopes.json, transcribed from NIST Standard
// Reference Database 144 (Atomic Weights and Isotopic Compositions). Nothing in
// this file hard-codes an isotope mass, so the table can be corrected without
// touching code. A missing or malformed pack yields ok == false with a warning
// naming the file - never a silently empty table, because a zero mass reads as a
// computed answer.
//
// WHY MEASURED. Isotope masses and terrestrial abundances are measured physical
// constants with published uncertainties, so a mass computed from them is
// Provenance::Measured sourced to the NIST table. That is not a prediction and
// must not be tiered as one.
//
// SAFETY SCOPE. Balancing an equation the user wrote, finding its limiting
// reagent, its theoretical yield and its atom economy is arithmetic on a stated
// composition. There is deliberately nothing here about how to run a reaction:
// no conditions, no route, no precursor selection, no scale-up.
#pragma once

#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "data/Ionization.h"

namespace biocad::chem {

// One term of a parsed formula. `massNumber == 0` means natural isotopic
// composition; a nonzero value is an explicit label ([13C], [2H]) and pins both
// the exact and the average mass of those atoms to that one isotope.
struct FormulaTerm {
    int z = 0;
    int massNumber = 0;
    int count = 0;
};

struct ParsedFormula {
    std::vector<FormulaTerm> terms;      // merged and sorted by (z, massNumber)
    int                      charge = 0; // net charge written on the formula
    std::vector<std::string> warnings;
    bool                     ok = false; // false => the numbers are meaningless
};

// Parses a formula string. Handles nested () [] {} with multipliers, hydrate
// dots (CuSO4.5H2O, CuSO4*5H2O), a leading multiplier (2H2O), isotope labels
// ([13C], 13C, [2H]) and a trailing charge (+, -, 2+, +2, ++). Returns
// std::nullopt only for input with no atoms at all; a recoverable oddity comes
// back as ok == false with the reason in `warnings` so the caller can show it.
std::optional<ParsedFormula> parseFormula(std::string_view text);

// Hill order: carbon first, hydrogen second, every other element alphabetically
// by symbol; labelled atoms are written [13C] and sort beside their element. The
// charge is appended as +, -, 2+, 3- when nonzero.
std::string hillFormula(const ParsedFormula& f);

// Sum of the most abundant isotope of each element (the labelled isotope where
// one was written). Zero with a warning if the table is unavailable.
double monoisotopicMass(const ParsedFormula& f);

// Sum of standard atomic weights. Labelled atoms contribute their isotope mass,
// because a labelled position is not an isotopic average.
double averageMass(const ParsedFormula& f);

// Total electron count: sum of Z minus the net charge.
int electronCount(const ParsedFormula& f);

// Rings-plus-double-bond equivalents, 1 + sum n_i (v_i - 2) / 2 over a table of
// conventional valences. Elements with no conventional valence (transition
// metals, noble gases) are excluded and named in the warning, because an RDBE
// that quietly ignores an atom is worse than no RDBE.
double ringPlusDoubleBondEquivalents(const ParsedFormula& f, std::vector<std::string>* warnings);

// The whole FormulaMass DTO. `mz` is notComputed("charge is zero") for a neutral
// species; otherwise it is (M - charge * m_e) / |charge|, i.e. the electron mass
// is removed for a cation and added for an anion, because an ion is the neutral
// molecule minus (or plus) that many electrons.
FormulaMass toFormulaMass(const ParsedFormula& f);
FormulaMass toFormulaMass(std::string_view text);

// Sparse convolution of the per-element isotope distributions, peaks aggregated
// by nominal mass shift with an abundance-weighted exact mass. Combinations
// below `minIntensity` relative to the running maximum are pruned during the
// convolution, so the cost stays linear in the surviving peak count. Peaks come
// back ascending in mass, normalized so the base peak is exactly 1.0.
IsotopeEnvelope isotopeEnvelope(const ParsedFormula& f, double minIntensity = 1e-6);

// Every formula inside `elementBounds` whose monoisotopic mass is within
// `toleranceDa` of `targetMass`. Bounds are {symbol -> {min, max}} and are the
// ONLY source of candidate elements: no element is guessed, and no valence,
// senior-rule or "chemically sensible" filter is applied here. Results ascend by
// absolute mass error.
std::vector<ParsedFormula> findFormulas(double targetMass, double toleranceDa,
                                        const std::map<std::string, std::pair<int, int>>& elementBounds);

// Balances the equation by taking the integer null space of the element-
// conservation matrix (Gaussian elimination over exact rationals, scaled by the
// LCM of the denominators and reduced by the GCD). A null space that is not
// one-dimensional, or a solution with a non-positive coefficient, sets
// balanced == false with a warning naming why: an equation is either balanced
// exactly or not at all, and least-squares "nearly balanced" coefficients would
// be a fabrication.
//
// With `reactantGrams` supplied (parallel to `reactants`, non-positive entries
// ignored) the limiting reagent is the reactant with the smallest
// moles/coefficient, the theoretical yield is that of the FIRST product, and
// atom economy is 100 * (mass of the first product) / (total reactant mass) at
// the balanced stoichiometry. Amounts only - nothing here recommends conditions,
// a route, or a precursor.
BalancedEquation balanceEquation(const std::vector<std::string>& reactants,
                                 const std::vector<std::string>& products,
                                 const std::vector<double>&      reactantGrams = {});

// The exact attribution string every Quantity built from the isotope table
// carries as its source.
const char* formulaCitation();

// Whether the isotope pack loaded, and the note (pack path, or why it failed)
// that a panel should surface.
bool        isotopeTableOk();
const char* isotopeTableNote();

}  // namespace biocad::chem
