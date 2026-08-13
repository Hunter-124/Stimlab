// assay/Design.h - 10.4, prospective assay design by forward simulation.
//
// WHY THIS FILE EXISTS. Every other file in src/assay/ analyses data that already
// exists. This one answers the question a chemist asks BEFORE pipetting: with the
// noise I actually have, will this plate layout recover the EC50 I care about, and
// will the confidence interval I report mean what it says? That second half is the
// point. A fitter whose nominal 95% interval covers the truth 70% of the time is
// not reporting a 95% interval, and the only way to find that out is to simulate a
// known truth many times and count.
//
// WHY IT REUSES THE PRODUCTION PATH. Every simulated plate is serialised to the
// same long CSV the instrument export is turned into, then pushed through
// assay::importText -> assay::plateQc -> assay::fitSeries. A design simulator with
// a private analysis path measures the private path, not the product; it would
// stay green while the real importer silently mangled a column. There is exactly
// one analysis path in BioCAD and this file is a client of it.
//
// SAFETY SCOPE: a concentration ladder for a plate is not a dose for a person, and
// a dilution volume is not a synthesis procedure. Nothing here selects a reagent,
// a route, or a condition.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "data/Assay.h"

namespace biocad::assay {

// PCG64-DXSM, implemented here rather than taken from <random>.
//
// WHY NOT std::mt19937_64: the engine itself is specified bit-exactly, but
// std::normal_distribution and std::lognormal_distribution are NOT - the standard
// fixes their distribution, not their algorithm, so the same seed yields different
// draws on a different standard library. A design report whose coverage number
// changes with the toolchain is not reproducible, and reproducibility is the whole
// claim being made. So the bit stream (PCG64-DXSM, O'Neill 2014 with the 2019 DXSM
// output permutation: 128-bit LCG state, cheap-multiplier double-xorshift-multiply
// output) and the normal transform (Box-Muller on two 53-bit uniforms) both live
// in BioCAD and are byte-identical everywhere.
class DesignRng {
public:
    explicit DesignRng(std::uint64_t seed);

    std::uint64_t next();            // raw 64 bits
    double        uniform();         // [0,1), 53 significant bits
    double        normal();          // standard normal, Box-Muller, cached pair
    double        lognormal(double sigma);  // exp(sigma * normal()), median 1

private:
    __uint128_t state_ = 0;
    __uint128_t inc_ = 0;
    double      spare_ = 0.0;
    bool        hasSpare_ = false;
};

// Evaluates the stated truth model at one concentration (mol/L).
//
// Parameter order is the letter order of the 5PL in data/Assay.h,
// y = D + (A - D) / [1 + (x/C)^B]^G, so that a 4PL is exactly a 5PL with G = 1:
//   4PL:                  {A, B, C, D}          A = signal at x -> 0, D = plateau
//   5PL:                  {A, B, C, D, G}
//   MichaelisMenten:      {Vmax, Km}
//   Hill:                 {Vmax, K, n}
//   SubstrateInhibition:  {Vmax, Km, Ki}
// Returns false when the model is not one this simulator can generate from, which
// is honest: the biophysics models need a trace, not a concentration ladder.
bool truthValue(AssayModel model, const std::vector<double>& p, double concentration,
                double& out);

// The concentration the model's half-maximal point sits at, for the models that
// have one (4PL C, 5PL C*(2^(1/G)-1)^(1/B), Hill K, MM Km). Returns 0 when the
// model has none.
double truthEc50(AssayModel model, const std::vector<double>& p);

// One simulated plate as the long CSV an import would receive. The layout is
// column 1 = negative control (the fully-affected plateau D), column 2 = positive
// control (the untreated signal A), remaining wells = the concentration series,
// row-major, `replicates` series per plate. Controls are generated first so that a
// spec with no room for both controls fails QC loudly instead of silently
// producing a plate whose Z-prime is NotComputed for a reason nobody can see.
std::string simulatePlateCsv(const AssayDesignSpec& spec, std::uint64_t runSeed,
                             const std::string& plateId);

// Fedorov coordinate-exchange D-optimal concentration placement, constrained to a
// serial-dilution ladder.
//
// The candidate set is NOT a continuous interval: a bench ladder can only reach
// top * f^(-j/2) for integer j, i.e. the entered ladder plus the one intermediate
// point a second half-step dilution can make. Optimising over a continuum would
// return concentrations nobody can pipette.
struct DOptimalLadder {
    std::vector<double> candidates;       // the achievable ladder points
    std::vector<double> concentrations;   // the chosen design, ascending
    double              dEfficiencyGain = 1.0;  // (det_opt/det_entered)^(1/p)
    std::vector<std::string> notes;
};
DOptimalLadder dOptimalLadder(const AssayDesignSpec& spec);

// The whole 10.4 measurement: `replicateRuns` seeded repetitions, each one a fresh
// plate through import -> QC -> fit, reporting what the design would recover.
// Empirical CI coverage is the headline: the fraction of runs whose OWN reported
// interval contained the truth.
AssayDesignReport simulateDesign(const AssayDesignSpec& spec);

// ---------------------------------------------------------------------------
// Dilution and mass arithmetic. Exact, not modelled: these are unit conversions
// with a molar mass, and they carry Measured provenance for that reason.
// ---------------------------------------------------------------------------

// Grams of solid needed for `volumeL` of a `molarity` mol/L stock.
double massForStock(double molarMassGPerMol, double molarity, double volumeL);

struct DilutionStep {
    int    index = 0;              // 0 = the top well, filled from the stock
    double concentration = 0.0;    // mol/L in the well after the transfer
    double transferVolumeUl = 0.0; // taken from the previous well (or the stock)
    double diluentVolumeUl = 0.0;  // added to reach `wellVolumeUl`
    double fold = 0.0;             // dilution factor achieved by this transfer
};

struct DilutionPlan {
    std::vector<DilutionStep> steps;
    double                    stockMolarity = 0.0;
    double                    wellVolumeUl = 0.0;
    // The compounded pipetting CV at the LAST step, sqrt(n) * cv for n transfers:
    // the reason a long ladder is less trustworthy at the bottom than the top.
    double                    compoundedCvAtLastStep = 0.0;
    std::vector<std::string>  warnings;
};

// A serial ladder from `stockMolarity` down, `fold` per transfer, `steps` wells of
// `wellVolumeUl` each. Warns when a transfer volume falls below `minTransferUl`,
// because that is where pipetting CV stops being the number you assumed.
DilutionPlan serialDilution(double stockMolarity, double topMolarity, double fold, int steps,
                            double wellVolumeUl, double minTransferUl = 2.0);

}  // namespace biocad::assay
