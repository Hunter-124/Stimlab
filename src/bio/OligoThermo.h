#pragma once

// Nearest-neighbour duplex thermodynamics and secondary-structure scans.
//
// Everything numeric here comes out of assets/packs/nucleic/nn-thermodynamics.json
// (the SantaLucia unified parameter set), including the gas constant, the salt
// coefficient and the initiation terms, so a Tm is reproducible from the pack
// alone. Three conventions are made explicit rather than buried:
//
//   * Tm = dH / (dS + R ln(Ct/x)) - 273.15, where Ct is the TOTAL strand
//     concentration and x is 1 for a self-complementary strand and 4 for
//     non-self-complementary strands at equal concentration. Getting x wrong
//     shifts a self-complementary oligo's Tm by several degrees, so both
//     branches are exercised by the tests.
//   * Only monovalent salt is corrected, by the published entropy-per-phosphate
//     formula. A supplied Mg2+ or dNTP concentration is echoed on the result and
//     reported as UNUSED, because no divalent equivalence is transcribed in the
//     pack and inventing one would be worse than not applying it.
//   * Secondary-structure dG37 values are 1 M Na+ standard-state values with no
//     salt correction and no loop-closure penalty; the pack carries no loop
//     table, so a hairpin number here is the stem's dG37 and is labelled as such.

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include <nlohmann/json.hpp>

#include "data/Nucleic.h"

namespace biocad::bio {

struct ThermoOptions {
    double naMolar = 0.05;       // monovalent cation, mol/L
    double mgMolar = 0.0;        // echoed and reported unused
    double oligoMolar = 5e-8;    // TOTAL strand concentration Ct, mol/L
    double dntpMolar = 0.0;      // echoed and reported unused
};

class NnParameters {
public:
    // dH in kcal/mol and dS in cal/(mol*K) for a 5'->3' dinucleotide step.
    struct Term { double dH = 0, dS = 0; };

    const Term* step(char a, char b) const;
    Term        initiation(char terminalBase) const;   // G/C or A/T end term
    Term        symmetry() const { return symmetry_; }
    double      gasConstant() const { return gasConstant_; }
    double      saltCoefficient() const { return saltCoefficient_; }
    double      referenceTemperatureK() const { return referenceK_; }
    const std::string& source() const { return source_; }
    const std::string& saltMethod() const { return saltMethod_; }
    const std::string& magnesiumNote() const { return magnesiumNote_; }
    bool        empty() const { return steps_.empty(); }

    // dG at the pack's reference temperature (37 C) for one step.
    static double deltaG(const Term& t, double temperatureK) {
        return t.dH - temperatureK * t.dS / 1000.0;
    }

    friend NnParameters parseNnParameters(const nlohmann::json& j);

private:
    struct Entry { char a, b; Term term; };
    std::vector<Entry> steps_;
    Term        initGC_{}, initAT_{}, symmetry_{};
    double      gasConstant_ = 1.9872;
    double      saltCoefficient_ = 0.368;
    double      referenceK_ = 310.15;
    std::string source_, saltMethod_, magnesiumNote_;
};

NnParameters parseNnParameters(const nlohmann::json& j);
NnParameters loadNnParameters(const std::filesystem::path& file);
const NnParameters& builtinNnParameters();

// dH, dS, dG37, Tm and the conditions they were computed under. An ambiguity code
// anywhere in the sequence yields notComputed quantities naming the symbol - there
// are no nearest-neighbour parameters for a degenerate base, and averaging over
// its expansions would be an invention.
OligoThermo oligoThermo(std::string_view seq, const NnParameters& params,
                        const ThermoOptions& options);
OligoThermo oligoThermo(std::string_view seq, const ThermoOptions& options = {});

// True when the sequence is its own reverse complement, i.e. the x = 1 branch of
// the Tm equation applies.
bool isSelfComplementary(std::string_view seq);

// O(n^2) scans over all register shifts, scored with the same nearest-neighbour
// table. Only structures at or below `maxDeltaG37` (kcal/mol, so a negative
// cutoff asks for stable ones) are returned, most stable first.
std::vector<SecondaryStructure> hairpins(std::string_view seq, const NnParameters& params,
                                         double maxDeltaG37 = 0.0, int minLoop = 3);
std::vector<SecondaryStructure> selfDimers(std::string_view seq, const NnParameters& params,
                                           double maxDeltaG37 = 0.0);
std::vector<SecondaryStructure> heteroDimers(std::string_view a, std::string_view b,
                                             const NnParameters& params,
                                             double maxDeltaG37 = 0.0);
std::vector<SecondaryStructure> hairpins(std::string_view seq, double maxDeltaG37 = 0.0);
std::vector<SecondaryStructure> selfDimers(std::string_view seq, double maxDeltaG37 = 0.0);
std::vector<SecondaryStructure> heteroDimers(std::string_view a, std::string_view b,
                                             double maxDeltaG37 = 0.0);

}  // namespace biocad::bio
