// chem/Analysis.h - functional-group perception + Morgan fingerprint similarity.
#pragma once

#include <cstdint>
#include <vector>

#include "chem/Molecule.h"

namespace biocad::chem {

// Structural liabilities / motifs the pharmacology heuristics key on.
struct FunctionalGroups {
    bool aromaticRing = false;
    bool phenol = false;
    bool catechol = false;        // two adjacent ring -OH
    bool ester = false;
    bool carboxylicAcid = false;
    bool amide = false;
    bool ketone = false;
    bool arylKetone = false;      // beta-keto / acylarene
    bool aldehyde = false;
    bool ether = false;
    bool methylenedioxy = false;  // -O-CH2-O- bridging an arene
    bool primaryAmine = false;
    bool secondaryAmine = false;
    bool tertiaryAmine = false;
    bool basicAmine = false;      // any non-amide sp3 amine
    bool nitrile = false;
    bool nitro = false;
    bool halogen = false;
    bool sulfoxide = false;
    bool sulfone = false;
    bool phenethylamine = false;  // Ar-C-C-N core
    bool catecholamine = false;   // catechol + phenethylamine
};

FunctionalGroups detectGroups(const Molecule& m);

// Circular (Morgan/ECFP-like) fingerprint -> sorted unique 32-bit features.
std::vector<std::uint32_t> morganFingerprint(const Molecule& m, int radius = 2);

// Tanimoto coefficient over two sorted unique feature sets (0..1).
double tanimoto(const std::vector<std::uint32_t>& a, const std::vector<std::uint32_t>& b);

}  // namespace biocad::chem
