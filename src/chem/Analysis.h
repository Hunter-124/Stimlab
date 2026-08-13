// chem/Analysis.h - functional-group perception + Morgan fingerprint similarity.
//
// Group perception is DATA, not code: every flag below is one SMARTS pattern in
// assets/packs/rules/functional-groups.json, keyed by the member name. A group
// definition can therefore be read, cited and corrected without a rebuild, and
// the definitions are testable in isolation instead of being buried in bespoke
// graph walks.
#pragma once

#include <cstdint>
#include <string>
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
    bool catecholamine = false;   // phenethylamine core on a 3,4-dihydroxyphenyl ring
    bool anilide = false;         // N carrying both an arene and an acyl group (NAPQI-type)
    bool maoLabileAmine = false;  // primary/secondary amine on an unsubstituted alpha CH2
};

FunctionalGroups detectGroups(const Molecule& m);

// Rules that failed to load: a missing pack, an unknown group key, or a SMARTS
// that would not parse, each naming itself. A flag left false because its rule
// is broken is indistinguishable from a flag that is honestly false, so the
// breakage has to be surfaced rather than swallowed. Non-empty means some flag
// in FunctionalGroups is permanently false in this run.
const std::vector<std::string>& groupPackErrors();

// Circular (Morgan/ECFP-like) fingerprint -> sorted unique 32-bit features.
std::vector<std::uint32_t> morganFingerprint(const Molecule& m, int radius = 2);

// Tanimoto coefficient over two sorted unique feature sets (0..1).
double tanimoto(const std::vector<std::uint32_t>& a, const std::vector<std::uint32_t>& b);

}  // namespace biocad::chem
