// bio/Liabilities.h - sequence liabilities, developability descriptors, biologics
// mass ladders and peptide mapping.
//
// WHAT A LIABILITY FLAG IS. A motif match with a citation, and nothing else. The
// RATE of deamidation, isomerization, oxidation or clipping depends on local
// structure, formulation pH, excipients, light exposure and temperature, none of
// which a sequence contains. So a flag says "this site is the one the literature
// watches", never "this molecule degrades". The deamidation rules are ordered
// NG > NS > NT because that ordering is published; no rate is attached to it.
//
// WHY EXPOSURE IS A SEPARATE BIT. Weighting a Met or Trp oxidation flag by how
// solvent-exposed the side chain is only means something with coordinates. Without
// a structure, SequenceLiability::exposureKnown stays false and relativeSasa stays
// zero - it does NOT default to "exposed", which would turn every buried Met into a
// finding.
//
// WHY EVERY MASS COMES FROM chem::parseFormula. Composition arithmetic is exact and
// the isotope masses are NIST SRD 144 measurements, so a mass computed here is
// Provenance::Measured. Writing 145531.5 as a literal would be a number with no
// derivation and no way to correct it; instead the residue COMPOSITIONS are the
// data and every displayed mass is monoisotopicMass()/averageMass() over them. The
// same applies to the deltas: the disulfide -2.0157 is two hydrogens removed from
// the formula, and the 0.9840 / 1.0034 pair behind the required resolving power is
// (O - NH) and (13C - 12C) evaluated against the same table.
#pragma once

#include <optional>
#include <string>
#include <vector>

#include "bio/Structure.h"
#include "data/Biologics.h"

namespace biocad::bio {

// ---------------------------------------------------------------------------
// Liabilities
// ---------------------------------------------------------------------------

struct LiabilityRule {
    std::string ruleId;
    std::string motif;
    std::string label;
    std::string pattern;          // ECMAScript regex over the one-letter sequence
    int         risk = 0;         // rank index inside a motif family, never a rate
    bool        requiresExposure = false;
    bool        oddCountOnly = false;
    std::string citation;
};

struct LiabilityPack {
    bool                       ok = false;
    std::vector<LiabilityRule> rules;
    double                     exposureThreshold = 0.20;
    std::string                exposureNote;
    std::vector<std::string>   errors;
};

const LiabilityPack& liabilityPack();

// Scans a numbered domain. `structure`/`chainId` are optional; when they are absent
// (or the chain is not in the structure) every flag reports exposureKnown == false.
std::vector<SequenceLiability> scanLiabilities(const AbDomain& domain,
                                               const Structure* structure = nullptr,
                                               const std::string& chainId = {});

// ---------------------------------------------------------------------------
// Developability descriptors
// ---------------------------------------------------------------------------

// The amino-acid scales, loaded from assets/packs/biologics/descriptors.json.
struct DescriptorPack {
    bool ok = false;
    std::vector<std::string> errors;
};
const DescriptorPack& descriptorPack();

// Net charge at one pH from the Bjellqvist pKa set, sign included.
double netCharge(const std::string& sequence, double ph);

// (pH, charge) samples over [2, 12] in 0.1 steps, for the charge curve.
std::vector<std::pair<double, double>> netChargeCurve(const std::string& sequence);

struct DevelopabilityInput {
    // One entry per chain. Two chains are read as VH then VL for the Fv charge
    // symmetry parameter; any other count leaves that one NotComputed.
    std::vector<std::string> chains;
    const Structure*         structure = nullptr;
    // "crystal 1N8Z" or "homology model, <protocol>". EMPTY means the TAP metrics
    // are NotComputed: a TAP number without its structure origin invites a
    // homology-model threshold to be applied to a crystal structure.
    std::string              structureOrigin;
    std::vector<std::string> chainIds;   // parallel to `chains`, for the structure
};

DevelopabilityReport developability(const DevelopabilityInput& in);

// ---------------------------------------------------------------------------
// Mass ladders
// ---------------------------------------------------------------------------

// A glycan as a monosaccharide count, so its formula is composed rather than
// looked up: HexNAc/Hex/Fuc/NeuAc residue formulas plus the peptide.
struct Glycan {
    std::string name;
    int hexNAc = 0, hex = 0, fuc = 0, neuAc = 0;
};

// Common N-glycoforms of a therapeutic IgG, per heavy chain.
std::vector<Glycan> defaultGlycoforms();

struct MassLadderInput {
    std::vector<std::string> heavyChains;
    std::vector<std::string> lightChains;
    int  interchainDisulfides = 0;    // H-H and H-L bonds
    int  intrachainDisulfides = 0;    // domain bonds, counted once per chain set
    bool includeGlycoforms = true;
    bool pyroGlutamate = false;       // N-terminal Gln -> pyroGlu on every heavy chain
    bool cTerminalLysClipped = false;
    std::vector<Glycan> glycoforms;   // empty -> defaultGlycoforms()
};

MassLadder massLadder(const MassLadderInput& in);

// The Hill formula of a polypeptide, from the residue compositions plus one water.
// Empty when the sequence contains a residue with no composition (B, Z, J, X, U, O).
std::string peptideFormula(const std::string& sequence);

// ---------------------------------------------------------------------------
// Peptide mapping
// ---------------------------------------------------------------------------

// A variable modification expressed as formula arithmetic, never as a mass literal.
struct VariableMod {
    std::string name;
    std::string addFormula;      // e.g. "O" for oxidation
    std::string removeFormula;   // e.g. "NH" for deamidation (net +0.98)
    std::string residues;        // one-letter codes it applies to; empty = N-term
};

struct DigestOptions {
    std::string protease = "trypsin";   // trypsin | lysc | gluc | aspn | chymotrypsin
    int         maxMissedCleavages = 2;
    int         minLength = 4;
    bool        bAndYIons = true;
    std::vector<VariableMod> mods;
};

PeptideMap digest(const std::string& chain, const DigestOptions& opts);

// The proteases this build knows, for the UI and the agent tool schema.
std::vector<std::string> proteaseNames();

}  // namespace biocad::bio
