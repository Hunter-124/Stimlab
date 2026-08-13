// data/Biologics.h - JSON-serializable DTOs for antibody and protein-biologics
// analysis: numbering, regions, sequence liabilities, developability descriptors,
// mass ladders, peptide maps, interfaces and geometric scanning.
//
// SAFETY AND HONESTY SCOPE, permanent:
//  - "Closest germline set" is a sequence-similarity result. It is NEVER rendered
//    as species identification, and there is no species field here to render.
//  - There is no humanness score, no humanization suggestion, no immunogenicity
//    prediction, no affinity-maturation proposal, no aggregation free energy, no
//    viscosity, no expression titre and no shelf life. Those are the fields whose
//    absence is the feature.
//  - The alanine scan is GEOMETRIC: it truncates beyond C-beta and reports lost
//    buried area, contacts, hydrogen bonds and salt bridges. It is therefore
//    unit-free Provenance::Heuristic and is NEVER a delta-delta-G.
#pragma once

#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "data/Domain.h"

namespace biocad {

enum class NumberingScheme { Imgt, Kabat, Chothia, Martin, Aho };

NLOHMANN_JSON_SERIALIZE_ENUM(NumberingScheme, {
    {NumberingScheme::Imgt, "imgt"},
    {NumberingScheme::Kabat, "kabat"},
    {NumberingScheme::Chothia, "chothia"},
    {NumberingScheme::Martin, "martin"},
    {NumberingScheme::Aho, "aho"},
})

enum class AbChainType { Unknown, HeavyVh, LightVKappa, LightVLambda, Vhh, TcrBeta, TcrAlpha };

NLOHMANN_JSON_SERIALIZE_ENUM(AbChainType, {
    {AbChainType::Unknown, "unknown"},
    {AbChainType::HeavyVh, "VH"},
    {AbChainType::LightVKappa, "VK"},
    {AbChainType::LightVLambda, "VL"},
    {AbChainType::Vhh, "VHH"},
    {AbChainType::TcrBeta, "TRB"},
    {AbChainType::TcrAlpha, "TRA"},
})

// One numbered residue. `insertionCode` carries IMGT's insertion letters so a
// long CDR3 numbers without renumbering its neighbours.
struct NumberedResidue {
    int         position = 0;
    std::string insertionCode;
    char        aminoAcid = 'X';
    int         sequenceIndex = 0;   // 0-based index into the input sequence
    std::string region;              // FR1, CDR1, FR2, ...
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(NumberedResidue, position, insertionCode, aminoAcid,
                                   sequenceIndex, region)

// A numbered antibody domain. Anchor failures are FIRST-CLASS: IMGT 23 Cys,
// 41 Trp, 89 hydrophobic, 104 Cys and 118 Phe/Trp must all be satisfied, and a
// failure returns NotComputed numbering rather than plausible-looking wrong
// numbering, which is worse than none.
struct AbDomain {
    std::string                  sequence;
    AbChainType                  chain = AbChainType::Unknown;
    NumberingScheme              scheme = NumberingScheme::Imgt;
    std::vector<NumberedResidue> residues;
    std::string                  closestGermlineSet;   // NOT a species identification
    std::string                  runnerUpGermlineSet;
    double                       bestBitScore = 0;
    double                       runnerUpBitScore = 0;
    std::string                  vGene;
    std::string                  jGene;
    std::vector<int>             cdrLengths;           // CDR1, CDR2, CDR3
    std::vector<std::string>     anchorFailures;       // empty when numbering is valid
    bool                         numbered = false;
    std::vector<std::string>     warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(AbDomain, sequence, chain, scheme, residues,
                                   closestGermlineSet, runnerUpGermlineSet, bestBitScore,
                                   runnerUpBitScore, vGene, jGene, cdrLengths, anchorFailures,
                                   numbered, warnings)

// One pack-defined, cited sequence liability. `relativeSasa` is only populated
// when a structure was supplied; otherwise `exposureKnown` is false and the flag
// says exposure is unknown rather than assuming the site is exposed.
struct SequenceLiability {
    std::string ruleId;
    std::string motif;         // N-X-S/T, NG, DG, DP, free C...
    std::string label;
    int         sequenceIndex = 0;
    std::string imgtPosition;
    std::string region;
    bool        exposureKnown = false;
    double      relativeSasa = 0;
    std::string citation;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SequenceLiability, ruleId, motif, label, sequenceIndex,
                                   imgtPosition, region, exposureKnown, relativeSasa, citation)

// Sequence-derived developability descriptors. Every one is exact arithmetic on
// the sequence, hence Measured - which is precisely why none of them is a
// developability VERDICT.
struct DevelopabilityReport {
    Quantity                 isoelectricPoint;
    Quantity                 netChargeAtPh74;
    Quantity                 extinctionCoefficient280;   // M^-1 cm^-1
    Quantity                 gravy;
    Quantity                 aliphaticIndex;
    Quantity                 instabilityIndex;
    Quantity                 fvChargeSymmetry;
    // TAP metrics are only meaningful with the structure origin and modelling
    // protocol attached, so they are NotComputed unless `structureOrigin` is set.
    Quantity                 tapPsh;
    Quantity                 tapPpc;
    Quantity                 tapPnc;
    Quantity                 tapSfvcsp;
    std::string              structureOrigin;    // "crystal 1IGT" | "homology model, <protocol>"
    std::vector<std::string> assumptions;
    std::vector<std::string> warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(DevelopabilityReport, isoelectricPoint, netChargeAtPh74,
                                   extinctionCoefficient280, gravy, aliphaticIndex,
                                   instabilityIndex, fvChargeSymmetry, tapPsh, tapPpc, tapPnc,
                                   tapSfvcsp, structureOrigin, assumptions, warnings)

// ---------------------------------------------------------------------------
// Mass and peptide mapping.
// ---------------------------------------------------------------------------

struct MassLadderEntry {
    std::string species;        // "intact", "reduced HC", "Fab", "G0F/G0F"...
    Quantity    average;        // Da - the meaningful one above ~10 kDa
    Quantity    monoisotopic;   // Da
    int         disulfides = 0;
    std::string note;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MassLadderEntry, species, average, monoisotopic, disulfides,
                                   note)

struct MassLadder {
    std::vector<MassLadderEntry> entries;
    // The resolving power needed to separate deamidation (+0.984016 Da) from the
    // 13C isotope spacing (1.003355 Da) at the reported mass. Showing it is how a
    // mass difference stops being mistaken for an identification.
    Quantity                     requiredResolvingPower;
    std::vector<std::string>     assumptions;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MassLadder, entries, requiredResolvingPower, assumptions)

struct PeptideFragment {
    std::string sequence;
    int         begin = 0;
    int         end = 0;
    int         missedCleavages = 0;
    Quantity    monoisotopic;
    Quantity    average;
    std::vector<std::string> modifications;
    std::vector<double>      bIons;
    std::vector<double>      yIons;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PeptideFragment, sequence, begin, end, missedCleavages,
                                   monoisotopic, average, modifications, bIons, yIons)

struct PeptideMap {
    std::string                  protease;
    std::vector<PeptideFragment> peptides;
    double                       coveragePct = 0;
    std::vector<std::string>     warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PeptideMap, protease, peptides, coveragePct, warnings)

// ---------------------------------------------------------------------------
// Interfaces and geometric scanning.
// ---------------------------------------------------------------------------

struct ResidueContact {
    std::string chainA;
    std::string residueA;
    std::string chainB;
    std::string residueB;
    double      minDistance = 0;    // Angstrom
    int         atomContacts = 0;   // pairs within the 4.5 A cutoff
    bool        hydrogenBond = false;
    bool        saltBridge = false;
    bool        hydrophobic = false;
    bool        piStacking = false;
    bool        cationPi = false;
    bool        disulfide = false;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ResidueContact, chainA, residueA, chainB, residueB,
                                   minDistance, atomContacts, hydrogenBond, saltBridge,
                                   hydrophobic, piStacking, cationPi, disulfide)

// Levy's support/core/rim classification, which needs both the isolated and the
// complexed SASA - so the report carries both rather than a single burial number.
struct InterfaceReport {
    std::string                 chainsA;
    std::string                 chainsB;
    Quantity                    buriedSurfaceArea;   // A^2, SASA_A + SASA_B - SASA_AB
    Quantity                    sasaA;
    Quantity                    sasaB;
    Quantity                    sasaComplex;
    std::vector<ResidueContact> contacts;
    std::vector<std::string>    supportResidues;
    std::vector<std::string>    coreResidues;
    std::vector<std::string>    rimResidues;
    std::vector<std::string>    cdrContacts;    // antibody complexes only
    std::vector<std::string>    epitope;
    std::vector<std::string>    paratope;
    std::string                 sasaParameters;  // algorithm, probe, points, radii, H policy
    std::vector<std::string>    warnings;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(InterfaceReport, chainsA, chainsB, buriedSurfaceArea, sasaA,
                                   sasaB, sasaComplex, contacts, supportResidues, coreResidues,
                                   rimResidues, cdrContacts, epitope, paratope, sasaParameters,
                                   warnings)

// One geometric alanine-scan position. `impact` is unit-free Heuristic: it is a
// rank-ordering of how much interface each side chain contributes, and calling it
// an energy would be a lie with a plausible number attached.
struct AlanineScanPosition {
    std::string chain;
    std::string residue;
    double      lostBuriedAreaA2 = 0;
    int         lostContacts = 0;
    int         lostHydrogenBonds = 0;
    int         lostSaltBridges = 0;
    Quantity    impact;    // Heuristic, no unit, never kcal/mol
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(AlanineScanPosition, chain, residue, lostBuriedAreaA2,
                                   lostContacts, lostHydrogenBonds, lostSaltBridges, impact)

// The scan, WITH its calibration. Reporting the Spearman correlation against a
// measured benchmark subset is what separates "a geometric proxy that ranks
// hotspots this well" from "an energy".
struct AlanineScanReport {
    std::vector<AlanineScanPosition> positions;   // descending impact
    Quantity                         benchmarkSpearman;
    std::string                      benchmarkName;   // e.g. the SKEMPI antibody/antigen subset
    std::string                      disclaimer;      // the explicit not-a-ddG sentence
    std::vector<std::string>         assumptions;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(AlanineScanReport, positions, benchmarkSpearman,
                                   benchmarkName, disclaimer, assumptions)

}  // namespace biocad
