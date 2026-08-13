#pragma once

// Structure quality and geometry scores that are EXACT over the coordinates they are
// given. Nothing here predicts anything, so every Quantity is Provenance::Measured.
//
// lDDT (Mariani et al. 2013, Bioinformatics 29:2722) is superposition-free: it compares
// the set of interatomic distances in the reference with the same distances in the model
// and never fits one onto the other. That is why it is the right tool for LOCAL quality
// and why it is insensitive to domain motions - a hinge that swings a whole subdomain
// wrecks a global RMSD but barely touches lDDT, because the distances inside each domain
// are unchanged. It also means lDDT(X, X) is EXACTLY 1.0 by construction: every distance
// difference is identically zero, so every one of the four tolerance tests passes, and the
// mean of four exact ones is one - not 0.999...
//
// SASA is Shrake & Rupley 1973 (J Mol Biol 79:351) numerical quadrature. A bare SASA
// number is not reproducible: two tools disagree by ~10% purely from probe radius, point
// count, radii set and hydrogen policy. Every Quantity produced here therefore states all
// four in its `source` string.

#include <cstddef>
#include <string>
#include <vector>

#include "bio/Structure.h"
#include "data/Domain.h"

namespace biocad::bio {

// ---------------------------------------------------------------------------
// Residue pairing
// ---------------------------------------------------------------------------

// Residue identity for cross-structure comparison. The insertion code is part of the
// identity, not decoration: 100 and 100A are different residues in the same chain.
struct ResidueKey {
    std::string chainId;
    int         authSeqId = 0;
    char        insertionCode = ' ';

    [[nodiscard]] bool operator==(const ResidueKey& o) const {
        return authSeqId == o.authSeqId && insertionCode == o.insertionCode
               && chainId == o.chainId;
    }
    [[nodiscard]] std::string label() const;   // "A:100A" - always author numbering
};

struct ResiduePair {
    ResidueKey     key;
    const Residue* model = nullptr;
    const Residue* reference = nullptr;
};

// The result of pairing two models residue by residue. The unmatched counts exist so a
// caller can REFUSE a comparison that matched almost nothing rather than reporting a
// confident score computed over three residues.
struct ResiduePairing {
    std::vector<ResiduePair> pairs;
    std::size_t unmatchedModel = 0;
    std::size_t unmatchedReference = 0;
};

// Pairs residues by (chain id, authSeqId, insertionCode). Author numbering is what papers
// and both file formats agree on; label numbering exists only in mmCIF.
[[nodiscard]] ResiduePairing pairResidues(const Model& model, const Model& reference);

// ---------------------------------------------------------------------------
// lDDT
// ---------------------------------------------------------------------------

struct LddtOptions {
    double inclusionRadius = 15.0;               // R0: reference pairs farther apart are ignored
    std::vector<double> tolerances{0.5, 1.0, 2.0, 4.0};   // angstrom
    bool includeHydrogens = false;               // most crystal structures have none
    bool includeHetatm = false;                  // waters and ligands are not the protein
    // lDDT compares like with like. When true (the default) a differing total atom count is
    // refused outright instead of quietly scoring whatever happened to match.
    bool requireEqualAtomCounts = true;
};

struct ResidueLddt {
    ResidueKey  key;
    std::size_t distancePairs = 0;   // reference pairs this residue contributed to
    double      score = 0.0;         // mean preserved fraction over the four tolerances
};

struct LddtResult {
    biocad::Quantity            global;             // NotComputed when a prerequisite is missing
    std::vector<ResidueLddt>  perResidue;
    std::size_t               consideredPairs = 0;
    std::size_t               unmatchedModelResidues = 0;
    std::size_t               unmatchedReferenceResidues = 0;
};

// Model 1 of each structure is used. Atoms are corresponded by name inside paired residues.
[[nodiscard]] LddtResult lddt(const Structure& model, const Structure& reference,
                              const LddtOptions& opts = {});

// ---------------------------------------------------------------------------
// SASA (Shrake-Rupley)
// ---------------------------------------------------------------------------

struct SasaOptions {
    double probeRadius = 1.4;      // angstrom, the conventional water probe
    int    points = 92;            // golden-spiral test points per atom
    // Default is united-atom: hydrogens are IGNORED, because the overwhelming majority of
    // crystal structures do not have them and including them for some inputs and not others
    // would make two runs incomparable.
    bool   includeHydrogens = false;
    bool   includeHetatm = false;  // waters would bury the protein surface
};

struct ResidueSasa {
    ResidueKey     key;
    std::string    residueName;
    biocad::Quantity absolute;   // A^2
    // Fraction of the Tien et al. 2013 (PLoS ONE 8:e80635) theoretical maximum for this
    // residue type. NotComputed for anything not one of the 20 standard residues - there is
    // no published maximum for a ligand or a modified residue.
    biocad::Quantity relative;
};

struct SasaResult {
    biocad::Quantity           total;      // A^2, NotComputed when there are no atoms
    std::vector<ResidueSasa> perResidue;
    std::string              method;     // the same text embedded in every source string
};

[[nodiscard]] SasaResult sasa(const Structure& s, const SasaOptions& opts = {});

// van der Waals radius for an element symbol (case-insensitive). Unknown elements get a
// documented 1.80 A fallback rather than silently contributing zero area.
[[nodiscard]] double vdwRadius(const std::string& element);

// Tien et al. 2013 theoretical maximum accessibility, A^2. Returns 0 for non-standard names.
[[nodiscard]] double maxAccessibility(const std::string& residueName);

}  // namespace biocad::bio
