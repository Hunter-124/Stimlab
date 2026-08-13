// modules/BioModules.h - the protein core wired to the module contracts.
//
// The adapters here are the ONLY place bio:: geometry and alignment results are
// turned into domain Quantities. The free functions are inline and header-only on
// purpose: they are the single implementation of that mapping, so there is exactly
// one place where a unit, a provenance tier or a source string gets decided.
//
// Nothing here accepts a chem::Conformer. A small-molecule distance-geometry
// embedding is not a protein, and bio::Structure has no constructor from one, so
// feeding an embedded ligand into lDDT or RMSD is a compile error rather than a
// number that looks like a structure comparison.
#pragma once

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "bio/Align.h"
#include "bio/CifReader.h"
#include "bio/PdbReader.h"
#include "bio/Score.h"
#include "bio/Structure.h"
#include "bio/Superpose.h"
#include "contracts/IModules.h"
#include "core/Error.h"
#include "data/Domain.h"
#include "packs/Pack.h"

namespace biocad {
namespace bioadapt {

// How the matrix and the gap costs are described everywhere a number derived from
// them is rendered. An alignment score without its matrix and gap costs is not
// reproducible, so the string travels with every Quantity.
inline const char* kAlignSource() {
    return "Gotoh affine gaps, BLOSUM62 (half-bit units), gap open 11 extend 1";
}

inline bio::GapCost defaultGaps() { return bio::GapCost{11, 1}; }

// The BLOSUM62 matrix, loaded once from assets/packs/matrices/blosum62.json.
// Returns nullptr when the data file is missing or malformed - a missing matrix is
// a NotComputed alignment naming the file, never a silently substituted identity
// matrix, which would produce plausible-looking and wrong percent identities.
inline const bio::SubstitutionMatrix* blosum62() {
    struct Loaded {
        bool ok = false;
        bio::SubstitutionMatrix m;
        Loaded() {
            const auto dir = packs::builtinPackDir();
            if (dir.empty()) return;
            try {
                m = bio::loadSubstitutionMatrix(dir / "matrices" / "blosum62.json");
                ok = true;
            } catch (const Error&) {
            } catch (const std::exception&) {
            }
        }
    };
    static const Loaded loaded;
    return loaded.ok ? &loaded.m : nullptr;
}

// Fills the string/stat fields shared by the global and the local path.
inline void fillRows(SequenceAlignment& out, const bio::AlignedRows& rows,
                     const bio::AlignmentStats& st, int score) {
    out.aligned1 = rows.a;
    out.aligned2 = rows.b;
    out.midline = rows.midline;
    out.gapOpens = static_cast<int>(st.gapOpens);
    out.alignedLength = static_cast<int>(st.alignedColumns);
    // Exact combinatorics over the two inputs: nothing is predicted, so Measured.
    out.score = makeQuantity(static_cast<double>(score), "half-bits", 0.0,
                             Provenance::Measured, kAlignSource());
    out.identityPct = makeQuantity(st.percentIdentity, "%", 0.0, Provenance::Measured,
                                   "identical columns / aligned columns");
    out.similarityPct = makeQuantity(st.percentSimilarity, "%", 0.0, Provenance::Measured,
                                     "positive-scoring columns / aligned columns");
}

inline SequenceAlignment noMatrix() {
    SequenceAlignment out;
    out.score = notComputed("assets/packs/matrices/blosum62.json");
    out.identityPct = notComputed("assets/packs/matrices/blosum62.json");
    out.similarityPct = notComputed("assets/packs/matrices/blosum62.json");
    out.eValue = notComputed("assets/packs/matrices/blosum62.json");
    out.note = "The BLOSUM62 substitution matrix could not be loaded, so no alignment was "
               "attempted. Nothing is guessed in its place.";
    return out;
}

// Needleman-Wunsch, end to end. The E-value stays NotComputed by construction:
// Karlin-Altschul statistics describe the best LOCAL hit expected by chance and
// have no global analogue, so a global E-value would be a significance claim with
// no null distribution behind it.
inline SequenceAlignment alignGlobal(const std::string& a, const std::string& b) {
    const bio::SubstitutionMatrix* m = blosum62();
    if (!m) return noMatrix();
    SequenceAlignment out;
    const bio::GlobalAlignment g = bio::alignGlobal(a, b, *m, defaultGaps());
    fillRows(out, g.rows, g.stats, g.score);
    out.eValue = notComputed("a local alignment (E-values are undefined for global alignments)");
    out.note = "Needleman-Wunsch global alignment: every residue of both sequences is placed, "
               "including the terminal gaps that a local alignment would trim.";
    return out;
}

// Smith-Waterman best subsegment, with a Karlin-Altschul E-value when the matrix
// file publishes statistics for these gap costs. The "database" is the second
// sequence alone, which is what makes the E-value honest here: it is the chance of
// this hit in a one-sequence search, not in a real database sweep.
inline SequenceAlignment alignLocal(const std::string& a, const std::string& b) {
    const bio::SubstitutionMatrix* m = blosum62();
    if (!m) return noMatrix();
    SequenceAlignment out;
    const bio::LocalAlignment l = bio::alignLocal(a, b, *m, defaultGaps());
    fillRows(out, l.rows, l.stats, l.score);

    const bio::KarlinAltschul* ka = m->statisticsFor(defaultGaps().open, defaultGaps().extend);
    if (!ka) {
        out.eValue = notComputed("published Karlin-Altschul lambda/K for BLOSUM62 at gap 11/1");
    } else {
        const bio::Significance sig =
            bio::evalueOf(l, *ka, a.size(), static_cast<double>(b.size()), 1);
        // Dimensionless by definition (an expected count), hence the empty unit.
        out.eValue = makeQuantity(sig.evalue, "", 0.0, Provenance::Measured,
                                  "Karlin-Altschul, bit score " +
                                      std::to_string(static_cast<long long>(sig.bitScore)) +
                                      ", search space is the second sequence only (1 seq, " +
                                      std::to_string(b.size()) + " residues)");
    }
    out.note = "Smith-Waterman local alignment: only the best-scoring subsegment is reported, "
               "spanning query residues " + std::to_string(l.aBegin + 1) + "-" +
               std::to_string(l.aEnd) + " and subject residues " + std::to_string(l.bBegin + 1) +
               "-" + std::to_string(l.bEnd) + " (1-based, inclusive).";
    return out;
}

// Reads one local structure file. The extension selects the reader because the two
// formats are not distinguishable cheaply and guessing wrong yields an empty
// structure rather than an error. std::nullopt means "could not be read at all";
// recoverable problems arrive in Structure::warnings and the structure still opens.
inline std::optional<bio::Structure> loadStructureFile(const std::filesystem::path& file) {
    std::error_code ec;
    if (!std::filesystem::is_regular_file(file, ec)) return std::nullopt;
    std::string ext = file.extension().string();
    for (auto& c : ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    try {
        if (ext == ".pdb" || ext == ".ent") return bio::readPdbFile(file);
        if (ext == ".cif" || ext == ".mmcif") return bio::readCifFile(file);
    } catch (const Error&) {
        return std::nullopt;
    } catch (const std::exception&) {
        return std::nullopt;
    }
    return std::nullopt;
}

// Alpha carbons of residues present in BOTH structures, paired by (chain, author
// number, insertion code). CA-only is the conventional backbone superposition set
// and it sidesteps side-chain naming differences that would silently drop atoms.
inline std::size_t collectCaPairs(const bio::Model& refModel, const bio::Model& modelModel,
                                  std::vector<bio::Point3>& mobile,
                                  std::vector<bio::Point3>& reference,
                                  std::size_t& unmatched) {
    const bio::ResiduePairing pairing = bio::pairResidues(modelModel, refModel);
    unmatched = pairing.unmatchedModel + pairing.unmatchedReference;
    for (const auto& p : pairing.pairs) {
        if (!p.model || !p.reference) continue;
        const bio::Atom* ma = p.model->atom(" CA ");
        const bio::Atom* ra = p.reference->atom(" CA ");
        if (!ma || !ra) continue;
        mobile.push_back({ma->x, ma->y, ma->z});
        reference.push_back({ra->x, ra->y, ra->z});
    }
    return mobile.size();
}

inline StructureComparison compare(const bio::Structure& ref, const bio::Structure& model) {
    StructureComparison out;
    // Phase 5.8 vendors the reference TM-align implementation; until then this is
    // the landing site for it. An approximated TM-score is worse than none, because
    // the number is read as if the reference implementation produced it.
    out.tmScore = notComputed("vendored TM-align (Phase 5.8)");

    const bio::Model* refModel = ref.model(1);
    const bio::Model* modModel = model.model(1);
    if (!refModel || !modModel) {
        out.rmsd = notComputed("model 1 in both structures");
        out.lddt = notComputed("model 1 in both structures");
        out.alignedResidues = notComputed("model 1 in both structures");
        out.note = "One of the two structures has no model 1, so nothing was compared.";
        return out;
    }

    std::vector<bio::Point3> mobile, reference;
    std::size_t unmatched = 0;
    const std::size_t paired = collectCaPairs(*refModel, *modModel, mobile, reference, unmatched);
    out.unmatchedResidues = static_cast<int>(unmatched);
    out.alignedResidues = makeQuantity(static_cast<double>(paired), "residues", 0.0,
                                       Provenance::Measured,
                                       "CA atoms paired by chain + author number + insertion code");

    if (paired < 3) {
        out.rmsd = notComputed("at least 3 paired CA atoms (a rotation is undetermined below that)");
    } else {
        try {
            const bio::Superposition s = bio::kabsch(mobile, reference);
            out.rmsd = makeQuantity(s.rmsd, "A", 0.0, Provenance::Measured,
                                    "Kabsch superposition over " + std::to_string(s.pairs) +
                                        " CA pairs" +
                                        (s.reflectionCorrected ? ", reflection corrected" : ""));
        } catch (const Error& e) {
            out.rmsd = notComputed(e.message);
        }
    }

    const bio::LddtResult l = bio::lddt(model, ref);
    out.lddt = l.global;

    out.note = "RMSD is CA-only and superposition-dependent; lDDT is superposition-free and "
               "therefore insensitive to domain motions. Residues are matched on AUTHOR "
               "numbering, which is the numbering both file formats agree on.";
    return out;
}

inline Quantity sasa(const bio::Structure& s) {
    return bio::sasa(s).total;
}

}  // namespace bioadapt

// Pairwise sequence alignment over the packaged BLOSUM62 matrix.
class RealSequence final : public ISequenceModule {
public:
    SequenceAlignment alignGlobal(const std::string& a, const std::string& b) const override;
    SequenceAlignment alignLocal(const std::string& a, const std::string& b) const override;
};

// Structure I/O and comparison. `load` reads a local file - including one the
// docking cache or a previous download already placed on disk. It never fetches.
class RealStructure final : public IStructureModule {
public:
    std::optional<bio::Structure> load(const std::filesystem::path& file) const override;
    StructureComparison compare(const bio::Structure& ref,
                                const bio::Structure& model) const override;
    Quantity sasa(const bio::Structure& s) const override;
};

}  // namespace biocad
