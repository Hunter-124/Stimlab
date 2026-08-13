#pragma once

// Codon usage, CAI, and constraint-based codon optimization.
//
// Two things this file will not do, and cannot be made to do by the shape of its
// types: it does not predict expression level, protein yield or titre, and it does
// not emit a sequence that fails one of its own hard constraints. What it does is
// solve a constraint-satisfaction problem exactly - maximise CAI over synonymous
// codon choices such that the translation is preserved character for character
// and no forbidden pattern occurs anywhere in the result, including across codon
// boundaries and on the complementary strand. If no such assignment exists it
// reports the violations and returns no sequence.
//
// The optimizer is a dynamic program over codon positions whose state is a
// position in an Aho-Corasick automaton built from the forbidden patterns, so the
// across-boundary constraint is exact rather than a post-hoc filter that "usually"
// catches things. A GC window, when requested, is handled by a Lagrangian
// relaxation on top of that DP: it is reported as satisfied or not, never faked.

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include <nlohmann/json.hpp>

#include "bio/NucSeq.h"
#include "data/Nucleic.h"

namespace biocad::bio {

class CodonUsageTable {
public:
    CodonUsageTable() = default;
    CodonUsageTable(std::string id, std::string name, std::string source, int geneticCodeId,
                    std::vector<CodonUsageEntry> entries);

    const std::string& id() const { return id_; }
    const std::string& name() const { return name_; }
    const std::string& source() const { return source_; }
    int geneticCodeId() const { return geneticCodeId_; }
    const std::vector<CodonUsageEntry>& entries() const { return entries_; }
    const CodonUsageEntry* find(std::string_view codon) const;

    // The most frequent codon of a synonymous family, i.e. the one whose w is 1.
    std::string mostFrequentCodon(char aminoAcid) const;

private:
    std::string                  id_, name_, source_;
    int                          geneticCodeId_ = 1;
    std::vector<CodonUsageEntry> entries_;
};

class CodonUsageSet {
public:
    const CodonUsageTable* find(std::string_view id) const;
    std::vector<std::string> ids() const;
    void add(CodonUsageTable t) { tables_.push_back(std::move(t)); }
    bool empty() const { return tables_.empty(); }

private:
    std::vector<CodonUsageTable> tables_;
};

CodonUsageSet parseCodonUsage(const nlohmann::json& j);
CodonUsageSet loadCodonUsage(const std::filesystem::path& file);
const CodonUsageSet& builtinCodonUsage();          // ids: "ecoli-k12", "hsapiens"
const CodonUsageTable& builtinCodonUsageTable(std::string_view id);

// CAI is the geometric mean of the relative adaptiveness w of the coding codons.
// CONVENTION, stated because implementations differ and the number is not
// comparable across them: Met (ATG) and Trp (TGG) are excluded because a
// single-codon family has w = 1 by construction and carries no adaptation signal,
// and stop codons are excluded because the family is not under the same
// selection. A codon with w = 0 in the table would make the geometric mean 0, so
// it is excluded and counted in `warnings` instead of silently zeroing the score.
CodonMetrics codonMetrics(std::string_view cds, const CodonUsageTable& table);

struct OptimizeOptions {
    int geneticCodeId = 1;
    // IUPAC patterns that must not occur. Restriction sites and any other motif;
    // the reverse complement of each is forbidden too, since a site on either
    // strand is a site.
    std::vector<std::string> forbiddenPatterns;
    int    maxHomopolymer = 6;      // 0 disables; expressed internally as four patterns
    bool   enforceGcWindow = false;
    double minGcPercent = 40.0;
    double maxGcPercent = 60.0;
    std::string usageTableId;       // recorded in the result's assumptions
};

// Solves the constraint problem. On success `translationPreserved` is true,
// `optimized` translates to exactly `protein` under the requested table, and
// `remainingViolations` is empty. On failure `optimized` is empty and
// `remainingViolations` names what could not be satisfied - the function never
// returns a sequence that violates the translation or a forbidden pattern.
CodonOptimizationResult optimizeCodons(std::string_view protein, const CodonUsageTable& table,
                                       const OptimizeOptions& options);

// Every occurrence of a forbidden pattern (or its reverse complement) in `seq`,
// as "PATTERN@index" strings. Used by the optimizer's own output check and
// directly testable.
std::vector<std::string> findForbidden(std::string_view seq,
                                       const std::vector<std::string>& patterns,
                                       int maxHomopolymer);

}  // namespace biocad::bio
