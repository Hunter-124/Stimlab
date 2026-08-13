// modules/VariantModule.h - protein variant analysis wired to IVariantModule.
//
// bio/Conservation.* owns the statistics and bio/Rotamer.* owns the library, the
// geometry and the dead-end elimination. This adapter is the one place those turn
// into DTOs, and its whole job is to refuse in the right places:
//
//  - Below bio::kMinimumHomologs supplied homologs the profile is not usable and
//    EVERY quantity in the resulting VariantScore is notComputed() naming the
//    shortfall by number. The BLOSUM62 delta is the single exception: it is an
//    exact table lookup that does not depend on the homolog set at all.
//  - A rebuild with no rotamer pack on disk returns no chi angles and says which
//    file was missing. It never falls back to a guessed angle.
//  - StabilityPrediction is NOT produced here. See docs/variants.md: a ddG needs a
//    model actually running, and no model weights ship in this tree.
#pragma once

#include <string>
#include <vector>

#include "contracts/IModules.h"

namespace biocad {

class RealVariants final : public IVariantModule {
public:
    ConservationProfile conservation(const std::string& query,
                                     const std::vector<std::string>& homologs) const override;
    VariantScore score(const ConservationProfile& profile, int position,
                       char mutant) const override;
    RotamerRebuild rebuild(const bio::Structure& s, const std::string& chain, int residueNumber,
                           char mutant) const override;
};

// The one-letter code -> three-letter name mapping the rebuild path needs, and its
// inverse. 'X' / "UNK" for anything not one of the 20 standard residues.
[[nodiscard]] std::string threeLetterOf(char oneLetter);

// Neighbourhood radius for the repack, in angstrom, measured CB-to-CB (CA for
// glycine). Stated in every RotamerRebuild's assumptions rather than buried.
inline constexpr double kRepackRadiusAngstrom = 6.0;

// The repack is capped at this many neighbours, nearest first, so a buried
// position in a dense core cannot turn one mutation into an unbounded search.
inline constexpr int kMaxRepackedNeighbours = 8;

}  // namespace biocad
