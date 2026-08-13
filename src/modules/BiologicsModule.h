// modules/BiologicsModule.h - the antibody / protein-biologics adapter.
//
// bio/Imgt.*, bio/Liabilities.* and bio/Interface.* know numbering, motif rules,
// composition arithmetic and interface geometry. This adapter is the one place that
// composes them into the contract the UI and the agent code against, and it is
// deliberately thin: every refusal (a failed anchor, a T-cell receptor, a scheme with
// no published table, a TAP metric with no structure origin) already happens one
// layer down, and this file must never soften one of them into a default.
//
// The one judgement it does make is the heavy/light split for a mass ladder, because
// IBiologicsModule::massLadder() receives a flat chain list: a chain that numbers as
// VH or VHH is treated as heavy, everything else as light, and that is stated in the
// ladder's assumptions rather than assumed silently.
#pragma once

#include <string>
#include <vector>

#include "contracts/IModules.h"

namespace biocad {

class RealBiologics final : public IBiologicsModule {
public:
    AbDomain number(const std::string& sequence, NumberingScheme scheme) const override;
    AbDomain convertScheme(const AbDomain& domain, NumberingScheme to) const override;
    std::vector<SequenceLiability> liabilities(const AbDomain& domain,
                                               const bio::Structure* structure) const override;
    DevelopabilityReport developability(const std::vector<std::string>& chains,
                                       const bio::Structure* structure) const override;
    MassLadder massLadder(const std::vector<std::string>& chains,
                          int disulfideCount) const override;
    PeptideMap digest(const std::string& chain, const std::string& protease,
                      int maxMissedCleavages) const override;
    InterfaceReport interfaceOf(const bio::Structure& complex, const std::string& chainsA,
                                const std::string& chainsB) const override;
    AlanineScanReport alanineScan(const bio::Structure& complex, const std::string& chainsA,
                                  const std::string& chainsB) const override;
};

}  // namespace biocad
