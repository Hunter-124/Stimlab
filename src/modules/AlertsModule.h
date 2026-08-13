// modules/AlertsModule.h - the structural-alert screen wired to its contract.
//
// This adapter is the ONLY place a chem::AlertHit becomes an AlertFlag, so there
// is exactly one place where an alert's severity is decided. It maps the pack's
// two severities onto Verdict::Warn and Verdict::Info and nothing else:
// Verdict::Danger is unreachable from here because a substructure match cannot
// support a toxicity verdict. The summary text carries the same framing, so a
// caller that renders only the summary still says "liability flag", not "toxic".
#pragma once

#include "contracts/IModules.h"

namespace biocad {

class RealAlerts final : public IAlertsModule {
public:
    AlertReport screen(const Molecule& m) const override;
};

}  // namespace biocad
