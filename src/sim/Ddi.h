// sim/Ddi.h - mechanistic drug-interaction arithmetic.
//
// Three models, in increasing order of what they require and of what they claim:
//
//   1. The FDA basic-model R-values. R1 = 1 + [I]u/Ki (hepatic, reversible),
//      R1,gut = 1 + [I]g/Ki, R2 for time-dependent inactivation, and the induction
//      R. These are SCREENING numbers: they say whether an interaction can be ruled
//      out, not how large it is, and they need no victim parameters at all.
//   2. The mechanistic static AUCR, which combines reversible inhibition,
//      time-dependent inactivation and induction over the hepatic and gut terms.
//      This one needs fm, and without fm it is notComputed("fm") - assuming fm = 1
//      is precisely what turns a 1.3-fold interaction into a 5-fold one. Without Fg
//      only the hepatic term is reported and gutIncluded stays false.
//   3. The dynamic enzyme model, integrated through numeric::rk4Integrate. At a
//      constant inhibitor concentration its steady state must EQUAL the static
//      model, and EnzymeTimeCourse carries both numbers so the agreement is visible
//      rather than asserted in a comment.
//
// Every physiological constant (Qh, Qen, kdeg) comes from assets/packs/physiology.json
// through core::physiology(). There is no fallback: a missing pack yields
// NotComputed, because a built-in default is an uncited constant with better
// manners.
//
// SAFETY SCOPE: an AUCR is an exposure ratio. Nothing here converts one into a dose,
// a dose adjustment or a risk category, and the impairment scenarios are explicitly
// editable ratios rather than a Child-Pugh or creatinine-clearance formula. Whole-body
// PBPK is absent for the reason stated in data/Population.h.
#pragma once

#include "data/Population.h"

namespace biocad::sim {

// R-values plus the mechanistic static AUCR for one perpetrator against one victim.
InteractionReport interaction(const PerpetratorSpec& perpetrator, const VictimSpec& victim);

// Dynamic enzyme activity over `horizonH` hours at the perpetrator's stated constant
// unbound hepatic concentration:
//   dE/dt = kdeg*(1 + d*Emax*I/(EC50+I)) - kdeg*E - kinact*I/(KI+I)*E
//   CLint(t) = CLint0 * E / (1 + I/Ki)
EnzymeTimeCourse enzymeTimeCourse(const PerpetratorSpec& perpetrator, double horizonH);

// Renal and hepatic impairment as an explicit exposure ratio.
//   Renal:   AUC ratio = 1 / (1 - fe*(1 - renalFunctionRatio))
//   Hepatic: the well-stirred model re-evaluated at hepaticClintRatio * fu*CLint,
//            with Qh from the physiology pack.
ImpairmentScenario impairment(const VictimSpec& victim, double renalFunctionRatio,
                              double hepaticClintRatio);

}  // namespace biocad::sim
