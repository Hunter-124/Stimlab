// ui/Panels.h - feature panels rendered into the central Workspace window.
// Each takes the AppShell (for Services, UiState and the highlight bridge).
#pragma once

#include "data/Domain.h"

namespace biocad {

class AppShell;

namespace panels {

// Renders "label  value +/- error unit   (tier - source)" in the provenance
// colour. The error bar and the tier are part of the value, never a tooltip:
// a number whose provenance is hidden is a number that lies.
void drawQuantity(const char* label, const Quantity& q);

void dashboard(AppShell& shell);
void structureWorkbench(AppShell& shell);
void moleculeInput(AppShell& shell);
void compare(AppShell& shell);
void analogExplorer(AppShell& shell);
void stability(AppShell& shell);
void absorption(AppShell& shell);
void metabolism(AppShell& shell);
void alerts(AppShell& shell);
void metabolites(AppShell& shell);
void pkpd(AppShell& shell);
void popPk(AppShell& shell);
void interactionScenarios(AppShell& shell);
void ionization(AppShell& shell);
void assayWorkbench(AppShell& shell);
void assayDesign(AppShell& shell);
void sequenceCompare(AppShell& shell);
void proteinStructure(AppShell& shell);
void nucleicAcid(AppShell& shell);
void antibody(AppShell& shell);
void networks(AppShell& shell);
void flux(AppShell& shell);
void enrichment(AppShell& shell);
void similarity(AppShell& shell);
void legal(AppShell& shell);
void docking(AppShell& shell);
void workflows(AppShell& shell);
void library(AppShell& shell);
void runs(AppShell& shell);
void presets(AppShell& shell);
void settings(AppShell& shell);

}  // namespace panels
}  // namespace biocad
