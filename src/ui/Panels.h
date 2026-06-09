// ui/Panels.h - feature panels rendered into the central Workspace window.
// Each takes the AppShell (for Services, UiState and the highlight bridge).
#pragma once

namespace stimlab {

class AppShell;

namespace panels {

void dashboard(AppShell& shell);
void structureWorkbench(AppShell& shell);
void moleculeInput(AppShell& shell);
void compare(AppShell& shell);
void analogExplorer(AppShell& shell);
void stability(AppShell& shell);
void absorption(AppShell& shell);
void metabolism(AppShell& shell);
void similarity(AppShell& shell);
void legal(AppShell& shell);
void docking(AppShell& shell);
void library(AppShell& shell);
void runs(AppShell& shell);
void presets(AppShell& shell);
void settings(AppShell& shell);

}  // namespace panels
}  // namespace stimlab
