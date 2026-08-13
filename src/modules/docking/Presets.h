// modules/docking/Presets.h - receptor presets, resolved from the loaded data packs.
//
// The target catalog is DATA, not code: it comes from assets/packs/*.json plus the
// user overlay under %APPDATA%/BioCAD/packs (see src/packs/Pack.h). This header is
// the docking module's read-only view of it, cached once per process so a preset
// pointer stays valid for the lifetime of the app.
//
// Only targets that carry a real PDB reference AND a real binding-site box are
// dockable presets. Pack entries without a box are coverage gaps, deliberately
// visible in the Presets panel rather than filled in with an invented site.
//
// SAFETY SCOPE: a receptor box defines WHERE a ligand binds for a binding-affinity
// (target-engagement / pharmacology) prediction. It is never a synthesis signal.
#pragma once

#include <string>
#include <vector>

#include "contracts/IDockingBackend.h"
#include "packs/Pack.h"

namespace biocad::docking {

// Every dockable preset from the loaded packs, in pack order. Cached on first use.
const std::vector<ReceptorTarget>& targetPresets();

// The full load report behind targetPresets(), including packs that failed to load.
// The Presets panel renders its errors; nothing else should need it.
const packs::LoadReport& targetPackReport();

// Re-read the packs from disk (used after the user edits their pack directory).
void reloadTargetPacks();

// Resolve a preset by display name OR stable id (case-insensitive on id), e.g.
// "DAT", "DAT (dopamine transporter)", or "dat". Returns nullptr if unknown.
const ReceptorTarget* findPreset(const std::string& nameOrId);

// Display names of every preset, in table order (drives the UI dropdown).
std::vector<std::string> presetNames();

// The receptor targets a pack marked `"headline": true`: the ones prepared first
// during provisioning. Receptor prep for the rest happens on demand, which keeps
// first-launch provisioning bounded.
std::vector<ReceptorTarget> headlinePresets();

}  // namespace biocad::docking
