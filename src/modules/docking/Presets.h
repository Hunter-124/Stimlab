// modules/docking/Presets.h - built-in CNS receptor preset table (Phase D, WP-3).
//
// The 29 CNS docking targets relevant to stimulant pharmacology: monoamine
// transporters, aminergic GPCRs, trace-amine receptor, opioid/sigma/cannabinoid
// receptors, the NMDA ionotropic receptor, the major metabolic enzymes, and the
// hERG cardiac-safety channel. Each entry carries a real PDB reference and a
// binding-site box (Vina center+size convention) sized to the co-crystal ligand
// pocket. This C++ table is the SINGLE SOURCE OF TRUTH so there is no runtime
// file IO and no path-resolution failure mode; src/presets/targets.yaml mirrors
// it for documentation only.
//
// SAFETY SCOPE: a receptor box defines WHERE a ligand binds for a binding-affinity
// (target-engagement / pharmacology) prediction. It is never a synthesis signal.
#pragma once

#include <string>
#include <vector>

#include "contracts/IDockingBackend.h"

namespace biocad::docking {

// The full built-in table of 29 CNS receptor presets. Deterministic; never empty.
const std::vector<ReceptorTarget>& cnsPresets();

// Resolve a preset by display name OR stable id (case-insensitive on id), e.g.
// "DAT", "DAT (dopamine transporter)", or "dat". Returns nullptr if unknown.
const ReceptorTarget* findPreset(const std::string& nameOrId);

// Display names of every preset, in table order (drives the UI dropdown). The
// historical five names (DAT/NET/SERT/TAAR1/hERG) appear first as a subset.
std::vector<std::string> presetNames();

// The headline CNS transporters/receptor to prepare first during provisioning
// (DAT/NET/SERT/TAAR1): the highest-value stimulant targets. Receptor prep for the
// remaining presets happens on demand; this keeps first-launch provisioning bounded.
std::vector<ReceptorTarget> headlinePresets();

}  // namespace biocad::docking
