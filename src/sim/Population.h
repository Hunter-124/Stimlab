// sim/Population.h - the reproducible population layer over Phase 4's PK engine.
//
// This file adds VARIABILITY to a PK simulation and nothing else: every subject is
// integrated through pkpd::simulate, so there is exactly one PK model in BioCAD and
// a population band can never disagree with the typical-value profile drawn beside
// it about what the model is.
//
// The three layers are separately toggled because they answer different questions.
// Between-subject variability (eta ~ MVN(0, Omega), Pi = theta_i * exp(eta_i))
// asks how much people differ. Parameter uncertainty, sampled by Latin hypercube
// with Iman-Conover correlation induction from the fit covariance, asks how well
// the typical value itself is known. Residual error asks how noisy the assay is.
// A single merged band cannot answer any of the three.
//
// SAFETY SCOPE: the output is an exposure scenario with a percentile band. The band
// describes the variability that was ENTERED. It is not a prediction about an
// individual and there is no dose anywhere in this file.
#pragma once

#include "data/Population.h"

namespace biocad::sim {

// Simulate `variability.subjects` subjects and reduce them to percentile bands.
// The seed is part of the contract: the same VariabilitySpec produces a
// byte-identical PopulationProfile, and a different seed produces a different one.
PopulationProfile simulatePopulation(const PkModelSpec& model, const DoseRegimen& regimen,
                                     const VariabilitySpec& variability);

// The maximum number of individual trajectories kept for the faint-line overlay.
// Fifty is a rendering decision, not a statistical one: the band is computed from
// every subject, and this cap only limits what is drawn and serialised.
inline constexpr int kMaxStoredTrajectories = 50;

}  // namespace biocad::sim
