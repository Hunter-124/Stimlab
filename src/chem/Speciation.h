// chem/Speciation.h - acid/base speciation: the general component-tableau
// equilibrium solver, and the independent-group microspecies ladder.
//
// WHY A TABLEAU AND NOT A CLOSED FORM. A monoprotic acid has an algebraic
// solution; a polyprotic acid in a buffer with a counter-ion does not, and the
// closed forms that are quoted for it are all approximations that silently fail
// at low concentration or near a pKa. So one solver handles every case: species
// i is a product of component activities,
//     a_i = K_i * prod_j x_j^(a_ij),      c_i = a_i / gamma_i
// mass balance on component j is
//     R_j = sum_i a_ij * c_i - T_j
// and Newton's method in ln x has the ANALYTIC Jacobian
//     J_jk = d R_j / d ln x_k = sum_i a_ij * a_ik * c_i
// which is symmetric positive semi-definite by construction (it is B^T C B with
// C diagonal positive), hence the LDLT first and a full-pivot LU only as the
// fallback for a rank-deficient tableau. Iterating in ln x is what keeps every
// concentration positive without a barrier or a clamp on c itself.
//
// CONDITIONS. Everything here is 25 C and thermodynamic-constant based: the
// caller's logK values are thermodynamic formation constants, and water's
// self-ionization is fixed at log Kw = -14.00 (25 C, zero ionic strength). There
// is no temperature argument because there is no enthalpy data here to
// extrapolate with; a constant at another temperature must be supplied as such
// by the caller rather than corrected by this file.
//
// A pKa IS NEVER PREDICTED HERE. Every pKa is an input - typed by the user or
// read from a cited pack - and carries its own Provenance. This file computes
// distributions from pKa values; it does not and will not estimate one from
// structure. A group whose pKa is NotComputed makes every quantity that depends
// on it NotComputed too, naming that group.
#pragma once

#include <vector>

#include "data/Ionization.h"

namespace biocad::chem {

// log10 of water's ion product at 25 C, zero ionic strength (Kw = 1.0e-14).
inline constexpr double kLogKwAt25C = -14.00;

// Damped-Newton solve of the tableau at fixed totals. `p.fixedComponent`, when
// >= 0, holds that component's log10 activity at `p.fixedLog10Activity` and
// drops its mass balance from the system - which is what "at pH x" means.
//
// Water is always present as chemistry, not as the caller's responsibility: if
// `p.species` contains no free-proton species (stoichiometry +1 on the proton
// component alone) and no hydroxide species (-1 on the proton component alone),
// they are APPENDED, in that order, after the caller's species. The returned
// vectors are parallel to that augmented list, so `concentrations.size()` may
// exceed `p.species.size()` by up to two. The proton component is the one named
// "H", "H+" or "H3O+" (case-sensitive on the first character), else
// `p.fixedComponent`.
SpeciationResult solveSpeciation(const SpeciationProblem& p);

// pH of a solution whose component totals are given, by an outer scalar solve on
// pH = -log10 a_H closing on charge balance: the inner solve holds the proton
// component fixed and satisfies every other mass balance exactly, so the only
// residual left for the outer solve is sum_i z_i c_i. Bracketed over pH 0..14 on
// a coarse scan and then bisected, because the charge-balance function is
// monotone in pH for real solutions but not smooth enough to trust a bare Newton
// step near a buffer plateau. Counter-ions must appear as components with their
// own totals (e.g. Cl in 0.1 M NH4Cl); otherwise the solution is not electrically
// neutral and there is no root to find.
SpeciationResult solveSpeciationPh(const SpeciationProblem& p);

// Microspecies ladder, net charge, logD and the isoelectric point over a pH
// range, for one molecule's `groups`.
//
// APPROXIMATION, stated because it is invisible in the output: the groups are
// treated as INDEPENDENT, so an N-group molecule's 2^N microstates are the
// product of N two-state distributions. The error is exactly the interaction
// between titrating groups - a proton already bound at one site shifts its
// neighbour's pKa (electrostatically, and by direct H-bonding for adjacent
// groups), which for a carboxyl/amine pair a couple of bonds apart is of order a
// pKa unit. Correcting it requires microscopic (site-specific) constants, which
// are inputs this signature does not have; when they exist, the tableau solver
// above is the thing to feed them to.
//
// `logP` is the input partition coefficient; logD = logP + log10(f_neutral),
// where f_neutral sums every microstate of net charge zero. Its Provenance is
// inherited (weakest of logP and every pKa), never asserted.
SpeciationCurve titrateGroups(const std::vector<IonizableGroup>& groups,
                              const Quantity& logP, double pHmin, double pHmax,
                              double step);

// Citation for the water constant and the Davies expression, for Quantity::source.
const char* speciationCitation();

}  // namespace biocad::chem
