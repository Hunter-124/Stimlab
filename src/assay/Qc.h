// assay/Qc.h - plate quality control.
//
// WHY every number here is a Quantity: a Z-prime computed from "whatever the
// extreme wells happened to be" is worse than no Z-prime, so a plate without both
// control roles returns notComputed naming the missing control instead. Detected
// edge, row and column effects are REPORTED and never corrected: median-polishing
// a plate in the QC step hides the pipetting fault that produced the pattern.
#pragma once

#include <cstddef>
#include <vector>

#include "data/Assay.h"

namespace biocad::assay {

// FINDING (reported to the plan owner, not worked around): data/Assay.h's QcReport
// has no field for the plain descriptive statistics or the blank-subtraction
// variance that Phase 10.2 asks for. Rather than smuggling them into
// QcReport::warnings as prose, they are returned by the two functions below and
// the QC panel can render them beside the report.
struct RoleStats {
    WellRole    role = WellRole::Unknown;
    std::size_t n = 0;             // included wells only
    Quantity    mean;
    Quantity    sd;
    Quantity    sem;
    Quantity    cvPct;             // notComputed unless the data are positive ratio-scale
    Quantity    median;
    Quantity    mad;               // raw median absolute deviation, not scaled
    Quantity    robustSd;          // 1.4826 * MAD, the sigma estimate
};

RoleStats roleStatistics(const Plate& plate, WellRole role);

// Var(x - mean(blank)) = s_x^2 + s_b^2 / n_b: subtracting a shared blank mean adds
// its own uncertainty to every well, which is the argument for more blanks.
Quantity blankSubtractedVariance(const Plate& plate);

QcReport plateQc(const Plate& plate);

// Two-sided Mann-Whitney U, normal approximation with tie correction and a
// continuity correction. Accurate to roughly 1e-3 in p for group sizes >= 8 and
// conservative below that; the exact permutation distribution is not implemented,
// which is stated here because a plate ring/interior split is always >= 20 vs 20.
double mannWhitneyP(const std::vector<double>& a, const std::vector<double>& b);

// Kruskal-Wallis H with tie correction; p from the chi-square survival function
// with k-1 degrees of freedom, evaluated as the regularized upper incomplete gamma
// Q(df/2, H/2) by the standard series/continued-fraction pair (relative accuracy
// ~1e-14, i.e. far better than a Wilson-Hilferty approximation would be).
double kruskalWallisP(const std::vector<std::vector<double>>& groups);

// Chi-square upper tail, exposed because both the H statistic above and the
// tests want it.
double chiSquareSurvival(double x, double degreesOfFreedom);

}  // namespace biocad::assay
