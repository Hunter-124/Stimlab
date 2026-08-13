// Tests for plate QC. The point of every case is that a QC number is either right
// to the last digit or is notComputed naming what is missing - a Z-prime guessed
// from the extreme wells is worse than no Z-prime at all.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <string>
#include <vector>

#include "assay/Qc.h"

using namespace biocad;
using namespace biocad::assay;
using Catch::Matchers::WithinAbs;

namespace {

// Four positive and four negative control wells with hand-computable statistics:
//   pos = {100,102,98,100}: mean 100, sum of squared deviations 8, sd = sqrt(8/3)
//   neg = { 10, 12, 8, 10}: mean  10, identical sd
Plate controlPlate() {
    Plate p;
    p.id = "qc";
    p.rows = 8;
    p.columns = 12;
    p.readoutUnit = "RFU";
    const double pos[] = {100, 102, 98, 100};
    const double neg[] = {10, 12, 8, 10};
    for (int i = 0; i < 4; ++i) {
        Well w;
        w.plateId = "qc";
        w.row = i;
        w.role = WellRole::PositiveControl;
        w.column = 0;
        w.readout = pos[i];
        w.well = std::string(1, static_cast<char>('A' + i)) + "1";
        p.wells.push_back(w);
        w.role = WellRole::NegativeControl;
        w.column = 11;
        w.readout = neg[i];
        w.well = std::string(1, static_cast<char>('A' + i)) + "12";
        p.wells.push_back(w);
    }
    return p;
}

Plate patternPlate(bool plantEdgeEffect) {
    Plate p;
    p.id = "mw";
    p.rows = 8;
    p.columns = 12;
    for (int r = 0; r < 8; ++r) {
        for (int c = 0; c < 12; ++c) {
            Well w;
            w.role = WellRole::Sample;
            w.row = r;
            w.column = c;
            const bool ring = r == 0 || c == 0 || r == 7 || c == 11;
            w.readout = 100.0 + static_cast<double>((r * 12 + c) % 5 - 2) * 0.5 +
                        ((plantEdgeEffect && ring) ? 20.0 : 0.0);
            p.wells.push_back(w);
        }
    }
    return p;
}

}  // namespace

TEST_CASE("plateQc matches the hand-calculated Z-prime and SSMD", "[assay]") {
    const QcReport r = plateQc(controlPlate());
    const double sd = std::sqrt(8.0 / 3.0);
    // Z' = 1 - 3*(sd_p + sd_n)/|mean_p - mean_n| = 1 - 3*2*sqrt(8/3)/90
    //    = 1 - 9.797958971132712/90 = 0.8911337892096365
    const double zPrime = 1.0 - 3.0 * (2.0 * sd) / 90.0;
    // SSMD = (100-10)/sqrt(sd_p^2 + sd_n^2) = 90/sqrt(16/3) = 38.97114317029974
    const double ssmd = 90.0 / std::sqrt(16.0 / 3.0);

    REQUIRE(r.zPrime.provenance == Provenance::Measured);
    REQUIRE(r.zPrime.unit.empty());
    REQUIRE_THAT(r.zPrime.value, WithinAbs(zPrime, 1e-12));
    REQUIRE_THAT(r.ssmd.value, WithinAbs(ssmd, 1e-12));
    REQUIRE_THAT(r.signalToBackground.value, WithinAbs(10.0, 1e-12));
    REQUIRE_THAT(r.signalToNoise.value, WithinAbs(90.0 / sd, 1e-12));
    REQUIRE_THAT(r.positiveSd.value, WithinAbs(sd, 1e-12));
    REQUIRE_THAT(r.cvPositivePct.value, WithinAbs(100.0 * sd / 100.0, 1e-12));
    // The bands are stated rather than collapsed into a colour.
    REQUIRE(r.interpretation.find(">= 0.5 excellent") != std::string::npos);
    REQUIRE(r.interpretation.find("marginal") != std::string::npos);
    REQUIRE(r.interpretation.find("unusable") != std::string::npos);
}

TEST_CASE("plateQc refuses a Z-prime without both controls", "[assay]") {
    Plate p = controlPlate();
    std::vector<Well> kept;
    for (const Well& w : p.wells) {
        if (w.role != WellRole::PositiveControl) kept.push_back(w);
    }
    p.wells = kept;

    const QcReport r = plateQc(p);
    REQUIRE(r.zPrime.provenance == Provenance::NotComputed);
    REQUIRE(r.zPrime.source.find("positive control") != std::string::npos);
    REQUIRE(r.robustZPrime.provenance == Provenance::NotComputed);
    REQUIRE(r.ssmd.provenance == Provenance::NotComputed);
    REQUIRE(r.signalToNoise.provenance == Provenance::NotComputed);
    REQUIRE_FALSE(r.warnings.empty());
    // The negative controls are still summarised: only the comparison is refused.
    REQUIRE(r.negativeMean.provenance == Provenance::Measured);
}

TEST_CASE("%CV refuses signed data and MAD statistics do not", "[assay]") {
    Plate p = controlPlate();
    for (Well& w : p.wells) {
        if (w.role == WellRole::PositiveControl) w.readout -= 100.0;   // 0, 2, -2, 0
    }
    const QcReport r = plateQc(p);
    REQUIRE(r.cvPositivePct.provenance == Provenance::NotComputed);
    REQUIRE(r.cvPositivePct.source == "%CV requires positive ratio-scale data");

    const RoleStats stats = roleStatistics(p, WellRole::PositiveControl);
    REQUIRE(stats.n == 4);
    REQUIRE_THAT(stats.median.value, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(stats.mad.value, WithinAbs(1.0, 1e-12));         // |0|,|2|,|2|,|0| -> median 1
    REQUIRE_THAT(stats.sd.value, WithinAbs(std::sqrt(8.0 / 3.0), 1e-12));
    REQUIRE_THAT(stats.sem.value, WithinAbs(std::sqrt(8.0 / 3.0) / 2.0, 1e-12));
}

TEST_CASE("blankSubtractedVariance adds the blank mean's own uncertainty", "[assay]") {
    Plate p;
    p.id = "bl";
    p.rows = 1;
    p.columns = 6;
    const double samples[] = {100, 104, 96, 100};   // s_x^2 = 32/3
    for (int i = 0; i < 4; ++i) {
        Well w;
        w.role = WellRole::Sample;
        w.column = i;
        w.readout = samples[i];
        p.wells.push_back(w);
    }
    for (int i = 0; i < 2; ++i) {
        Well w;
        w.role = WellRole::Blank;
        w.column = 4 + i;
        w.readout = 2.0 + 2.0 * i;   // {2,4}: s_b^2 = 2, n_b = 2
        p.wells.push_back(w);
    }
    // 32/3 + 2/2 = 11.666666666666667
    const Quantity v = blankSubtractedVariance(p);
    REQUIRE(v.provenance == Provenance::Measured);
    REQUIRE_THAT(v.value, WithinAbs(32.0 / 3.0 + 1.0, 1e-12));

    Plate oneBlank = p;
    oneBlank.wells.pop_back();
    REQUIRE(blankSubtractedVariance(oneBlank).provenance == Provenance::NotComputed);
}

TEST_CASE("Mann-Whitney detects a planted edge effect and clears a flat plate", "[assay]") {
    const QcReport planted = plateQc(patternPlate(true));
    const QcReport flat = plateQc(patternPlate(false));

    REQUIRE(planted.edgeEffectP.provenance == Provenance::Measured);
    REQUIRE(planted.edgeEffectP.value < 0.01);
    REQUIRE(flat.edgeEffectP.value > 0.2);
    // Detected, reported, and NOT corrected: the readouts are untouched.
    bool warned = false;
    for (const std::string& w : planted.warnings) {
        if (w.find("edge effect detected") != std::string::npos) warned = true;
    }
    REQUIRE(warned);
    REQUIRE(flat.rowEffectP.value > 0.2);
    REQUIRE(flat.columnEffectP.value > 0.2);
}

TEST_CASE("Kruskal-Wallis finds a row gradient and the chi-square tail is exact", "[assay]") {
    Plate p;
    p.id = "kw";
    p.rows = 8;
    p.columns = 12;
    for (int r = 0; r < 8; ++r) {
        for (int c = 0; c < 12; ++c) {
            Well w;
            w.role = WellRole::Sample;
            w.row = r;
            w.column = c;
            w.readout = 100.0 + 10.0 * r + static_cast<double>((r * 12 + c) % 5 - 2) * 0.5;
            p.wells.push_back(w);
        }
    }
    const QcReport r = plateQc(p);
    REQUIRE(r.rowEffectP.value < 1e-6);
    REQUIRE(r.columnEffectP.value > 0.2);

    // Chi-square upper tail against published values: P(X^2_1 > 3.841459) = 0.05,
    // P(X^2_2 > 5.991465) = 0.05, P(X^2_10 > 18.307038) = 0.05.
    REQUIRE_THAT(chiSquareSurvival(3.8414588206941236, 1.0), WithinAbs(0.05, 1e-9));
    REQUIRE_THAT(chiSquareSurvival(5.9914645471079799, 2.0), WithinAbs(0.05, 1e-9));
    REQUIRE_THAT(chiSquareSurvival(18.307038053275146, 10.0), WithinAbs(0.05, 1e-9));
    REQUIRE(chiSquareSurvival(0.0, 3.0) == 1.0);

    // Mann-Whitney against a hand-checkable case: {1,2,3,4} vs {5,6,7,8} is the most
    // extreme rank split, U = 0.
    REQUIRE(mannWhitneyP({1, 2, 3, 4}, {5, 6, 7, 8}) < 0.05);
    REQUIRE(mannWhitneyP({1, 2, 3, 4}, {1, 2, 3, 4}) == 1.0);
}
