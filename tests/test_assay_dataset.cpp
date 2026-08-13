// Tests for plate import, normalization and the opt-in outlier rules. Each case is
// one that would catch a real bug: a silently dropped column, a well name clamped
// into its neighbour, a normalization that quietly returns zeros, or a rule that
// deletes a point instead of flagging it.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "assay/Dataset.h"
#include "assay/Stats.h"

using namespace biocad;
using namespace biocad::assay;
using Catch::Matchers::WithinAbs;

namespace {

// A reader-style grid export: an empty label cell, then 1..cols, then one row per
// row letter. Values are base + 10*row + column so every cell is identifiable.
std::string gridBlock(int rows, int cols, double base) {
    std::string s;
    for (int c = 1; c <= cols; ++c) s += "," + std::to_string(c);
    s += "\n";
    for (int r = 0; r < rows; ++r) {
        std::string label;
        int n = r + 1;
        while (n > 0) {
            label.insert(label.begin(), static_cast<char>('A' + (n - 1) % 26));
            n = (n - 1) / 26;
        }
        s += label;
        for (int c = 0; c < cols; ++c) s += "," + std::to_string(base + 10.0 * r + c);
        s += "\n";
    }
    return s;
}

}  // namespace

TEST_CASE("importText reads a long table order-free and keeps unknown columns", "[assay]") {
    const char* csv =
        "Plate_ID,Well,ROLE,sample_id,series_id,Concentration,conc unit,replicate,readout,"
        "readout_unit,Operator\n"
        "p1,A1,positive,ctrl-max,,0,M,1,100,RFU,jane\n"
        "p1,A2,negative,vehicle,,0,M,1,10,RFU,jane\n"
        "p1,B1,sample,cpd-1,s1,1e-6,M,1,42.5,RFU,jane\n"
        "p1,B2,sample,cpd-1,s1,1e-7,M,2,61.25,RFU,jane\n"
        "p1,B3,blank,,,0,M,1,2,RFU,jane\n";
    std::string error;
    const AssayDataset ds = importText(csv, &error);

    REQUIRE(error.empty());
    REQUIRE(ds.detectedLayout == "long");
    REQUIRE(ds.plates.size() == 1);
    const Plate& p = ds.plates[0];
    REQUIRE(p.id == "p1");
    REQUIRE(p.wells.size() == 5);
    // Deduced from the observed labels A1..B3, not from the well count.
    REQUIRE(p.rows == 2);
    REQUIRE(p.columns == 3);
    REQUIRE(p.wells[0].role == WellRole::PositiveControl);
    REQUIRE(p.wells[1].role == WellRole::NegativeControl);
    REQUIRE(p.wells[2].role == WellRole::Sample);
    REQUIRE(p.wells[4].role == WellRole::Blank);
    REQUIRE(p.wells[2].concentration == 1e-6);
    REQUIRE(p.wells[3].concentration == 1e-7);
    REQUIRE(p.wells[3].replicate == 2);
    REQUIRE(p.wells[2].readoutUnit == "RFU");
    // The instrument's own column survives instead of being dropped.
    REQUIRE(ds.metadata.size() == 1);
    REQUIRE(ds.metadata[0].first == "Operator");
    REQUIRE(ds.metadata[0].second == "jane");
}

TEST_CASE("importText melts a 96-well grid block", "[assay]") {
    std::string error;
    const AssayDataset ds = importText(gridBlock(8, 12, 100.0), &error);

    REQUIRE(error.empty());
    REQUIRE(ds.detectedLayout == "96-well grid (8x12)");
    REQUIRE(ds.plates.size() == 1);
    REQUIRE(ds.plates[0].wells.size() == 96);
    REQUIRE(ds.plates[0].rows == 8);
    REQUIRE(ds.plates[0].columns == 12);
    const Well& a1 = ds.plates[0].wells.front();
    const Well& h12 = ds.plates[0].wells.back();
    REQUIRE(a1.well == "A1");
    REQUIRE(a1.row == 0);
    REQUIRE(a1.column == 0);
    REQUIRE(a1.readout == 100.0);
    REQUIRE(h12.well == "H12");
    REQUIRE(h12.row == 7);
    REQUIRE(h12.column == 11);
    REQUIRE(h12.readout == 181.0);
}

TEST_CASE("importText melts a 384-well grid and parses two-letter rows", "[assay]") {
    std::string error;
    const AssayDataset ds = importText(gridBlock(16, 24, 50.0), &error);

    REQUIRE(error.empty());
    REQUIRE(ds.detectedLayout == "384-well grid (16x24)");
    REQUIRE(ds.plates[0].wells.size() == 384);
    const Well& p24 = ds.plates[0].wells.back();
    REQUIRE(p24.well == "P24");
    REQUIRE(p24.row == 15);
    REQUIRE(p24.column == 23);

    int row = 0;
    int column = 0;
    // Bijective base 26: AF is the 32nd row, not the 6th.
    REQUIRE(parseWellName("AF48", &row, &column));
    REQUIRE(row == 31);
    REQUIRE(column == 47);
    REQUIRE(layoutName(standardLayout(32, 48)) == std::string("1536-well grid (32x48)"));
}

TEST_CASE("importText warns instead of clamping a well outside the declared layout", "[assay]") {
    ImportOptions options;
    options.expectedRows = 8;
    options.expectedColumns = 12;
    std::string error;
    const AssayDataset ds =
        importText("well,role,readout\nA1,sample,100\nA13,sample,101\n", options, &error);

    REQUIRE(error.empty());
    REQUIRE(ds.plates[0].wells.size() == 2);
    REQUIRE(ds.plates[0].wells[1].well == "A13");
    // Column 12 (0-based), NOT clamped to the last column of a 96-well plate.
    REQUIRE(ds.plates[0].wells[1].column == 12);
    REQUIRE(ds.warnings.size() == 1);
    REQUIRE(ds.warnings[0].find("A13") != std::string::npos);
    REQUIRE(ds.warnings[0].find("clamped") != std::string::npos);

    int row = 0;
    int column = 0;
    REQUIRE_FALSE(parseWellName("A0", &row, &column));
    const AssayDataset zero = importText("well,readout\nA0,5\n", options, &error);
    REQUIRE(zero.warnings.size() == 1);
    REQUIRE(zero.warnings[0].find("A0") != std::string::npos);
}

TEST_CASE("normalizations are parallel to the wells and refuse missing controls", "[assay]") {
    Plate plate;
    plate.id = "n";
    plate.rows = 1;
    plate.columns = 4;
    const WellRole roles[] = {WellRole::PositiveControl, WellRole::NegativeControl,
                              WellRole::Sample, WellRole::Sample};
    const double values[] = {100.0, 10.0, 55.0, 30.0};
    for (int i = 0; i < 4; ++i) {
        Well w;
        w.role = roles[i];
        w.row = 0;
        w.column = i;
        w.readout = values[i];
        plate.wells.push_back(w);
    }

    const std::vector<double> poc = normalize(plate, Normalization::PercentOfControl);
    REQUIRE(poc.size() == plate.wells.size());
    REQUIRE_THAT(poc[0], WithinAbs(1000.0, 1e-12));   // 100 against a 10 vehicle
    REQUIRE_THAT(poc[1], WithinAbs(100.0, 1e-12));

    const std::vector<double> npi = normalize(plate, Normalization::NormalizedPercentInhibition);
    REQUIRE_THAT(npi[0], WithinAbs(100.0, 1e-12));    // positive control = full inhibition
    REQUIRE_THAT(npi[1], WithinAbs(0.0, 1e-12));      // vehicle = none

    // Without a vehicle there is no percent to be of: NaN, never a plausible zero.
    Plate noVehicle = plate;
    noVehicle.wells[1].role = WellRole::Sample;
    const std::vector<double> broken = normalize(noVehicle, Normalization::PercentOfControl);
    REQUIRE(broken.size() == noVehicle.wells.size());
    REQUIRE(std::isnan(broken[0]));
}

TEST_CASE("median-polish B-score removes a linear row gradient", "[assay]") {
    Plate plate;
    plate.id = "grad";
    plate.rows = 8;
    plate.columns = 12;
    for (int r = 0; r < 8; ++r) {
        for (int c = 0; c < 12; ++c) {
            Well w;
            w.role = WellRole::Sample;
            w.row = r;
            w.column = c;
            // A pure 10-per-row gradient plus a deterministic wobble, so the
            // residual MAD is non-zero and a B-score exists at all.
            w.readout = 100.0 + 10.0 * r + static_cast<double>((r * 12 + c) % 5 - 2) * 0.5;
            plate.wells.push_back(w);
        }
    }
    const std::vector<double> raw = normalize(plate, Normalization::None);
    const std::vector<double> bs = normalize(plate, Normalization::BScore);
    REQUIRE(bs.size() == plate.wells.size());

    // Row effect measured against the within-row noise it has to be seen over.
    // Dividing residuals by their MAD makes any purely relative measure
    // scale-invariant, so the denominator must be the within-row spread.
    auto rowEffectToNoise = [&](const std::vector<double>& v) {
        std::vector<double> means;
        double noise = 0.0;
        for (int r = 0; r < 8; ++r) {
            std::vector<double> row;
            for (std::size_t i = 0; i < plate.wells.size(); ++i) {
                if (plate.wells[i].row == r) row.push_back(v[i]);
            }
            means.push_back(mean(row));
            noise += stdDev(row) / 8.0;
        }
        const auto lohi = std::minmax_element(means.begin(), means.end());
        return (*lohi.second - *lohi.first) / noise;
    };
    const double before = rowEffectToNoise(raw);
    const double after = rowEffectToNoise(bs);
    REQUIRE(before > 50.0);
    REQUIRE(1.0 - after / before > 0.90);
}

TEST_CASE("outlier rules flag without deleting", "[assay]") {
    Plate plate;
    plate.id = "gr";
    plate.rows = 1;
    plate.columns = 10;
    const double values[10] = {100, 101, 99, 100, 100, 101, 99, 100, 100, 130};
    for (int i = 0; i < 10; ++i) {
        Well w;
        w.role = WellRole::PositiveControl;
        w.row = 0;
        w.column = i;
        w.well = "A" + std::to_string(i + 1);
        w.readout = values[i];
        plate.wells.push_back(w);
    }

    SECTION("Grubbs at n=10 flags exactly the planted outlier") {
        REQUIRE(applyGrubbs(plate, WellRole::PositiveControl) == 1);
        REQUIRE(plate.wells.size() == 10);          // nothing removed from the plate
        REQUIRE(plate.wells[9].excluded);
        REQUIRE(plate.wells[9].exclusionRule == "grubbs.two-sided.alpha05");
        std::size_t included = 0;
        for (const Well& w : plate.wells) {
            if (!w.excluded) ++included;
        }
        REQUIRE(included == 9);
    }

    SECTION("Dixon and Tukey flag the same well with their own rule ids") {
        REQUIRE(applyDixonQ(plate, WellRole::PositiveControl) == 1);
        REQUIRE(plate.wells[9].exclusionRule == "dixon.q.alpha05");
        Plate fences = plate;
        for (Well& w : fences.wells) {
            w.excluded = false;
            w.exclusionRule.clear();
        }
        REQUIRE(applyTukeyFences(fences, WellRole::PositiveControl) == 1);
        REQUIRE(fences.wells[9].exclusionRule == "tukey.fence.1.5iqr");
    }

    SECTION("a clean set of controls is left alone") {
        Plate clean = plate;
        clean.wells[9].readout = 100.5;
        REQUIRE(applyGrubbs(clean, WellRole::PositiveControl) == 0);
        REQUIRE(applyDixonQ(clean, WellRole::PositiveControl) == 0);
        for (const Well& w : clean.wells) {
            REQUIRE_FALSE(w.excluded);
        }
    }

    SECTION("an excluded well no longer votes in a normalization") {
        applyGrubbs(plate, WellRole::PositiveControl);
        const std::vector<double> z = normalize(plate, Normalization::ZScore);
        REQUIRE(z.size() == plate.wells.size());
        // No sample wells at all, so the sample-population z-score is undefined
        // rather than silently computed over the controls.
        REQUIRE(std::isnan(z[0]));
    }
}
