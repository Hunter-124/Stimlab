// assay/Dataset.h - plate import, normalization and opt-in outlier rules.
//
// WHY this file is separate from the fitters: everything here is bookkeeping on
// real measurements, and bookkeeping is where experimental data is usually lost.
// Two rules follow from that and are enforced by the code below rather than by
// convention:
//
//   1. Nothing is silently dropped. An input column BioCAD does not understand
//      survives in AssayDataset::metadata; a well name that does not fit the
//      declared layout is kept with its parsed coordinates and a warning, never
//      clamped into a neighbouring well.
//   2. No rule deletes a point. Grubbs, Dixon and the Tukey fences set
//      Well::excluded plus Well::exclusionRule, so the plot can hollow the marker
//      and the reader can see what was removed and by which rule.
#pragma once

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

#include "data/Assay.h"

namespace biocad::assay {

// Layout the caller knows to be true, e.g. from the plate barcode or the reader
// method. Zero means "deduce from the observed well names". Declaring it is what
// lets the importer say "A13 is not a well on a 96-well plate" instead of
// promoting the whole plate to 384 because one label was mistyped.
struct ImportOptions {
    int         expectedRows = 0;
    int         expectedColumns = 0;
    std::string datasetId;
};

// Parse a long CSV/TSV table or a 96/384/1536 grid export. Returns a dataset
// whose `detectedLayout` names the shape that was recognised ("long",
// "96-well grid (8x12)", ...). `error` receives a message and the result is empty
// only when the text cannot be interpreted at all (no header, no well column);
// recoverable problems become AssayDataset::warnings.
AssayDataset importText(std::string_view text, std::string* error);
AssayDataset importText(std::string_view text, const ImportOptions& options, std::string* error);

// Bijective base-26 row letters ("A" -> 0, "Z" -> 25, "AA" -> 26, "AF" -> 31) and
// a column index, both 0-based on output ("P24" -> row 15, column 23). Returns
// false when the label is not <letters><digits> or the column number is zero.
bool parseWellName(std::string_view name, int* row, int* column);

// The smallest standard microplate layout that contains (rows, columns), or
// {0,0} when the extents exceed 1536.
struct PlateLayout { int rows = 0; int columns = 0; };
PlateLayout standardLayout(int rowsNeeded, int columnsNeeded);
const char* layoutName(const PlateLayout& layout);   // "96-well grid (8x12)" etc.

// ---------------------------------------------------------------------------
// Normalization. Every mode returns one value per Plate::wells entry, in order.
// A well whose normalization has no defined value (a missing control, a zero
// denominator) yields NaN: the caller renders that as notComputed rather than as
// a zero, because a zero here looks like a real measurement.
// ---------------------------------------------------------------------------
std::vector<double> normalize(const Plate& plate, Normalization mode);

// ---------------------------------------------------------------------------
// Outlier rules. All opt-in, all non-destructive; each returns the number of
// wells it newly excluded. `scope` restricts the test to one role (the usual
// case: test replicate controls, not the whole plate).
// ---------------------------------------------------------------------------

// Two-sided Grubbs at alpha = 0.05, SINGLE pass: the classic test detects one
// outlier, and iterating it inflates the type-I error rate well past 5%.
// Critical values are TABULATED for n = 3..30 (the standard two-sided 5% Grubbs
// table) rather than computed from a t quantile, so the numbers are exactly the
// published ones; n outside that range returns 0 with no exclusions.
std::size_t applyGrubbs(Plate& plate, WellRole scope);

// Dixon's Q at alpha = 0.05, tabulated for n = 3..10 (Dixon's original two-sided
// 95% table); outside that range the test does not apply and nothing is flagged.
std::size_t applyDixonQ(Plate& plate, WellRole scope);

// Tukey fences: below Q1 - 1.5*IQR or above Q3 + 1.5*IQR, quartiles by the
// linear-interpolation definition (R type 7).
std::size_t applyTukeyFences(Plate& plate, WellRole scope);

}  // namespace biocad::assay
