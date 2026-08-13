#include "assay/Dataset.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <map>
#include <sstream>
#include <unordered_map>

#include "assay/Stats.h"

namespace biocad::assay {
namespace {

std::string trim(std::string_view s) {
    std::size_t b = 0;
    std::size_t e = s.size();
    while (b < e && (s[b] == ' ' || s[b] == '\t' || s[b] == '\r' || s[b] == '"')) ++b;
    while (e > b && (s[e - 1] == ' ' || s[e - 1] == '\t' || s[e - 1] == '\r' || s[e - 1] == '"')) --e;
    return std::string(s.substr(b, e - b));
}

// Column names are matched case-insensitively and order-free, and '-'/' ' are
// folded to '_' so "Conc Unit", "conc-unit" and "CONC_UNIT" are one column.
std::string canonicalKey(std::string_view s) {
    std::string out;
    out.reserve(s.size());
    for (char c : trim(s)) {
        if (c == ' ' || c == '-') c = '_';
        out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    return out;
}

std::vector<std::string> splitLine(const std::string& line, char delim) {
    std::vector<std::string> out;
    std::string field;
    bool quoted = false;
    for (char c : line) {
        if (c == '"') {
            quoted = !quoted;
            continue;
        }
        if (c == delim && !quoted) {
            out.push_back(trim(field));
            field.clear();
            continue;
        }
        field.push_back(c);
    }
    out.push_back(trim(field));
    return out;
}

std::vector<std::string> splitLines(std::string_view text) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : text) {
        if (c == '\n') {
            out.push_back(cur);
            cur.clear();
        } else if (c != '\r') {
            cur.push_back(c);
        }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

bool blank(const std::string& line) { return trim(line).empty(); }

char detectDelimiter(const std::vector<std::string>& lines) {
    for (const std::string& l : lines) {
        if (blank(l)) continue;
        if (l.find('\t') != std::string::npos) return '\t';
        if (l.find(',') != std::string::npos) return ',';
        if (l.find(';') != std::string::npos) return ';';
    }
    return ',';
}

bool parseNumber(const std::string& s, double* out) {
    const std::string t = trim(s);
    if (t.empty()) return false;
    const char* begin = t.c_str();
    char* end = nullptr;
    const double v = std::strtod(begin, &end);
    if (end == begin) return false;
    while (*end == ' ') ++end;
    if (*end != '\0') return false;
    *out = v;
    return true;
}

bool parseInteger(const std::string& s, int* out) {
    double v = 0.0;
    if (!parseNumber(s, &v)) return false;
    *out = static_cast<int>(std::lround(v));
    return true;
}

WellRole parseRole(const std::string& raw, bool* known) {
    const std::string k = canonicalKey(raw);
    *known = true;
    if (k.empty() || k == "unknown" || k == "na") return WellRole::Unknown;
    if (k == "sample" || k == "compound" || k == "test") return WellRole::Sample;
    if (k == "positive" || k == "pos" || k == "positive_control" || k == "max" || k == "high")
        return WellRole::PositiveControl;
    if (k == "negative" || k == "neg" || k == "negative_control" || k == "vehicle" || k == "min" ||
        k == "low" || k == "dmso")
        return WellRole::NegativeControl;
    if (k == "blank" || k == "buffer") return WellRole::Blank;
    if (k == "reference" || k == "ref") return WellRole::Reference;
    if (k == "empty") return WellRole::Empty;
    *known = false;
    return WellRole::Unknown;
}

// A grid header is the row of column numbers a reader export puts above the
// block: an optional leading label cell, then 1..12, 1..24 or 1..48 in order.
bool isGridHeader(const std::vector<std::string>& fields, int* columns) {
    std::size_t start = 0;
    if (!fields.empty() && trim(fields[0]).empty()) start = 1;
    // Trailing empty cells are common in exports padded to a fixed width.
    std::size_t end = fields.size();
    while (end > start && trim(fields[end - 1]).empty()) --end;
    const std::size_t count = end - start;
    if (count != 12 && count != 24 && count != 48) return false;
    for (std::size_t i = start; i < end; ++i) {
        int v = 0;
        if (!parseInteger(fields[i], &v)) return false;
        if (v != static_cast<int>(i - start) + 1) return false;
    }
    *columns = static_cast<int>(count);
    return true;
}

const std::vector<std::string>& longColumnKeys() {
    static const std::vector<std::string> keys = {
        "plate_id", "well", "role", "sample_id", "series_id", "concentration", "conc_unit",
        "replicate", "readout", "readout_unit", "time_s", "temperature_c", "excluded",
        "exclusion_rule"};
    return keys;
}

bool truthy(const std::string& s) {
    const std::string k = canonicalKey(s);
    return k == "1" || k == "true" || k == "yes" || k == "y" || k == "t";
}

void finishPlates(AssayDataset& ds, const ImportOptions& options) {
    for (Plate& p : ds.plates) {
        int maxRow = -1;
        int maxCol = -1;
        for (const Well& w : p.wells) {
            maxRow = std::max(maxRow, w.row);
            maxCol = std::max(maxCol, w.column);
            if (p.readoutUnit.empty()) p.readoutUnit = w.readoutUnit;
        }
        // Deduced from the observed labels, not guessed from the well count: a
        // partially filled 384-well plate is not a 96-well plate.
        p.rows = maxRow + 1;
        p.columns = maxCol + 1;
        if (options.expectedRows > 0 && options.expectedColumns > 0) {
            for (const Well& w : p.wells) {
                if (w.row >= options.expectedRows || w.column >= options.expectedColumns) {
                    ds.warnings.push_back(
                        "well " + w.well + " on plate " + p.id + " lies outside the declared " +
                        std::to_string(options.expectedRows) + "x" +
                        std::to_string(options.expectedColumns) +
                        " layout; kept with its parsed coordinates rather than clamped");
                }
            }
        }
    }
}

Plate& plateFor(AssayDataset& ds, const std::string& id) {
    for (Plate& p : ds.plates) {
        if (p.id == id) return p;
    }
    ds.plates.push_back(Plate{});
    ds.plates.back().id = id;
    return ds.plates.back();
}

// --- collection helpers shared by normalize() and the outlier rules ----------

std::vector<double> roleValues(const Plate& plate, WellRole role) {
    std::vector<double> out;
    for (const Well& w : plate.wells) {
        if (w.role == role && !w.excluded) out.push_back(w.readout);
    }
    return out;
}

// Sample and Unknown wells are the population a plate-level normalization is
// computed over: controls define the scale, they are not part of it. Unknown is
// included because a bare grid export carries no roles at all.
bool isPopulation(const Well& w) {
    return !w.excluded && (w.role == WellRole::Sample || w.role == WellRole::Unknown);
}

std::vector<double> populationValues(const Plate& plate) {
    std::vector<double> out;
    for (const Well& w : plate.wells) {
        if (isPopulation(w)) out.push_back(w.readout);
    }
    return out;
}

std::vector<double> bScore(const Plate& plate) {
    const double nan = std::nan("");
    std::vector<double> out(plate.wells.size(), nan);
    if (plate.rows <= 0 || plate.columns <= 0) return out;

    const std::size_t rows = static_cast<std::size_t>(plate.rows);
    const std::size_t cols = static_cast<std::size_t>(plate.columns);
    std::vector<double> resid(rows * cols, nan);
    for (const Well& w : plate.wells) {
        if (!isPopulation(w)) continue;
        if (w.row < 0 || w.column < 0 || w.row >= plate.rows || w.column >= plate.columns) continue;
        resid[static_cast<std::size_t>(w.row) * cols + static_cast<std::size_t>(w.column)] =
            w.readout;
    }

    // Tukey median polish: alternate removing row and column medians until the
    // largest removed median stops moving. Missing cells simply do not vote,
    // which is why a partial plate still polishes.
    std::vector<double> line;
    for (int iter = 0; iter < 200; ++iter) {
        double largest = 0.0;
        for (std::size_t r = 0; r < rows; ++r) {
            line.clear();
            for (std::size_t c = 0; c < cols; ++c) {
                const double v = resid[r * cols + c];
                if (!std::isnan(v)) line.push_back(v);
            }
            if (line.empty()) continue;
            const double m = median(line);
            for (std::size_t c = 0; c < cols; ++c) {
                if (!std::isnan(resid[r * cols + c])) resid[r * cols + c] -= m;
            }
            largest = std::max(largest, std::abs(m));
        }
        for (std::size_t c = 0; c < cols; ++c) {
            line.clear();
            for (std::size_t r = 0; r < rows; ++r) {
                const double v = resid[r * cols + c];
                if (!std::isnan(v)) line.push_back(v);
            }
            if (line.empty()) continue;
            const double m = median(line);
            for (std::size_t r = 0; r < rows; ++r) {
                if (!std::isnan(resid[r * cols + c])) resid[r * cols + c] -= m;
            }
            largest = std::max(largest, std::abs(m));
        }
        if (largest < 1e-12) break;
    }

    std::vector<double> flat;
    for (double v : resid) {
        if (!std::isnan(v)) flat.push_back(v);
    }
    const double scale = kMadToSigma * medianAbsoluteDeviation(flat);
    if (!(scale > 0.0)) return out;   // a perfectly flat plate has no B-score

    for (std::size_t i = 0; i < plate.wells.size(); ++i) {
        const Well& w = plate.wells[i];
        if (!isPopulation(w)) continue;
        if (w.row < 0 || w.column < 0 || w.row >= plate.rows || w.column >= plate.columns) continue;
        out[i] = resid[static_cast<std::size_t>(w.row) * cols +
                       static_cast<std::size_t>(w.column)] / scale;
    }
    return out;
}

// Indices into plate.wells that a rule may test: one role, still included.
std::vector<std::size_t> scopeIndices(const Plate& plate, WellRole scope) {
    std::vector<std::size_t> idx;
    for (std::size_t i = 0; i < plate.wells.size(); ++i) {
        if (plate.wells[i].role == scope && !plate.wells[i].excluded) idx.push_back(i);
    }
    return idx;
}

void exclude(Plate& plate, std::size_t index, const char* rule) {
    plate.wells[index].excluded = true;
    plate.wells[index].exclusionRule = rule;
}

}  // namespace

// ---------------------------------------------------------------------------
// Well names and layouts
// ---------------------------------------------------------------------------

bool parseWellName(std::string_view name, int* row, int* column) {
    const std::string t = trim(name);
    std::size_t i = 0;
    int r = 0;
    while (i < t.size() && std::isalpha(static_cast<unsigned char>(t[i]))) {
        // Bijective base 26: A..Z, AA..AZ, ... so "AF" is 31, not 5.
        r = r * 26 + (std::toupper(static_cast<unsigned char>(t[i])) - 'A' + 1);
        ++i;
    }
    if (i == 0 || i == t.size()) return false;
    int c = 0;
    for (; i < t.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(t[i]))) return false;
        c = c * 10 + (t[i] - '0');
    }
    if (c <= 0) return false;
    *row = r - 1;
    *column = c - 1;
    return true;
}

PlateLayout standardLayout(int rowsNeeded, int columnsNeeded) {
    static const PlateLayout kLayouts[] = {{8, 12}, {16, 24}, {32, 48}};
    for (const PlateLayout& l : kLayouts) {
        if (rowsNeeded <= l.rows && columnsNeeded <= l.columns) return l;
    }
    return PlateLayout{0, 0};
}

const char* layoutName(const PlateLayout& layout) {
    if (layout.rows == 8 && layout.columns == 12) return "96-well grid (8x12)";
    if (layout.rows == 16 && layout.columns == 24) return "384-well grid (16x24)";
    if (layout.rows == 32 && layout.columns == 48) return "1536-well grid (32x48)";
    return "non-standard grid";
}

// ---------------------------------------------------------------------------
// Import
// ---------------------------------------------------------------------------

AssayDataset importText(std::string_view text, std::string* error) {
    return importText(text, ImportOptions{}, error);
}

AssayDataset importText(std::string_view text, const ImportOptions& options, std::string* error) {
    AssayDataset ds;
    ds.id = options.datasetId.empty() ? "assay" : options.datasetId;
    if (error) error->clear();

    const std::vector<std::string> lines = splitLines(text);
    const char delim = detectDelimiter(lines);

    // Grid path first: a grid header is unambiguous, and a long table can never
    // produce one because its header row is column names, not 1..N.
    for (std::size_t i = 0; i < lines.size(); ++i) {
        if (blank(lines[i])) continue;
        int gridColumns = 0;
        if (isGridHeader(splitLine(lines[i], delim), &gridColumns)) {
            std::string pendingLabel;
            int plateCounter = 0;
            int observedRows = 0;
            int observedCols = 0;
            for (std::size_t j = 0; j < lines.size(); ++j) {
                if (blank(lines[j])) continue;
                std::vector<std::string> fields = splitLine(lines[j], delim);
                int cols = 0;
                if (isGridHeader(fields, &cols)) {
                    ++plateCounter;
                    std::string id = pendingLabel.empty()
                                         ? "plate" + std::to_string(plateCounter)
                                         : pendingLabel;
                    pendingLabel.clear();
                    Plate& plate = plateFor(ds, id);
                    // Read the block that follows: rows keep coming while the first
                    // cell is a well row label.
                    for (std::size_t k = j + 1; k < lines.size(); ++k) {
                        if (blank(lines[k])) break;
                        std::vector<std::string> row = splitLine(lines[k], delim);
                        int probeCols = 0;
                        if (isGridHeader(row, &probeCols)) break;
                        const std::string label = trim(row.empty() ? std::string() : row[0]);
                        int r = 0;
                        int ignored = 0;
                        if (label.empty() || !parseWellName(label + "1", &r, &ignored)) break;
                        for (int c = 0; c < cols && static_cast<std::size_t>(c) + 1 < row.size();
                             ++c) {
                            double v = 0.0;
                            if (!parseNumber(row[static_cast<std::size_t>(c) + 1], &v)) {
                                if (!trim(row[static_cast<std::size_t>(c) + 1]).empty()) {
                                    ds.warnings.push_back("non-numeric cell at " + label +
                                                          std::to_string(c + 1) + " on plate " +
                                                          id + "; well omitted");
                                }
                                continue;
                            }
                            Well w;
                            w.plateId = id;
                            w.well = label + std::to_string(c + 1);
                            w.row = r;
                            w.column = c;
                            w.readout = v;
                            plate.wells.push_back(w);
                            observedRows = std::max(observedRows, r + 1);
                            observedCols = std::max(observedCols, c + 1);
                        }
                        j = k;
                    }
                    continue;
                }
                if (fields.size() <= 2) pendingLabel = trim(fields[0]);
            }
            const PlateLayout layout = standardLayout(observedRows, std::max(observedCols, gridColumns));
            ds.detectedLayout = layoutName(layout);
            finishPlates(ds, options);
            if (ds.plates.empty() && error) *error = "grid header found but no data rows followed";
            return ds;
        }
    }

    // Long path.
    std::size_t headerLine = lines.size();
    for (std::size_t i = 0; i < lines.size(); ++i) {
        if (!blank(lines[i])) {
            headerLine = i;
            break;
        }
    }
    if (headerLine == lines.size()) {
        if (error) *error = "input is empty";
        return ds;
    }

    const std::vector<std::string> header = splitLine(lines[headerLine], delim);
    std::unordered_map<std::string, std::size_t> known;
    std::vector<std::pair<std::string, std::size_t>> unknown;
    for (std::size_t i = 0; i < header.size(); ++i) {
        const std::string key = canonicalKey(header[i]);
        if (key.empty()) continue;
        const std::vector<std::string>& keys = longColumnKeys();
        if (std::find(keys.begin(), keys.end(), key) != keys.end()) {
            known[key] = i;
        } else {
            unknown.emplace_back(trim(header[i]), i);
        }
    }
    if (known.find("well") == known.end()) {
        if (error) *error = "no recognised layout: long input needs a 'well' column";
        return ds;
    }

    auto field = [&](const std::vector<std::string>& row, const char* key) -> std::string {
        auto it = known.find(key);
        if (it == known.end() || it->second >= row.size()) return std::string();
        return row[it->second];
    };

    // Unknown columns are accumulated per row and collapsed at the end: a single
    // distinct value is stored as-is (the run-level constants an instrument
    // exports), otherwise every row value is joined with '|' in the same order as
    // the wells, so the association is recoverable. Nothing is dropped.
    std::vector<std::vector<std::string>> extras(unknown.size());

    for (std::size_t i = headerLine + 1; i < lines.size(); ++i) {
        if (blank(lines[i])) continue;
        const std::vector<std::string> row = splitLine(lines[i], delim);
        const std::string wellName = field(row, "well");
        int r = 0;
        int c = 0;
        if (!parseWellName(wellName, &r, &c)) {
            ds.warnings.push_back("unparseable well label '" + wellName + "' on line " +
                                  std::to_string(i + 1) + "; row skipped");
            continue;
        }
        Well w;
        w.plateId = field(row, "plate_id");
        if (w.plateId.empty()) w.plateId = "plate1";
        w.well = trim(wellName);
        w.row = r;
        w.column = c;
        bool roleKnown = true;
        w.role = parseRole(field(row, "role"), &roleKnown);
        if (!roleKnown) {
            ds.warnings.push_back("unrecognised role '" + trim(field(row, "role")) + "' at " +
                                  w.well + "; recorded as unknown");
        }
        w.sampleId = field(row, "sample_id");
        w.seriesId = field(row, "series_id");
        parseNumber(field(row, "concentration"), &w.concentration);
        w.concUnit = field(row, "conc_unit");
        parseInteger(field(row, "replicate"), &w.replicate);
        if (!parseNumber(field(row, "readout"), &w.readout)) {
            ds.warnings.push_back("missing or non-numeric readout at " + w.well +
                                  "; well kept with readout 0 and excluded");
            w.excluded = true;
            w.exclusionRule = "import.no-readout";
        }
        w.readoutUnit = field(row, "readout_unit");
        parseNumber(field(row, "time_s"), &w.timeS);
        parseNumber(field(row, "temperature_c"), &w.temperatureC);
        if (truthy(field(row, "excluded"))) {
            w.excluded = true;
            w.exclusionRule = field(row, "exclusion_rule");
            if (w.exclusionRule.empty()) w.exclusionRule = "import.excluded";
        }
        plateFor(ds, w.plateId).wells.push_back(w);

        for (std::size_t u = 0; u < unknown.size(); ++u) {
            const std::size_t col = unknown[u].second;
            extras[u].push_back(col < row.size() ? row[col] : std::string());
        }
    }

    for (std::size_t u = 0; u < unknown.size(); ++u) {
        const std::vector<std::string>& vals = extras[u];
        if (vals.empty()) continue;
        const bool constant = std::all_of(vals.begin(), vals.end(),
                                          [&](const std::string& v) { return v == vals.front(); });
        std::string joined = vals.front();
        if (!constant) {
            joined.clear();
            for (std::size_t k = 0; k < vals.size(); ++k) {
                if (k) joined.push_back('|');
                joined += vals[k];
            }
        }
        ds.metadata.emplace_back(unknown[u].first, joined);
    }

    ds.detectedLayout = "long";
    finishPlates(ds, options);
    if (ds.plates.empty() && error) *error = "header parsed but no data rows";
    return ds;
}

// ---------------------------------------------------------------------------
// Normalization
// ---------------------------------------------------------------------------

std::vector<double> normalize(const Plate& plate, Normalization mode) {
    const double nan = std::nan("");
    std::vector<double> out(plate.wells.size(), nan);
    if (mode == Normalization::BScore) return bScore(plate);

    const std::vector<double> pos = roleValues(plate, WellRole::PositiveControl);
    const std::vector<double> neg = roleValues(plate, WellRole::NegativeControl);
    const std::vector<double> pop = populationValues(plate);

    double a = nan;   // offset
    double b = nan;   // divisor
    switch (mode) {
        case Normalization::None:
            for (std::size_t i = 0; i < plate.wells.size(); ++i) out[i] = plate.wells[i].readout;
            return out;
        case Normalization::PercentOfControl:
            // 100 * x / mean(vehicle). The vehicle (negative) control is the 100%
            // reference; without it there is no percent to be of.
            a = 0.0;
            b = neg.empty() ? nan : mean(neg) / 100.0;
            break;
        case Normalization::NormalizedPercentInhibition:
            // 100 * (vehicle - x) / (vehicle - max inhibition), so 0% is the
            // vehicle and 100% is the positive (fully inhibited) control.
            if (pos.empty() || neg.empty()) return out;
            a = mean(neg);
            b = (mean(neg) - mean(pos)) / 100.0;
            for (std::size_t i = 0; i < plate.wells.size(); ++i) {
                if (b != 0.0) out[i] = (a - plate.wells[i].readout) / b;
            }
            return out;
        case Normalization::ZScore:
            a = mean(pop);
            b = stdDev(pop);
            break;
        case Normalization::RobustZScore:
            a = median(pop);
            b = kMadToSigma * medianAbsoluteDeviation(pop);
            break;
        case Normalization::BScore:
            break;   // handled above
    }
    if (std::isnan(a) || std::isnan(b) || b == 0.0) return out;
    for (std::size_t i = 0; i < plate.wells.size(); ++i) {
        out[i] = (plate.wells[i].readout - a) / b;
    }
    return out;
}

// ---------------------------------------------------------------------------
// Outlier rules
// ---------------------------------------------------------------------------

std::size_t applyGrubbs(Plate& plate, WellRole scope) {
    // Two-sided Grubbs critical values at alpha = 0.05, n = 3..30, as published
    // (Grubbs 1969 / ISO 5725). Tabulated rather than derived so the thresholds
    // are exactly the reference numbers a reviewer can look up.
    static const double kG[] = {1.1543, 1.4812, 1.7150, 1.8871, 2.0200, 2.1266, 2.2150, 2.2900,
                                2.3547, 2.4116, 2.4620, 2.5073, 2.5483, 2.5857, 2.6200, 2.6516,
                                2.6809, 2.7082, 2.7337, 2.7576, 2.7801, 2.8013, 2.8214, 2.8404,
                                2.8586, 2.8759, 2.8924, 2.9082};
    const std::vector<std::size_t> idx = scopeIndices(plate, scope);
    const std::size_t n = idx.size();
    if (n < 3 || n > 30) return 0;

    std::vector<double> vals;
    vals.reserve(n);
    for (std::size_t i : idx) vals.push_back(plate.wells[i].readout);
    const double m = mean(vals);
    const double sd = stdDev(vals);
    if (!(sd > 0.0)) return 0;

    std::size_t worst = 0;
    double worstG = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        const double g = std::abs(vals[k] - m) / sd;
        if (g > worstG) {
            worstG = g;
            worst = k;
        }
    }
    if (worstG <= kG[n - 3]) return 0;
    exclude(plate, idx[worst], "grubbs.two-sided.alpha05");
    return 1;
}

std::size_t applyDixonQ(Plate& plate, WellRole scope) {
    // Dixon's Q at 95% confidence, n = 3..10, tabulated (Dixon 1951 / Rorabacher
    // 1991 two-sided 95% column). Beyond n = 10 the test is superseded by Grubbs.
    static const double kQ[] = {0.970, 0.829, 0.710, 0.625, 0.568, 0.526, 0.493, 0.466};
    const std::vector<std::size_t> idx = scopeIndices(plate, scope);
    const std::size_t n = idx.size();
    if (n < 3 || n > 10) return 0;

    std::vector<std::size_t> order = idx;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return plate.wells[a].readout < plate.wells[b].readout;
    });
    const double lo = plate.wells[order.front()].readout;
    const double hi = plate.wells[order.back()].readout;
    const double range = hi - lo;
    if (!(range > 0.0)) return 0;
    const double qLow = (plate.wells[order[1]].readout - lo) / range;
    const double qHigh = (hi - plate.wells[order[n - 2]].readout) / range;
    const double crit = kQ[n - 3];
    std::size_t flagged = 0;
    if (qLow > crit) {
        exclude(plate, order.front(), "dixon.q.alpha05");
        ++flagged;
    }
    if (qHigh > crit) {
        exclude(plate, order.back(), "dixon.q.alpha05");
        ++flagged;
    }
    return flagged;
}

std::size_t applyTukeyFences(Plate& plate, WellRole scope) {
    const std::vector<std::size_t> idx = scopeIndices(plate, scope);
    if (idx.size() < 4) return 0;
    std::vector<double> vals;
    vals.reserve(idx.size());
    for (std::size_t i : idx) vals.push_back(plate.wells[i].readout);
    const double q1 = quantileType7(vals, 0.25);
    const double q3 = quantileType7(vals, 0.75);
    const double iqr = q3 - q1;
    if (!(iqr > 0.0)) return 0;
    const double loFence = q1 - 1.5 * iqr;
    const double hiFence = q3 + 1.5 * iqr;
    std::size_t flagged = 0;
    for (std::size_t i : idx) {
        const double v = plate.wells[i].readout;
        if (v < loFence || v > hiFence) {
            exclude(plate, i, "tukey.fence.1.5iqr");
            ++flagged;
        }
    }
    return flagged;
}

}  // namespace biocad::assay
