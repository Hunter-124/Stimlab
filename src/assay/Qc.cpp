#include "assay/Qc.h"

#include <algorithm>
#include <cmath>
#include <string>

#include "assay/Stats.h"

namespace biocad::assay {
namespace {

constexpr const char* kZPrimeSource = "Zhang, Chung & Oldenburg 1999, J Biomol Screen 4:67";
constexpr const char* kSsmdSource = "Zhang 2007 strictly standardized mean difference";
constexpr const char* kPlateSource = "plate readout statistics";

std::vector<double> included(const Plate& plate, WellRole role) {
    std::vector<double> out;
    for (const Well& w : plate.wells) {
        if (w.role == role && !w.excluded) out.push_back(w.readout);
    }
    return out;
}

bool isPopulation(const Well& w) {
    return !w.excluded && (w.role == WellRole::Sample || w.role == WellRole::Unknown);
}

// A coefficient of variation is only meaningful on a ratio scale with a true zero
// and strictly positive values; on a signed readout (a delta, a normalized score)
// the mean can pass through zero and %CV explodes. That is a data question, not a
// formatting question, so it returns notComputed.
Quantity cvPercent(const std::vector<double>& v) {
    if (v.size() < 2) return notComputed("at least two included wells");
    for (double x : v) {
        if (!(x > 0.0)) return notComputed("%CV requires positive ratio-scale data");
    }
    const double m = mean(v);
    return makeQuantity(100.0 * stdDev(v) / m, "%", 0.0, Provenance::Measured, kPlateSource);
}

// Regularized lower incomplete gamma by series, and the upper one by the modified
// Lentz continued fraction. Numerical Recipes' pair; both converge to ~1e-15 here.
double gammaSeries(double a, double x) {
    const double lg = std::lgamma(a);
    double ap = a;
    double sum = 1.0 / a;
    double del = sum;
    for (int n = 0; n < 500; ++n) {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if (std::abs(del) < std::abs(sum) * 1e-16) break;
    }
    return sum * std::exp(-x + a * std::log(x) - lg);
}

double gammaContinuedFraction(double a, double x) {
    const double lg = std::lgamma(a);
    const double tiny = 1e-300;
    double b = x + 1.0 - a;
    double c = 1.0 / tiny;
    double d = 1.0 / b;
    double h = d;
    for (int i = 1; i < 500; ++i) {
        const double an = -static_cast<double>(i) * (static_cast<double>(i) - a);
        b += 2.0;
        d = an * d + b;
        if (std::abs(d) < tiny) d = tiny;
        c = b + an / c;
        if (std::abs(c) < tiny) c = tiny;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (std::abs(del - 1.0) < 1e-16) break;
    }
    return std::exp(-x + a * std::log(x) - lg) * h;
}

// Ranks with ties averaged; also returns sum(t^3 - t) over tie groups, which is
// what both the U variance and the H statistic need for their tie correction.
std::vector<double> midRanks(const std::vector<double>& v, double* tieTerm) {
    const std::size_t n = v.size();
    std::vector<std::size_t> order(n);
    for (std::size_t i = 0; i < n; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return v[a] < v[b]; });
    std::vector<double> rank(n, 0.0);
    double ties = 0.0;
    std::size_t i = 0;
    while (i < n) {
        std::size_t j = i;
        while (j + 1 < n && v[order[j + 1]] == v[order[i]]) ++j;
        const double r = 0.5 * (static_cast<double>(i) + static_cast<double>(j)) + 1.0;
        for (std::size_t k = i; k <= j; ++k) rank[order[k]] = r;
        const double t = static_cast<double>(j - i + 1);
        ties += t * t * t - t;
        i = j + 1;
    }
    if (tieTerm) *tieTerm = ties;
    return rank;
}

}  // namespace

double chiSquareSurvival(double x, double degreesOfFreedom) {
    if (!(degreesOfFreedom > 0.0)) return std::nan("");
    if (x <= 0.0) return 1.0;
    const double a = 0.5 * degreesOfFreedom;
    const double z = 0.5 * x;
    if (z < a + 1.0) return 1.0 - gammaSeries(a, z);
    return gammaContinuedFraction(a, z);
}

double mannWhitneyP(const std::vector<double>& a, const std::vector<double>& b) {
    const double n1 = static_cast<double>(a.size());
    const double n2 = static_cast<double>(b.size());
    if (a.empty() || b.empty()) return std::nan("");
    std::vector<double> all = a;
    all.insert(all.end(), b.begin(), b.end());
    double ties = 0.0;
    const std::vector<double> rank = midRanks(all, &ties);
    double r1 = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) r1 += rank[i];

    const double u = n1 * n2 + n1 * (n1 + 1.0) / 2.0 - r1;
    const double mu = 0.5 * n1 * n2;
    const double n = n1 + n2;
    const double var = (n1 * n2 / 12.0) * ((n + 1.0) - ties / (n * (n - 1.0)));
    if (!(var > 0.0)) return 1.0;
    const double z = (std::abs(u - mu) - 0.5) / std::sqrt(var);
    if (z <= 0.0) return 1.0;
    return std::erfc(z / std::sqrt(2.0));   // two-sided
}

double kruskalWallisP(const std::vector<std::vector<double>>& groups) {
    std::vector<double> all;
    std::size_t k = 0;
    for (const std::vector<double>& g : groups) {
        if (g.size() < 2) continue;
        ++k;
        all.insert(all.end(), g.begin(), g.end());
    }
    if (k < 2 || all.size() < 3) return std::nan("");
    double ties = 0.0;
    const std::vector<double> rank = midRanks(all, &ties);
    const double n = static_cast<double>(all.size());
    double h = 0.0;
    std::size_t offset = 0;
    for (const std::vector<double>& g : groups) {
        if (g.size() < 2) continue;
        double sum = 0.0;
        for (std::size_t i = 0; i < g.size(); ++i) sum += rank[offset + i];
        h += sum * sum / static_cast<double>(g.size());
        offset += g.size();
    }
    h = 12.0 / (n * (n + 1.0)) * h - 3.0 * (n + 1.0);
    const double correction = 1.0 - ties / (n * n * n - n);
    if (correction > 0.0) h /= correction;
    return chiSquareSurvival(h, static_cast<double>(k) - 1.0);
}

RoleStats roleStatistics(const Plate& plate, WellRole role) {
    RoleStats s;
    s.role = role;
    const std::vector<double> v = included(plate, role);
    s.n = v.size();
    if (v.empty()) {
        const Quantity none = notComputed("at least one included well of this role");
        s.mean = s.sd = s.sem = s.cvPct = s.median = s.mad = s.robustSd = none;
        return s;
    }
    const std::string unit = plate.readoutUnit;
    s.mean = makeQuantity(mean(v), unit, 0.0, Provenance::Measured, kPlateSource);
    s.median = makeQuantity(median(v), unit, 0.0, Provenance::Measured, kPlateSource);
    s.mad = makeQuantity(medianAbsoluteDeviation(v), unit, 0.0, Provenance::Measured, kPlateSource);
    s.robustSd = makeQuantity(kMadToSigma * medianAbsoluteDeviation(v), unit, 0.0,
                              Provenance::Measured, "1.4826 * MAD");
    if (v.size() < 2) {
        const Quantity none = notComputed("at least two included wells of this role");
        s.sd = s.sem = none;
    } else {
        s.sd = makeQuantity(stdDev(v), unit, 0.0, Provenance::Measured, kPlateSource);
        s.sem = makeQuantity(stdError(v), unit, 0.0, Provenance::Measured, kPlateSource);
    }
    s.cvPct = cvPercent(v);
    return s;
}

Quantity blankSubtractedVariance(const Plate& plate) {
    const std::vector<double> blanks = included(plate, WellRole::Blank);
    if (blanks.size() < 2) return notComputed("at least two blank wells");
    std::vector<double> samples;
    for (const Well& w : plate.wells) {
        if (isPopulation(w)) samples.push_back(w.readout);
    }
    if (samples.size() < 2) return notComputed("at least two sample wells");
    const double sx = stdDev(samples);
    const double sb = stdDev(blanks);
    const double var = sx * sx + sb * sb / static_cast<double>(blanks.size());
    const std::string unit = plate.readoutUnit.empty() ? std::string()
                                                       : plate.readoutUnit + "^2";
    return makeQuantity(var, unit, 0.0, Provenance::Measured,
                        "Var(x - mean(blank)) = s_x^2 + s_b^2/n_b");
}

QcReport plateQc(const Plate& plate) {
    QcReport r;
    r.plateId = plate.id;

    const std::vector<double> pos = included(plate, WellRole::PositiveControl);
    const std::vector<double> neg = included(plate, WellRole::NegativeControl);
    const RoleStats ps = roleStatistics(plate, WellRole::PositiveControl);
    const RoleStats ns = roleStatistics(plate, WellRole::NegativeControl);
    r.positiveMean = ps.mean;
    r.positiveSd = ps.sd;
    r.negativeMean = ns.mean;
    r.negativeSd = ns.sd;
    r.cvPositivePct = ps.cvPct;
    r.cvNegativePct = ns.cvPct;

    const bool bothControls = pos.size() >= 2 && neg.size() >= 2;
    if (!bothControls) {
        const char* missing = pos.size() < 2 ? "at least two positive control wells"
                                            : "at least two negative control wells";
        r.zPrime = notComputed(missing);
        r.robustZPrime = notComputed(missing);
        r.ssmd = notComputed(missing);
        r.signalToNoise = notComputed(missing);
        r.warnings.push_back(std::string("Z-prime, robust Z-prime, SSMD and S/N need ") + missing +
                             "; they are not estimated from the extreme wells");
    } else {
        const double mp = mean(pos);
        const double mn = mean(neg);
        const double sp = stdDev(pos);
        const double sn = stdDev(neg);
        const double sep = std::abs(mp - mn);
        if (sep > 0.0) {
            r.zPrime = makeQuantity(1.0 - 3.0 * (sp + sn) / sep, "", 0.0, Provenance::Measured,
                                    kZPrimeSource);
        } else {
            r.zPrime = notComputed("separated control means (|mean_p - mean_n| > 0)");
        }
        const double rp = kMadToSigma * medianAbsoluteDeviation(pos);
        const double rn = kMadToSigma * medianAbsoluteDeviation(neg);
        const double rsep = std::abs(median(pos) - median(neg));
        if (rsep > 0.0) {
            r.robustZPrime = makeQuantity(1.0 - 3.0 * (rp + rn) / rsep, "", 0.0,
                                          Provenance::Measured,
                                          "robust Z-prime from median and 1.4826*MAD");
        } else {
            r.robustZPrime = notComputed("separated control medians");
        }
        const double pooled = std::sqrt(sp * sp + sn * sn);
        if (pooled > 0.0) {
            r.ssmd = makeQuantity((mp - mn) / pooled, "", 0.0, Provenance::Measured, kSsmdSource);
        } else {
            r.ssmd = notComputed("non-zero control variance");
        }
        if (sn > 0.0) {
            r.signalToNoise = makeQuantity((mp - mn) / sn, "", 0.0, Provenance::Measured,
                                           "(mean_p - mean_n) / sd_n");
        } else {
            r.signalToNoise = notComputed("non-zero negative control variance");
        }
    }

    // Signal over background: the positive control is the assay signal well and the
    // negative control is its background. Needs only the two means.
    if (pos.empty() || neg.empty()) {
        r.signalToBackground = notComputed("positive and negative control wells");
    } else if (mean(neg) == 0.0) {
        // A ratio against a zero background is not "no signal", it is undefined;
        // saying which is missing is the difference between the two.
        r.signalToBackground = notComputed("a non-zero background (negative control) mean");
    } else {
        r.signalToBackground = makeQuantity(mean(pos) / mean(neg), "", 0.0, Provenance::Measured,
                                            "mean_p / mean_n");
    }

    // Edge effect: outer ring against the interior, Mann-Whitney. Reported only.
    std::vector<double> ring;
    std::vector<double> interior;
    std::vector<std::vector<double>> byRow(static_cast<std::size_t>(std::max(plate.rows, 0)));
    std::vector<std::vector<double>> byCol(static_cast<std::size_t>(std::max(plate.columns, 0)));
    for (const Well& w : plate.wells) {
        if (!isPopulation(w)) continue;
        if (w.row >= 0 && w.row < plate.rows) byRow[static_cast<std::size_t>(w.row)].push_back(w.readout);
        if (w.column >= 0 && w.column < plate.columns)
            byCol[static_cast<std::size_t>(w.column)].push_back(w.readout);
        const bool onRing = w.row == 0 || w.column == 0 || w.row == plate.rows - 1 ||
                            w.column == plate.columns - 1;
        (onRing ? ring : interior).push_back(w.readout);
    }
    if (plate.rows >= 3 && plate.columns >= 3 && ring.size() >= 2 && interior.size() >= 2) {
        const double p = mannWhitneyP(ring, interior);
        r.edgeEffectP = makeQuantity(p, "", 0.0, Provenance::Measured,
                                     "Mann-Whitney U, outer ring vs interior, normal "
                                     "approximation with tie correction");
        if (p < 0.01) {
            r.warnings.push_back("edge effect detected (Mann-Whitney p = " + std::to_string(p) +
                                 "); reported, not corrected");
        }
    } else {
        r.edgeEffectP = notComputed("a plate of at least 3x3 with interior sample wells");
    }

    const double rowP = kruskalWallisP(byRow);
    r.rowEffectP = std::isnan(rowP) ? notComputed("at least two rows with two sample wells each")
                                    : makeQuantity(rowP, "", 0.0, Provenance::Measured,
                                                   "Kruskal-Wallis H across rows, tie-corrected");
    const double colP = kruskalWallisP(byCol);
    r.columnEffectP = std::isnan(colP)
                          ? notComputed("at least two columns with two sample wells each")
                          : makeQuantity(colP, "", 0.0, Provenance::Measured,
                                         "Kruskal-Wallis H across columns, tie-corrected");
    if (!std::isnan(rowP) && rowP < 0.01) {
        r.warnings.push_back("row effect detected (Kruskal-Wallis p = " + std::to_string(rowP) +
                             "); reported, not corrected");
    }
    if (!std::isnan(colP) && colP < 0.01) {
        r.warnings.push_back("column effect detected (Kruskal-Wallis p = " + std::to_string(colP) +
                             "); reported, not corrected");
    }

    // The bands are the published ones and are stated, not collapsed into a colour.
    std::string band = "Z-prime not available";
    if (r.zPrime.provenance == Provenance::Measured) {
        if (r.zPrime.value >= 0.5) {
            band = "Z-prime >= 0.5: excellent separation";
        } else if (r.zPrime.value > 0.0) {
            band = "Z-prime between 0 and 0.5: marginal, hit calls are unreliable";
        } else {
            band = "Z-prime <= 0: unusable, the control distributions overlap";
        }
    }
    r.interpretation = band +
                       " (bands: >= 0.5 excellent, 0 to 0.5 marginal, <= 0 unusable). Edge, row "
                       "and column effects are reported as p-values and never auto-corrected.";
    return r;
}

}  // namespace biocad::assay
