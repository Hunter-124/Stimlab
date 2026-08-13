#include "pkpd/PkEngine.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "numeric/Ode.h"

namespace biocad::pkpd {
namespace {

constexpr double kLn2 = 0.6931471805599453;

// Formats one number the way the assumptions block reads best: enough digits to be
// reproducible, not so many that a stated default looks measured.
std::string fmt(double v, int decimals = 3) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*f", decimals, v);
    return buf;
}

// A parameter whose provenance is not Measured must be visible as an assumption.
// `why` explains the tier, not the number.
void noteAssumption(std::vector<std::string>& out, const char* symbol, const Quantity& q,
                    const char* why, int decimals = 3) {
    if (q.provenance == Provenance::Measured) return;
    std::string line = std::string(symbol) + " = " + fmt(q.value, decimals);
    if (!q.unit.empty()) line += " " + q.unit;
    line += " (" + std::string(provenanceLabel(q.provenance)) + "; " + why + ")";
    if (!q.source.empty()) line += " [" + q.source + "]";
    out.push_back(line);
    // A predicted clearance carries the only honest generalisation figure available.
    if (q.provenance == Provenance::Predicted && std::string(symbol) == "CL") {
        out.push_back("Predicted clearance generalises poorly: only 10.4% of "
                      "high-clearance compounds are predicted within 2-fold "
                      "(PKSmart, Seal et al. 2025, J Cheminform 17:147).");
    }
}

bool isTwoCompartment(PkModel m) { return m == PkModel::OralTwoCompartment; }
bool isOral(PkModel m) {
    return m == PkModel::OralOneCompartment || m == PkModel::OralTwoCompartment;
}

// Index of the central compartment in the state vector.
std::size_t centralIndex(PkModel m) { return isOral(m) ? 1u : 0u; }

std::size_t stateSize(PkModel m) {
    if (isTwoCompartment(m)) return 3;   // depot, central, peripheral
    if (isOral(m)) return 2;             // depot, central
    return 1;                            // central
}

}  // namespace

double batemanConcentration(double dose, double F, double ka, double ke, double V,
                            double t) {
    if (V <= 0.0 || t < 0.0) return 0.0;
    const double amount = F * dose;
    // The 0/0 boundary. The limit of the Bateman function as ka -> ke is finite, so
    // the flip-flop boundary is a value, not an error.
    if (std::fabs(ka - ke) < 1e-12) {
        return amount * ke * t * std::exp(-ke * t) / V;
    }
    return amount * ka / (V * (ka - ke)) * (std::exp(-ke * t) - std::exp(-ka * t));
}

double accumulationRatio(double ke, double tauH) {
    if (ke <= 0.0 || tauH <= 0.0) return 0.0;
    return 1.0 / (1.0 - std::exp(-ke * tauH));
}

double steadyStateAverage(double F, double doseMg, double clearance, double tauH) {
    if (clearance <= 0.0 || tauH <= 0.0) return 0.0;
    return F * doseMg / (clearance * tauH);
}

PkProfile simulate(const PkModelSpec& spec, const DoseRegimen& regimen) {
    PkProfile profile;

    const double V = spec.volume.value;
    const double CL = spec.clearance.value;
    const double ka = spec.absorptionRate.value;
    const double F = isOral(spec.model) ? spec.bioavailability.value : 1.0;
    const double fu = spec.unboundFraction.value;
    const double V2 = spec.volumePeripheral.value;
    const double Q = spec.intercompartmental.value;
    const bool michaelisMenten = spec.vmax.value > 0.0;
    const double vmax = spec.vmax.value;
    const double km = spec.km.value;

    if (V <= 0.0) {
        profile.cmax = notComputed("central volume of distribution V");
        profile.tmax = notComputed("central volume of distribution V");
        profile.auc = notComputed("central volume of distribution V");
        profile.halfLife = notComputed("central volume of distribution V");
        profile.accumulation = notComputed("central volume of distribution V");
        profile.note = "No profile: the central volume V is required to convert an "
                       "amount into a concentration.";
        return profile;
    }
    const double ke = (CL > 0.0) ? CL / V : 0.0;

    // ---- integration -----------------------------------------------------
    const double step = spec.stepH > 0.0 ? spec.stepH : 0.01;
    const double horizon = spec.horizonH > 0.0 ? spec.horizonH : 24.0;
    const std::size_t n = stateSize(spec.model);
    const std::size_t ci = centralIndex(spec.model);

    // Infusion windows are rate terms; bolus/oral doses are state additions applied
    // at a segment boundary. Splitting the horizon at every event time is what keeps
    // a bolus from being averaged into an RK4 step.
    struct Infusion { double start, end, rateMgPerH; };
    std::vector<Infusion> infusions;
    std::vector<std::pair<double, double>> instant;  // (time, mg)
    for (const auto& d : regimen.doses) {
        if (d.timeH > horizon || d.amountMg <= 0.0) continue;
        if (spec.model == PkModel::IvInfusion && d.durationH > 0.0) {
            infusions.push_back({d.timeH, d.timeH + d.durationH,
                                 d.amountMg / d.durationH});
        } else {
            instant.emplace_back(d.timeH, d.amountMg);
        }
    }
    std::sort(instant.begin(), instant.end());

    std::vector<double> breakpoints{0.0, horizon};
    for (const auto& iv : infusions) {
        breakpoints.push_back(iv.start);
        breakpoints.push_back(std::min(iv.end, horizon));
    }
    for (const auto& d : instant) breakpoints.push_back(d.first);
    std::sort(breakpoints.begin(), breakpoints.end());
    breakpoints.erase(std::remove_if(breakpoints.begin(), breakpoints.end(),
                                     [&](double t) { return t < 0.0 || t > horizon; }),
                      breakpoints.end());
    breakpoints.erase(std::unique(breakpoints.begin(), breakpoints.end(),
                                  [](double a, double b) {
                                      return std::fabs(a - b) < 1e-12;
                                  }),
                      breakpoints.end());

    // The infusion rate is held constant across a segment rather than tested against
    // `t` inside the derivative: the horizon is already split at every window edge, and
    // a time test would let a boundary RK4 stage sample the wrong side of the switch
    // (a k1 evaluated exactly at the window's end still saw the pump running, which
    // over-delivered one sixth of a step's worth of drug).
    double infusionRate = 0.0;
    const numeric::OdeDerivative f = [&](double /*t*/, const std::vector<double>& y,
                                         std::vector<double>& dydt) {
        std::fill(dydt.begin(), dydt.end(), 0.0);
        const double amountCentral = y[ci];
        const double conc = amountCentral / V;
        const double elimination = michaelisMenten
                                       ? (km + conc > 0.0 ? vmax * conc / (km + conc) : 0.0)
                                       : CL * conc;
        dydt[ci] -= elimination;
        if (isOral(spec.model)) {
            const double absorbed = ka * y[0];
            dydt[0] -= absorbed;
            dydt[ci] += absorbed;
        }
        if (isTwoCompartment(spec.model) && V2 > 0.0 && Q > 0.0) {
            const double flux = Q * (conc - y[2] / V2);
            dydt[ci] -= flux;
            dydt[2] += flux;
        }
        dydt[ci] += infusionRate;
    };

    std::vector<double> y(n, 0.0);
    auto applyDosesAt = [&](double t) {
        for (const auto& d : instant) {
            if (std::fabs(d.first - t) < 1e-12) {
                // Oral input lands in the depot; an IV bolus lands in the central
                // compartment already dissolved, so F does not apply to it.
                if (isOral(spec.model)) y[0] += F * d.second;
                else                    y[ci] += d.second;
            }
        }
    };

    auto record = [&](double t) {
        profile.timeH.push_back(t);
        const double c = y[ci] / V;
        profile.concentrationMgPerL.push_back(c);
        profile.unboundMgPerL.push_back(fu > 0.0 ? fu * c : c);
    };

    applyDosesAt(0.0);
    record(0.0);
    for (std::size_t i = 0; i + 1 < breakpoints.size(); ++i) {
        const double t0 = breakpoints[i];
        const double t1 = breakpoints[i + 1];
        infusionRate = 0.0;
        const double mid = 0.5 * (t0 + t1);
        for (const auto& iv : infusions) {
            if (mid >= iv.start && mid <= iv.end) infusionRate += iv.rateMgPerH;
        }
        bool first = true;
        numeric::rk4Integrate(t0, t1, step, y, f,
                              [&](double t, const std::vector<double>&) {
                                  if (first) { first = false; return; }  // t0 already recorded
                                  profile.timeH.push_back(t);
                                  const double c = y[ci] / V;
                                  profile.concentrationMgPerL.push_back(c);
                                  profile.unboundMgPerL.push_back(fu > 0.0 ? fu * c : c);
                              });
        applyDosesAt(t1);
        if (!instant.empty()) {
            // A dose applied at t1 changes the concentration discontinuously; record
            // the post-dose value so the plotted curve shows the jump rather than
            // interpolating through it.
            bool dosed = false;
            for (const auto& d : instant) dosed = dosed || std::fabs(d.first - t1) < 1e-12;
            if (dosed) record(t1);
        }
    }

    // ---- summary ---------------------------------------------------------
    // A simulated exposure is a CONSTRUCTED ARTEFACT, so its tier is Model even when
    // every input parameter was measured: no plasma was assayed. It drops to
    // NotComputed when a required input is missing. It is deliberately never
    // Heuristic, because a concentration carries mg/L and the domain (rightly)
    // forbids a heuristic from carrying a physical unit - a parameter that is only
    // rank-ordering information is reported in `assumptions`, not laundered into the
    // tier of a number that has units.
    Provenance inputs = weakest(spec.clearance.provenance, spec.volume.provenance);
    if (isOral(spec.model)) {
        inputs = weakest(inputs, weakest(spec.bioavailability.provenance,
                                         spec.absorptionRate.provenance));
    }
    if (michaelisMenten) {
        inputs = weakest(inputs, weakest(spec.vmax.provenance, spec.km.provenance));
    }
    const Provenance tier =
        inputs == Provenance::NotComputed ? Provenance::NotComputed : Provenance::Model;
    const std::string engine =
        "RK4, step " + fmt(step, 4) + " h over " + fmt(horizon, 2) + " h";

    std::size_t peak = 0;
    for (std::size_t i = 1; i < profile.concentrationMgPerL.size(); ++i) {
        if (profile.concentrationMgPerL[i] > profile.concentrationMgPerL[peak]) peak = i;
    }
    if (!profile.concentrationMgPerL.empty()) {
        profile.cmax = makeQuantity(profile.concentrationMgPerL[peak], "mg/L", 0.0, tier,
                                    engine);
        profile.tmax = makeQuantity(profile.timeH[peak], "h", 0.0, tier, engine);
        profile.auc = makeQuantity(
            numeric::trapezoid(profile.timeH, profile.concentrationMgPerL), "mg*h/L", 0.0,
            tier, "trapezoidal over " + fmt(horizon, 2) + " h; " + engine);
    } else {
        profile.cmax = notComputed("a simulated concentration series");
        profile.tmax = notComputed("a simulated concentration series");
        profile.auc = notComputed("a simulated concentration series");
    }

    if (michaelisMenten) {
        // A saturable system has no single half-life: the apparent one changes with
        // concentration, so reporting a number here would be a fabrication.
        profile.halfLife = notComputed("linear elimination");
    } else if (CL > 0.0) {
        profile.halfLife = makeQuantity(kLn2 * V / CL, "h", 0.0, tier, "ln2 * V / CL");
    } else {
        profile.halfLife = notComputed("a positive clearance CL");
    }

    // Rac needs a repeated, evenly spaced regimen and a first-order ke.
    double tau = 0.0;
    bool evenlySpaced = regimen.doses.size() >= 2;
    if (evenlySpaced) {
        std::vector<double> times;
        for (const auto& d : regimen.doses) times.push_back(d.timeH);
        std::sort(times.begin(), times.end());
        tau = times[1] - times[0];
        for (std::size_t i = 1; i + 1 < times.size(); ++i) {
            if (std::fabs((times[i + 1] - times[i]) - tau) > 1e-9) evenlySpaced = false;
        }
        if (tau <= 0.0) evenlySpaced = false;
    }
    if (!evenlySpaced) {
        profile.accumulation = notComputed("a repeated regimen of evenly spaced doses");
    } else if (michaelisMenten) {
        profile.accumulation = notComputed("linear elimination");
    } else if (ke <= 0.0) {
        profile.accumulation = notComputed("a positive elimination rate constant ke");
    } else {
        profile.accumulation = makeQuantity(
            accumulationRatio(ke, tau), "", 0.0, tier,
            "Rac = 1 / (1 - exp(-ke*tau)), tau = " + fmt(tau, 2) + " h");
    }

    profile.flipFlop = isOral(spec.model) && ka > 0.0 && ke > 0.0 && ka < ke;

    // ---- assumptions -----------------------------------------------------
    if (isOral(spec.model)) {
        noteAssumption(profile.assumptions, "F", spec.bioavailability,
                       "bioavailability has no credible structure-only predictor", 2);
        noteAssumption(profile.assumptions, "ka", spec.absorptionRate,
                       "the absorption rate constant has no credible structure-only "
                       "predictor");
    }
    noteAssumption(profile.assumptions, "CL", spec.clearance, "clearance was not measured "
                                                              "in this subject");
    noteAssumption(profile.assumptions, "V", spec.volume,
                   "the central volume was not measured in this subject");
    if (isTwoCompartment(spec.model)) {
        noteAssumption(profile.assumptions, "V2", spec.volumePeripheral,
                       "the peripheral volume was not measured in this subject");
        noteAssumption(profile.assumptions, "Q", spec.intercompartmental,
                       "intercompartmental clearance was not measured in this subject");
    }
    noteAssumption(profile.assumptions, "fu", spec.unboundFraction,
                   "the unbound fraction has no credible structure-only predictor");
    if (michaelisMenten) {
        noteAssumption(profile.assumptions, "Vmax", spec.vmax,
                       "saturable elimination capacity was not measured in this subject");
        noteAssumption(profile.assumptions, "Km", spec.km,
                       "the half-saturation concentration was not measured in this "
                       "subject");
    }
    profile.assumptions.push_back(
        "This is an exposure scenario under the parameters above, not a dose "
        "recommendation.");

    profile.note = "Integrated with " + engine + ".";
    if (michaelisMenten) {
        profile.note += " Elimination is Michaelis-Menten (Vmax/Km), so exposure is "
                        "not proportional to dose and there is no single half-life.";
    }
    if (profile.flipFlop) {
        profile.note += " Flip-flop kinetics: ka < ke, so the terminal phase reflects "
                        "absorption, not elimination, and a half-life read off that "
                        "slope is an absorption half-life.";
    }
    return profile;
}

OccupancyCurve occupancy(const PkProfile& profile, const Quantity& kd) {
    OccupancyCurve curve;
    if (kd.provenance == Provenance::NotComputed || kd.value <= 0.0) {
        curve.peakOccupancy = notComputed("a positive measured Kd");
        curve.timeAbove50Pct = notComputed("a positive measured Kd");
        curve.note = "No occupancy curve: a positive Kd for this target is required, "
                     "and occupancy computed from a guessed Kd would be fiction.";
        return curve;
    }
    if (profile.unboundMgPerL.empty()) {
        curve.peakOccupancy = notComputed("a simulated unbound concentration series");
        curve.timeAbove50Pct = notComputed("a simulated unbound concentration series");
        curve.note = "No occupancy curve: the profile has no unbound concentration "
                     "series.";
        return curve;
    }

    const double kdValue = kd.value;
    curve.timeH = profile.timeH;
    curve.occupancy.reserve(profile.unboundMgPerL.size());
    double peak = 0.0;
    for (double cu : profile.unboundMgPerL) {
        const double theta = cu / (kdValue + cu);
        curve.occupancy.push_back(theta);
        peak = std::max(peak, theta);
    }

    // Time above 50% is measured with linear interpolation of the crossings rather
    // than by counting samples, so it does not depend on the integration step.
    double above = 0.0;
    for (std::size_t i = 0; i + 1 < curve.occupancy.size(); ++i) {
        const double t0 = curve.timeH[i], t1 = curve.timeH[i + 1];
        const double a = curve.occupancy[i], b = curve.occupancy[i + 1];
        const double dt = t1 - t0;
        if (dt <= 0.0) continue;
        if (a >= 0.5 && b >= 0.5) {
            above += dt;
        } else if (a >= 0.5 || b >= 0.5) {
            const double frac = (0.5 - a) / (b - a);   // crossing fraction of the step
            above += (a >= 0.5) ? frac * dt : (1.0 - frac) * dt;
        }
    }

    // Same rule as the profile: occupancy is a constructed artefact, not a measurement,
    // and it carries units, so it is Model rather than the Kd's own tier.
    const Provenance tier = Provenance::Model;
    curve.peakOccupancy = makeQuantity(peak, "fraction", 0.0, tier,
                                       "theta = Cu / (Kd + Cu); Kd " + kd.source);
    curve.timeAbove50Pct = makeQuantity(above, "h", 0.0, tier,
                                        "linear interpolation of the 50% crossings");
    curve.note = "Fractional occupancy from the unbound concentration and the supplied "
                 "Kd. Occupancy is a target-engagement scenario, not an effect and not "
                 "a dose recommendation.";
    return curve;
}

}  // namespace biocad::pkpd
