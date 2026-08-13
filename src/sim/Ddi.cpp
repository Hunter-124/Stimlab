#include "sim/Ddi.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <limits>
#include <vector>

#include "core/Physiology.h"
#include "numeric/Ode.h"

namespace biocad::sim {
namespace {

std::string fmt(double v) {
    char b[48];
    std::snprintf(b, sizeof b, "%.6g", v);
    return b;
}

// The pack key for an enzyme's hepatic turnover. The pack keys are explicit
// ("CYP3A4_hepatic") because intestinal CYP3A4 turns over faster than hepatic and a
// single "CYP3A4" key would silently pick one of them.
std::string hepaticKey(const std::string& enzyme) { return enzyme + "_hepatic"; }
std::string intestinalKey(const std::string& enzyme) { return enzyme + "_intestinal"; }

// kobs, the observed inactivation rate: kinact*[I]/(KI+[I]). Zero when either TDI
// parameter is absent, which makes the TDI factor exactly 1 and leaves the other
// mechanisms untouched.
double kObs(double kinact, double kI, double conc) {
    if (kinact <= 0 || kI <= 0 || conc <= 0) return 0.0;
    return kinact * conc / (kI + conc);
}

// The induction factor 1 + d*Emax*[I]/(EC50+[I]).
double inductionFactor(double d, double emax, double ec50, double conc) {
    if (emax <= 0 || ec50 <= 0 || conc <= 0) return 1.0;
    return 1.0 + d * emax * conc / (ec50 + conc);
}

// The reversible factor 1/(1 + [I]/Ki).
double reversibleFactor(double ki, double conc) {
    if (ki <= 0 || conc <= 0) return 1.0;
    return 1.0 / (1.0 + conc / ki);
}

}  // namespace

InteractionReport interaction(const PerpetratorSpec& perpetrator, const VictimSpec& victim) {
    InteractionReport r;
    r.perpetrator = perpetrator;
    r.victim = victim;

    const auto& phys = core::physiology();
    const double kdegHep = core::enzymeDegradationRate(hepaticKey(perpetrator.enzyme));
    const double kdegGut = core::enzymeDegradationRate(intestinalKey(perpetrator.enzyme));

    const double iHep = perpetrator.unboundHepaticInletUM;
    const double iGut = perpetrator.enterocyteUM;
    const double ki = perpetrator.ki;

    // ---- 1. FDA basic-model screening R-values -----------------------------
    if (ki > 0 && iHep > 0)
        r.r1 = makeQuantity(1.0 + iHep / ki, "", 0, Provenance::Predicted,
                            "FDA basic model R1 = 1 + [I]h,u/Ki, from the in vitro Ki in " +
                                perpetrator.source);
    else
        r.r1 = notComputed(ki > 0 ? "an unbound hepatic inlet concentration [I]h,u"
                                  : "a reversible inhibition constant Ki");

    if (ki > 0 && iGut > 0)
        r.r1Gut = makeQuantity(1.0 + iGut / ki, "", 0, Provenance::Predicted,
                               "FDA basic model R1,gut = 1 + [I]g/Ki");
    else
        r.r1Gut = notComputed(ki > 0 ? "an enterocyte concentration [I]g"
                                     : "a reversible inhibition constant Ki");

    const double iSys = perpetrator.unboundSystemicUM > 0 ? perpetrator.unboundSystemicUM : iHep;
    if (perpetrator.kinact > 0 && perpetrator.kI > 0 && iSys > 0 && kdegHep > 0) {
        const double r2 = (kdegHep + kObs(perpetrator.kinact, perpetrator.kI, iSys)) / kdegHep;
        r.r2 = makeQuantity(r2, "", 0, Provenance::Predicted,
                            "FDA basic model R2 = (kdeg + kobs)/kdeg with kdeg = " +
                                fmt(kdegHep) + " /h from assets/packs/physiology.json");
    } else {
        r.r2 = notComputed(kdegHep > 0 ? "time-dependent inactivation parameters kinact and KI"
                                       : "a kdeg for " + hepaticKey(perpetrator.enzyme) +
                                             " in assets/packs/physiology.json");
    }

    if (perpetrator.indEmax > 0 && perpetrator.indEc50 > 0 && iSys > 0)
        r.rInduction = makeQuantity(
            1.0 / inductionFactor(perpetrator.indD, perpetrator.indEmax, perpetrator.indEc50,
                                  iSys),
            "", 0, Provenance::Predicted,
            "FDA basic model R3 = 1/(1 + d*Emax*[I]/(EC50+[I])); a value at or below 0.8 is "
            "the screening cut for induction");
    else
        r.rInduction = notComputed("induction parameters Emax and EC50");

    // ---- 2. The mechanistic static AUCR ------------------------------------
    const double fm = victim.fractionMetabolizedByEnzyme;
    const double fg = victim.intestinalAvailability;

    if (fm < 0 || fm > 1) {
        // The whole point of the check: without fm there IS no AUCR, and a default
        // of 1 would silently claim the victim is cleared entirely by this enzyme.
        r.aucRatio = notComputed("fm");
        r.aucRatioHepaticOnly = notComputed("fm");
        r.theoreticalCeiling = notComputed("fm");
        r.warnings.push_back("fm, the fraction of the victim's clearance carried by " +
                             perpetrator.enzyme +
                             ", was not supplied: the mechanistic static AUCR is not computed. "
                             "Assuming fm = 1 is what turns a 1.3-fold interaction into a "
                             "5-fold one.");
        r.assumptions.push_back("physiology: " + phys.path);
        return r;
    }

    const double aHep = (kdegHep > 0)
                            ? kdegHep / (kdegHep + kObs(perpetrator.kinact, perpetrator.kI, iHep))
                            : 1.0;
    const double bHep = inductionFactor(perpetrator.indD, perpetrator.indEmax,
                                        perpetrator.indEc50, iHep);
    const double cHep = reversibleFactor(ki, iHep);
    const double hepaticTerm = aHep * bHep * cHep;
    const double aucrHepatic = 1.0 / (hepaticTerm * fm + (1.0 - fm));

    r.aucRatioHepaticOnly =
        makeQuantity(aucrHepatic, "", 0, Provenance::Predicted,
                     "mechanistic static model, hepatic term only: "
                     "1/((A*B*C)*fm + (1-fm)) with A = " +
                         fmt(aHep) + " (TDI), B = " + fmt(bHep) + " (induction), C = " +
                         fmt(cHep) + " (reversible)");
    r.theoreticalCeiling =
        makeQuantity(fm < 1.0 ? 1.0 / (1.0 - fm) : std::numeric_limits<double>::infinity(), "",
                     0, Provenance::Predicted,
                     "1/(1-fm): the AUC ratio at complete inhibition of this pathway, which no "
                     "inhibitor of it can exceed");

    if (fg >= 0 && fg <= 1) {
        const double aGut =
            (kdegGut > 0) ? kdegGut / (kdegGut + kObs(perpetrator.kinact, perpetrator.kI, iGut))
                          : 1.0;
        const double bGut = inductionFactor(perpetrator.indD, perpetrator.indEmax,
                                            perpetrator.indEc50, iGut);
        const double cGut = reversibleFactor(ki, iGut);
        const double gutTerm = aGut * bGut * cGut;
        const double aucr = aucrHepatic / (gutTerm * (1.0 - fg) + fg);
        r.gutIncluded = true;
        r.aucRatio = makeQuantity(aucr, "", 0, Provenance::Predicted,
                                  "mechanistic static model, hepatic x gut: the hepatic term "
                                  "above times 1/((Ag*Bg*Cg)*(1-Fg) + Fg)");
        if (kdegGut <= 0 && perpetrator.kinact > 0)
            r.warnings.push_back("no intestinal kdeg for " + perpetrator.enzyme +
                                 " in the physiology pack, so the gut TDI term was left at 1");
        r.assumptions.push_back("Qen = " + fmt(phys.enterocyteBloodFlowLPerH) +
                                " L/h from the physiology pack bounds the gut-wall term");
    } else {
        r.gutIncluded = false;
        r.aucRatio = notComputed("Fg, the victim's intestinal availability");
        r.warnings.push_back("Fg was not supplied, so only the hepatic AUC ratio is reported; "
                             "for a CYP3A4 substrate with a low Fg the gut term can be the "
                             "larger of the two");
    }

    // The dominant mechanism is named by comparing the three factors' distance from
    // 1, because an AUCR alone cannot tell a reader whether to expect a slow onset
    // (inactivation) or an immediate one (reversible).
    const double dRev = std::fabs(1.0 - cHep), dTdi = std::fabs(1.0 - aHep),
                 dInd = std::fabs(1.0 - bHep);
    if (dRev == 0 && dTdi == 0 && dInd == 0)
        r.dominantMechanism = "none: no mechanism had usable in vitro parameters";
    else if (dRev >= dTdi && dRev >= dInd)
        r.dominantMechanism = "reversible inhibition (immediate onset and offset)";
    else if (dTdi >= dInd)
        r.dominantMechanism = "time-dependent inactivation (onset and offset over several "
                              "enzyme half-lives, ln2/kdeg = " +
                              fmt(kdegHep > 0 ? std::log(2.0) / kdegHep : 0.0) + " h)";
    else
        r.dominantMechanism = "induction (onset and offset over several enzyme half-lives)";

    r.assumptions.push_back("Qh = " + fmt(phys.hepaticBloodFlowLPerH) +
                            " L/h, kdeg(" + hepaticKey(perpetrator.enzyme) + ") = " +
                            fmt(kdegHep) + " /h, both from " + phys.path);
    r.assumptions.push_back("fm = " + fmt(fm) + " and Fg = " +
                            (fg >= 0 ? fmt(fg) : std::string("not supplied")) +
                            " are victim INPUTS, not predictions");
    r.assumptions.push_back(
        "this is an exposure ratio under the stated in vitro parameters; it is not a dose, a "
        "dose adjustment or a risk category");
    if (!phys.loaded) r.warnings.push_back(phys.error);
    return r;
}

EnzymeTimeCourse enzymeTimeCourse(const PerpetratorSpec& perpetrator, double horizonH) {
    EnzymeTimeCourse tc;
    const double kdeg = core::enzymeDegradationRate(hepaticKey(perpetrator.enzyme));
    tc.kdegUsed = kdeg;
    tc.kdegSource = core::physiology().path + " (" + hepaticKey(perpetrator.enzyme) + ")";
    if (kdeg <= 0) {
        tc.assumptions.push_back("no kdeg for " + hepaticKey(perpetrator.enzyme) +
                                 " in the physiology pack, so no dynamic enzyme model was run");
        return tc;
    }

    const double conc = perpetrator.unboundHepaticInletUM > 0 ? perpetrator.unboundHepaticInletUM
                                                              : 0.0;
    const double kobs = kObs(perpetrator.kinact, perpetrator.kI, conc);
    const double induction = inductionFactor(perpetrator.indD, perpetrator.indEmax,
                                             perpetrator.indEc50, conc);
    const double reversible = reversibleFactor(perpetrator.ki, conc);

    // dE/dt = kdeg*(1 + d*Emax*I/(EC50+I)) - kdeg*E - kobs*E, so E(0) = 1 and the
    // approach to steady state has rate constant (kdeg + kobs).
    const numeric::OdeDerivative f = [&](double, const std::vector<double>& y,
                                         std::vector<double>& dydt) {
        dydt[0] = kdeg * induction - kdeg * y[0] - kobs * y[0];
    };

    std::vector<double> y{1.0};
    const double rate = kdeg + kobs;
    // A step small against the fastest time constant; 1/200th of it keeps the RK4
    // error far below the 1e-9 agreement the static comparison is checked at.
    const double step = std::min(std::max(horizonH, 1.0) / 200.0, 1.0 / (200.0 * rate));
    if (horizonH > 0) {
        numeric::rk4Integrate(0.0, horizonH, step, y, f,
                     [&](double t, const std::vector<double>& state) {
                         tc.timeH.push_back(t);
                         tc.relativeActivity.push_back(state[0]);
                         tc.clearanceRatio.push_back(state[0] * reversible);
                     });
    }

    // The steady state is integrated separately, out to 40 time constants, because
    // the requested horizon is a DISPLAY window: a 24 h plot of a 36 h-turnover
    // enzyme has not reached steady state, and comparing an unconverged value with
    // the static model would report a disagreement that is the plot's fault.
    std::vector<double> ySs{1.0};
    const double tSs = 40.0 / rate;
    numeric::rk4Integrate(0.0, tSs, std::min(tSs / 20000.0, 1.0 / (200.0 * rate)), ySs, f,
                 [](double, const std::vector<double>&) {});
    tc.steadyStateActivity = ySs[0];

    // The static model's equivalent: E_ss = B / (1 + kobs/kdeg) = induction * A.
    tc.staticModelActivity = induction * kdeg / (kdeg + kobs);
    tc.agreement = std::fabs(tc.steadyStateActivity - tc.staticModelActivity);

    tc.assumptions.push_back("constant unbound hepatic inlet concentration [I] = " + fmt(conc) +
                             " uM: a real perpetrator's concentration oscillates, and this "
                             "scenario deliberately does not");
    tc.assumptions.push_back("kdeg = " + fmt(kdeg) + " /h (turnover half-life " +
                             fmt(std::log(2.0) / kdeg) + " h) from " + tc.kdegSource);
    tc.assumptions.push_back(
        "CLint(t) = CLint0 * E(t) / (1 + [I]/Ki); the reversible term is instantaneous and "
        "the E(t) term is what carries the slow onset");
    tc.assumptions.push_back("steady-state activity " + fmt(tc.steadyStateActivity) +
                             " versus the static model's " + fmt(tc.staticModelActivity) +
                             ", difference " + fmt(tc.agreement));
    return tc;
}

ImpairmentScenario impairment(const VictimSpec& victim, double renalFunctionRatio,
                              double hepaticClintRatio) {
    ImpairmentScenario s;
    s.label = victim.label;
    s.renalFunctionRatio = renalFunctionRatio;
    s.hepaticClintRatio = hepaticClintRatio;
    s.boundaryStatement =
        "This is an exposure ratio for an editable scenario, not a dose, a dose adjustment or "
        "a statement about any individual's organ function. BioCAD does not convert a "
        "Child-Pugh class or a creatinine clearance into a clearance or a dose.";

    const double fe = victim.fractionExcretedUnchanged;
    if (fe < 0 || fe > 1) {
        s.exposureRatio = notComputed("fe, the fraction of the dose excreted unchanged");
        s.assumptions.push_back(
            "without fe the renal and non-renal routes cannot be separated, so no exposure "
            "ratio is computed");
        return s;
    }
    if (renalFunctionRatio <= 0 || hepaticClintRatio <= 0) {
        s.exposureRatio = notComputed("positive renal-function and hepatic-CLint ratios");
        return s;
    }

    // CL_impaired/CL_normal = fe*RF + (1-fe)*r under the low-extraction limit of the
    // well-stirred model, where CL_H tracks fu*CLint directly. For a HIGH-extraction
    // drug clearance is bounded by Qh and barely moves with CLint, which is why the
    // bound is quoted in the assumptions instead of being silently ignored.
    const double clRatio = fe * renalFunctionRatio + (1.0 - fe) * hepaticClintRatio;
    s.exposureRatio = makeQuantity(
        1.0 / clRatio, "", 0, Provenance::Predicted,
        "AUC ratio = 1 / (fe*renal function + (1-fe)*hepatic CLint ratio) = 1/" + fmt(clRatio));
    s.assumptions.push_back("fe = " + fmt(fe) + " is a victim input");
    s.assumptions.push_back(
        "low-extraction limit of the well-stirred model: hepatic clearance tracks fu*CLint. "
        "A high-extraction drug's clearance is bounded by Qh = " +
        fmt(core::physiology().hepaticBloodFlowLPerH) +
        " L/h and would move far less than this ratio suggests.");
    s.assumptions.push_back(
        "renal function ratio " + fmt(renalFunctionRatio) + " and hepatic CLint ratio " +
        fmt(hepaticClintRatio) + " are scenario inputs the user edits, not a diagnosis");
    return s;
}

}  // namespace biocad::sim
