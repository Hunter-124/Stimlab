// Tests for the pure pharmacodynamic fits. Every case checks a number against a closed
// form or an exact algebraic identity, not against the implementation's own output.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <vector>

#include "pkpd/Fits.h"

using namespace biocad;
using Catch::Matchers::ContainsSubstring;

namespace {

// Noiseless 4PL data, two points per decade over seven decades.
std::vector<DoseResponsePoint> syntheticCurve(double top, double bottom, double ec50, double nH) {
    std::vector<DoseResponsePoint> pts;
    for (int decade = -10; decade <= -4; ++decade) {
        for (double factor : {1.0, 3.0}) {
            const double c = factor * std::pow(10.0, decade);
            const double u = std::pow(ec50 / c, nH);  // u = (EC50/[A])^nH
            pts.push_back({c, bottom + (top - bottom) / (1.0 + u)});
        }
    }
    return pts;
}

// Schild data from a known KB: log10(DR-1) = slope * (log10[B] - log10 KB).
std::vector<SchildPoint> syntheticSchild(double kb, double slope) {
    std::vector<SchildPoint> pts;
    for (int decade = -8; decade <= -5; ++decade) {
        const double b = std::pow(10.0, decade);
        pts.push_back({b, 1.0 + std::pow(10.0, slope * (std::log10(b) - std::log10(kb)))});
    }
    return pts;
}

}  // namespace

TEST_CASE("4PL recovers known parameters from noiseless data", "[pkpd]") {
    const auto fit = pkpd::fitFourParameterLogistic(syntheticCurve(100.0, 0.0, 1e-7, 1.2));

    REQUIRE(fit.converged);
    REQUIRE(std::fabs(fit.top.value - 100.0) < 1e-9);
    REQUIRE(std::fabs(fit.bottom.value - 0.0) < 1e-9);
    REQUIRE(std::fabs(fit.ec50.value / 1e-7 - 1.0) < 1e-9);
    REQUIRE(std::fabs(fit.hillSlope.value - 1.2) < 1e-9);
    REQUIRE(fit.rSquared > 1.0 - 1e-12);

    // EC50 is a physical concentration; the asymptotes are on the assay's arbitrary
    // response scale, so they must be unitless Heuristic values.
    REQUIRE(fit.ec50.unit == "mol/L");
    REQUIRE(fit.ec50.provenance == Provenance::Predicted);
    REQUIRE(fit.top.unit.empty());
    REQUIRE(fit.top.provenance == Provenance::Heuristic);
    REQUIRE(fit.bottom.unit.empty());
    REQUIRE(fit.bottom.provenance == Provenance::Heuristic);
    REQUIRE(fit.hillSlope.unit.empty());
    REQUIRE(fit.hillSlope.provenance == Provenance::Predicted);
    REQUIRE_THAT(fit.ec50.source, ContainsSubstring("n=14 points"));
}

TEST_CASE("4PL with 1/y^2 weighting still recovers the truth", "[pkpd]") {
    const auto fit = pkpd::fitFourParameterLogistic(syntheticCurve(100.0, 1.0, 5e-8, 0.9), true);
    REQUIRE(fit.converged);
    REQUIRE(std::fabs(fit.ec50.value / 5e-8 - 1.0) < 1e-8);
    REQUIRE(std::fabs(fit.hillSlope.value - 0.9) < 1e-8);
    REQUIRE_THAT(fit.note, ContainsSubstring("1/y^2 weighting"));
}

TEST_CASE("4PL rejects an underdetermined or unloggable data set", "[pkpd]") {
    const auto full = syntheticCurve(100.0, 0.0, 1e-7, 1.2);

    SECTION("three points cannot determine four parameters") {
        const auto fit = pkpd::fitFourParameterLogistic({full[0], full[1], full[2]});
        REQUIRE_FALSE(fit.converged);
        REQUIRE(fit.top.provenance == Provenance::NotComputed);
        REQUIRE(fit.bottom.provenance == Provenance::NotComputed);
        REQUIRE(fit.ec50.provenance == Provenance::NotComputed);
        REQUIRE(fit.hillSlope.provenance == Provenance::NotComputed);
        REQUIRE_THAT(fit.note, ContainsSubstring("at least 4 points"));
    }

    SECTION("a zero concentration has no logarithm") {
        const auto fit = pkpd::fitFourParameterLogistic(
            {{0.0, 1.0}, {1e-9, 2.0}, {1e-8, 3.0}, {1e-7, 4.0}});
        REQUIRE_FALSE(fit.converged);
        REQUIRE(fit.ec50.provenance == Provenance::NotComputed);
        REQUIRE_THAT(fit.note, ContainsSubstring("non-positive concentration"));
    }
}

TEST_CASE("Cheng-Prusoff applies the modality's own equation", "[pkpd]") {
    ChengPrusoffInput in;
    in.ic50 = 1e-7;
    in.substrate = 1e-5;
    in.km = 1e-5;

    SECTION("competitive at [S] = Km halves the IC50") {
        in.modality = InhibitionModality::Competitive;
        const auto ki = pkpd::kiFromIc50(in);
        REQUIRE(ki.value == 5e-8);
        REQUIRE(ki.unit == "mol/L");
        REQUIRE(ki.provenance == Provenance::Predicted);
        REQUIRE_THAT(ki.source, ContainsSubstring("competitive"));
    }

    SECTION("uncompetitive at [S] = Km also halves it") {
        in.modality = InhibitionModality::Uncompetitive;
        REQUIRE(pkpd::kiFromIc50(in).value == 5e-8);
    }

    SECTION("noncompetitive needs no substrate correction") {
        in.modality = InhibitionModality::Noncompetitive;
        REQUIRE(pkpd::kiFromIc50(in).value == 1e-7);
    }

    SECTION("radioligand binding uses [L*] and Kd") {
        in.modality = InhibitionModality::RadioligandBinding;
        in.radioligand = 1e-9;
        in.kdRadioligand = 1e-9;
        REQUIRE(pkpd::kiFromIc50(in).value == 5e-8);
    }
}

TEST_CASE("Cheng-Prusoff modalities diverge by orders of magnitude", "[pkpd]") {
    // Ki_uncompetitive / Ki_competitive = (1 + [S]/Km) / (1 + Km/[S]): 10x at [S] = 10*Km
    // and 100x at [S] = 100*Km. Assuming the wrong modality therefore fabricates one to
    // two orders of magnitude, which is why the modality is required input.
    ChengPrusoffInput in;
    in.ic50 = 1e-7;
    in.km = 1e-5;
    in.substrate = 10.0 * in.km;

    in.modality = InhibitionModality::Competitive;
    const double competitive = pkpd::kiFromIc50(in).value;
    in.modality = InhibitionModality::Uncompetitive;
    const double uncompetitive = pkpd::kiFromIc50(in).value;

    REQUIRE(std::fabs(competitive - 1e-7 / 11.0) < 1e-20);
    REQUIRE(std::fabs(uncompetitive - 1e-7 / 1.1) < 1e-20);
    REQUIRE(std::fabs(uncompetitive / competitive - 10.0) < 1e-9);

    in.substrate = 100.0 * in.km;
    in.modality = InhibitionModality::Competitive;
    const double competitiveHigh = pkpd::kiFromIc50(in).value;
    in.modality = InhibitionModality::Uncompetitive;
    const double uncompetitiveHigh = pkpd::kiFromIc50(in).value;
    REQUIRE(std::fabs(uncompetitiveHigh / competitiveHigh - 100.0) < 1e-6);
}

TEST_CASE("Cheng-Prusoff names the missing prerequisite", "[pkpd]") {
    ChengPrusoffInput in;
    in.modality = InhibitionModality::Competitive;
    in.ic50 = 1e-7;
    in.substrate = 1e-5;  // km left absent

    const auto ki = pkpd::kiFromIc50(in);
    REQUIRE(ki.provenance == Provenance::NotComputed);
    REQUIRE_THAT(ki.source, ContainsSubstring("km"));

    SECTION("a missing substrate is named too") {
        in.km = 1e-5;
        in.substrate = -1.0;
        REQUIRE_THAT(pkpd::kiFromIc50(in).source, ContainsSubstring("substrate"));
    }

    SECTION("a radioligand fit names its own fields") {
        in.modality = InhibitionModality::RadioligandBinding;
        REQUIRE_THAT(pkpd::kiFromIc50(in).source, ContainsSubstring("radioligand"));
    }
}

TEST_CASE("Cheng-Prusoff shows the depletion-corrected Ki when binding is tight", "[pkpd]") {
    ChengPrusoffInput in;
    in.modality = InhibitionModality::Noncompetitive;
    in.ic50 = 1e-9;
    in.enzymeConc = 5e-10;

    const auto ki = pkpd::kiFromIc50(in);
    REQUIRE(ki.value == 1e-9);
    REQUIRE_THAT(ki.source, ContainsSubstring("[I] >> [E]t"));
    REQUIRE_THAT(ki.source, ContainsSubstring("depletion-corrected"));
    REQUIRE_THAT(ki.source, ContainsSubstring("7.5e-10"));
}

TEST_CASE("Schild recovers pA2 and only reports KB at unit slope", "[pkpd]") {
    SECTION("unit slope") {
        const auto r = pkpd::schild(syntheticSchild(1e-8, 1.0));
        REQUIRE(std::fabs(r.pA2.value - 8.0) < 1e-6);
        REQUIRE(std::fabs(r.slope.value - 1.0) < 1e-9);
        REQUIRE(r.slopeCiLow <= 1.0);
        REQUIRE(r.slopeCiHigh >= 1.0);
        REQUIRE(r.kbUsable);
        REQUIRE(r.kb.unit == "mol/L");
        REQUIRE(std::fabs(r.kb.value / 1e-8 - 1.0) < 1e-6);
    }

    SECTION("slope 0.5 makes KB meaningless") {
        const auto r = pkpd::schild(syntheticSchild(1e-8, 0.5));
        REQUIRE(std::fabs(r.slope.value - 0.5) < 1e-9);
        REQUIRE_FALSE(r.kbUsable);
        REQUIRE(r.kb.provenance == Provenance::NotComputed);
        REQUIRE(r.kb.source == "unit Schild slope");
        REQUIRE_THAT(r.note, ContainsSubstring("not simple competitive"));
    }
}

TEST_CASE("Schild rejects data it cannot regress", "[pkpd]") {
    SECTION("two points leave no residual degrees of freedom") {
        const auto r = pkpd::schild({{1e-8, 2.0}, {1e-7, 11.0}});
        REQUIRE(r.pA2.provenance == Provenance::NotComputed);
        REQUIRE_FALSE(r.kbUsable);
        REQUIRE_THAT(r.note, ContainsSubstring("at least 3 points"));
    }

    SECTION("a dose ratio of 1 has no log10(DR - 1)") {
        const auto r = pkpd::schild({{1e-8, 1.0}, {1e-7, 11.0}, {1e-6, 101.0}});
        REQUIRE(r.pA2.provenance == Provenance::NotComputed);
        REQUIRE_THAT(r.note, ContainsSubstring("dose ratio"));
    }

    SECTION("a non-positive antagonist concentration has no logarithm") {
        const auto r = pkpd::schild({{0.0, 2.0}, {1e-7, 11.0}, {1e-6, 101.0}});
        REQUIRE(r.pA2.provenance == Provenance::NotComputed);
        REQUIRE_THAT(r.note, ContainsSubstring("antagonist concentration"));
    }
}
