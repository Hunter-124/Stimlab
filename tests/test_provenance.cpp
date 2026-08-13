// tests/test_provenance.cpp - the honesty core.
//
// These guard the one rule the whole application rests on: a derived number
// carries the reason it may be trusted, and a rank-ordering score may never
// wear physical units.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

#include "chem/AdmetModel.h"
#include "contracts/IDockingBackend.h"
#include "data/Domain.h"

using namespace biocad;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("makeQuantity rejects a heuristic carrying a physical unit", "[data][provenance]") {
    REQUIRE_THROWS_AS(makeQuantity(-7.4, "kcal/mol", 0.0, Provenance::Heuristic, "descriptor"),
                      std::invalid_argument);

    // Every other tier may carry a unit, and a unitless heuristic is fine.
    REQUIRE_NOTHROW(makeQuantity(-7.4, "kcal/mol", 2.85, Provenance::Model, "Vina 1.2.5"));
    REQUIRE_NOTHROW(makeQuantity(0.62, "", 0.0, Provenance::Heuristic, "descriptor"));

    const auto q = makeQuantity(-7.4, "kcal/mol", 2.85, Provenance::Predicted, "model X");
    REQUIRE(q.value == -7.4);
    REQUIRE(q.unit == "kcal/mol");
    REQUIRE(q.error == 2.85);
    REQUIRE(q.provenance == Provenance::Predicted);
}

TEST_CASE("notComputed names the missing prerequisite", "[data][provenance]") {
    const auto q = notComputed("Km");
    REQUIRE(q.provenance == Provenance::NotComputed);
    REQUIRE(q.source == "Km");
    REQUIRE(std::string(provenanceLabel(q.provenance)) == "not computed");
}

TEST_CASE("weakest() inherits the least trustworthy input tier", "[data][provenance]") {
    REQUIRE(weakest(Provenance::Measured, Provenance::Heuristic) == Provenance::Heuristic);
    REQUIRE(weakest(Provenance::Model, Provenance::Predicted) == Provenance::Model);
    REQUIRE(weakest(Provenance::Measured, Provenance::Measured) == Provenance::Measured);
    REQUIRE(weakest(Provenance::Heuristic, Provenance::NotComputed) == Provenance::NotComputed);
}

TEST_CASE("Quantity round-trips through JSON with its tier", "[data][provenance][json]") {
    const auto q = makeQuantity(1.25, "nM", 0.1, Provenance::Measured, "ChEMBL 42");
    const nlohmann::json j = q;
    REQUIRE(j.at("provenance") == "measured");
    const auto back = j.get<Quantity>();
    REQUIRE(back.provenance == Provenance::Measured);
    REQUIRE(back.source == "ChEMBL 42");
    REQUIRE(back.unit == "nM");
}

TEST_CASE("dG <-> Kd conversion is exact and reversible", "[chem][thermo]") {
    // dG = -9.0 kcal/mol at 298.15 K is Kd ~ 2.5e-7 M.
    const double kd = chem::kdFromDeltaG(-9.0);
    REQUIRE_THAT(kd, WithinRel(2.5e-7, 0.02));
    REQUIRE_THAT(chem::deltaGFromKd(kd), WithinAbs(-9.0, 1e-12));

    // 1 M is the standard state: dG = 0.
    REQUIRE_THAT(chem::deltaGFromKd(1.0), WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(chem::kdFromDeltaG(0.0), WithinAbs(1.0, 1e-15));

    // Temperature enters through RT, not as a fudge factor.
    REQUIRE(chem::kdFromDeltaG(-9.0, 310.15) > chem::kdFromDeltaG(-9.0, 298.15));
}

TEST_CASE("Ligand-efficiency metrics match hand computation", "[chem][efficiency]") {
    // LE = -dG / HAC: -(-9.0) / 20 = 0.45 kcal/mol per heavy atom.
    REQUIRE_THAT(chem::ligandEfficiency(-9.0, 20), WithinAbs(0.45, 1e-12));
    REQUIRE(chem::ligandEfficiency(-9.0, 0) == 0.0);  // degenerate input, not a division by zero

    // LLE = pIC50 - logP: 8.0 - 3.2 = 4.8.
    REQUIRE_THAT(chem::lipophilicEfficiency(8.0, 3.2), WithinAbs(4.8, 1e-12));

    // LELP = logP / LE: 3.2 / 0.45 = 7.111...
    REQUIRE_THAT(chem::leLipophilicityPrice(0.45, 3.2), WithinAbs(3.2 / 0.45, 1e-12));

    // BEI = pActivity / (MW/1000): 8.0 / 0.35 = 22.857...
    REQUIRE_THAT(chem::bindingEfficiencyIndex(8.0, 350.0), WithinAbs(8.0 / 0.35, 1e-12));

    // SEI = pActivity / (TPSA/100): 8.0 / 0.60 = 13.333...
    REQUIRE_THAT(chem::surfaceEfficiencyIndex(8.0, 60.0), WithinAbs(8.0 / 0.60, 1e-12));

    // RT.ln(10) at 300 K converts a pActivity to a binding free energy.
    REQUIRE_THAT(chem::deltaGFromPActivity(8.0), WithinAbs(-10.96, 1e-12));
}

TEST_CASE("A docking score is never converted into an affinity", "[docking][provenance]") {
    // A real engine result is Model: a constructed pose, not a measurement.
    DockJobResult real;
    real.engine = "AutoDock Vina 1.2.5";
    real.provenance = Provenance::Model;
    real.poses.push_back({1, -9.3, 0.0, 0.0, {}});
    REQUIRE(real.fromEngine());

    const auto q = makeQuantity(real.bestAffinity(), "kcal/mol", 0.0, real.provenance, real.engine);
    REQUIRE(q.provenance == Provenance::Model);
    REQUIRE(q.unit == "kcal/mol");

    // The descriptor fallback is Heuristic, so the unit must be empty - constructing
    // it with kcal/mol is a hard error, which is exactly what stops a fallback score
    // from being presented as an energy.
    DockJobResult fallback;
    fallback.engine = "descriptor-estimate";
    fallback.provenance = Provenance::Heuristic;
    fallback.poses.push_back({1, -6.1, 0.0, 0.0, {}});
    REQUIRE_FALSE(fallback.fromEngine());
    REQUIRE_THROWS_AS(makeQuantity(fallback.bestAffinity(), "kcal/mol", 0.0, fallback.provenance,
                                   fallback.engine),
                      std::invalid_argument);
}

TEST_CASE("Hepatic availability follows the well-stirred model", "[chem][pk]") {
    chem::PkLiabilities none;

    // With an explicitly supplied fu.CLint the arithmetic is checkable by hand:
    // Q_H = 90 L/h, fu.CLint = 30 L/h -> F_H = 90/120 = 0.75, CL_H = 2700/120 = 22.5 L/h.
    chem::HepaticAssumptions measured;
    measured.unboundIntrinsicClearanceLPerH = 30.0;
    measured.clIntMeasured = true;
    const auto p = chem::predictBioavailability(300.0, 2.0, 50.0, 1, none, measured);
    REQUIRE_THAT(p.firstPassSurvival, WithinAbs(0.75, 1e-12));
    REQUIRE_THAT(p.hepaticClearanceLPerH, WithinAbs(22.5, 1e-12));
    REQUIRE(p.clIntMeasured);

    // A hepatic extraction ratio can never exceed blood flow: CL_H < Q_H always.
    chem::HepaticAssumptions huge;
    huge.unboundIntrinsicClearanceLPerH = 1e9;
    const auto ext = chem::predictBioavailability(300.0, 2.0, 50.0, 1, none, huge);
    REQUIRE(ext.hepaticClearanceLPerH < huge.hepaticBloodFlowLPerH);
    REQUIRE(ext.firstPassSurvival > 0.0);

    // Without a supplied CLint the value is ASSUMED from liabilities and is therefore
    // rank-ordering only: a catechol must rank below an otherwise identical scaffold.
    chem::PkLiabilities catechol;
    catechol.catechol = true;
    const auto plain = chem::predictBioavailability(300.0, 2.0, 50.0, 1, none);
    const auto cat = chem::predictBioavailability(300.0, 2.0, 50.0, 1, catechol);
    REQUIRE_FALSE(plain.clIntMeasured);
    REQUIRE_FALSE(cat.clIntMeasured);
    REQUIRE(cat.bioavailabilityPct < plain.bioavailabilityPct);
    REQUIRE(cat.limitingRoute == "catechol COMT/MAO first-pass metabolism");

    // The rehabilitated model reproduces the old 1/(1+burden) shape exactly, so the
    // rescope changed the units and the honesty, not the ranking.
    REQUIRE_THAT(plain.firstPassSurvival,
                 WithinAbs(1.0 / (1.0 + plain.unboundIntrinsicClearanceLPerH / 90.0), 1e-12));
}
