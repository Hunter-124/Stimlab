// Tests for the Phase 11 ionization adapter (modules/IonizationModule) and its
// cited input pack.
//
// The chemistry itself is tested in test_chem_formula / test_chem_speciation /
// test_chem_solubility. What is under test HERE is the honesty boundary this
// adapter exists to enforce, and it is the more important half:
//
// (1) mass and isotope envelope are ALWAYS available, because they are arithmetic
//     on measured isotope masses and depend on nothing the pack supplies;
// (2) every pKa- or melting-point-dependent Quantity is NotComputed, with a
//     reason NAMING the missing input, for any compound the pack does not carry -
//     and the curve is absent rather than computed from a default;
// (3) a measured logP in the pack beats the Wildman-Crippen estimate, and the
//     resulting Quantity's provenance says which one it was;
// (4) a malformed pack fails LOUDLY. Every rejection path is exercised, because a
//     pack that silently half-loaded is indistinguishable from a compound that
//     genuinely has no cited pKa, and those call for opposite responses.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cstdlib>
#include <filesystem>
#include <string>

#include "modules/IonizationModule.h"

using namespace biocad;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;

namespace {

// The ionization pack is data on disk and ctest runs from the build directory, so
// point the loader at the in-tree copy the same way the Crippen tests do.
void useInTreePack() {
    const std::string dir = std::string(BIOCAD_ASSETS_DIR) + "/packs/descriptors";
#if defined(_WIN32)
    _putenv_s("BIOCAD_DESCRIPTOR_DIR", dir.c_str());
#else
    setenv("BIOCAD_DESCRIPTOR_DIR", dir.c_str(), 1);
#endif
}

IonizationPack inTreePack() {
    useInTreePack();
    return loadIonizationPack(std::filesystem::path(BIOCAD_ASSETS_DIR) / "packs" /
                              "descriptors" / "ionization.json");
}

Molecule compound(std::string id, std::string name, std::string smiles, std::string formula,
                  double mw) {
    Molecule m;
    m.id = std::move(id);
    m.name = std::move(name);
    m.smiles = std::move(smiles);
    m.formula = std::move(formula);
    m.molWeight = mw;
    return m;
}

// Ibuprofen: one cited carboxyl pKa, a melting point, and a measured logP.
Molecule ibuprofen() {
    return compound("ibuprofen", "Ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "C13H18O2", 206.28);
}

// Omeprazole: cited pyridinium 4.06 (basic) and benzimidazole N-H 8.80 (acidic),
// so it is a genuine ampholyte with an isoelectric point.
Molecule omeprazole() {
    return compound("omeprazole", "Omeprazole", "CC1=CN=C(C(=C1OC)C)CS(=O)c1[nH]c2cc(OC)ccc2n1",
                    "C17H19N3O3S", 345.42);
}

// Sertraline: a real compound with a real formula that this pack does not carry.
Molecule absentCompound() {
    return compound("sertraline", "Sertraline", "CNC1CCC(c2ccccc12)c1ccc(Cl)c(Cl)c1",
                    "C17H17Cl2N", 306.23);
}

}  // namespace

TEST_CASE("The shipped ionization pack loads with no errors", "[ionization][pack]") {
    const IonizationPack pack = inTreePack();
    INFO(pack.sourcePath);
    for (const auto& e : pack.errors) INFO(e);
    CHECK(pack.errors.empty());
    CHECK(pack.schemaVersion == kIonizationSchemaVersion);
    CHECK(pack.entries.size() >= 20);
    // The statement that these values are inputs rather than predictions is a
    // required field, not decoration: the panel and the tools both quote it.
    CHECK_THAT(pack.note, ContainsSubstring("INPUT"));
    CHECK(!pack.sources.empty());

    // Every cited group must carry a describable source and a Measured tier. A
    // Predicted pKa would need a named model and a benchmark error, and no such
    // model exists in this build.
    std::size_t groups = 0;
    for (const auto& [id, entry] : pack.entries) {
        for (const auto& g : entry.groups) {
            INFO(id << " / " << g.label);
            CHECK(g.pKa.provenance == Provenance::Measured);
            CHECK(!g.pKa.source.empty());
            CHECK(g.pKa.source.size() > 20);   // resolved prose, not a bare key
            CHECK(g.label.size() > 2);
            ++groups;
        }
        if (entry.hasMeltingPoint) CHECK(!entry.meltingPointSource.empty());
        if (entry.hasLogPMeasured) CHECK(!entry.logPSource.empty());
    }
    CHECK(groups >= 25);
}

TEST_CASE("A malformed ionization pack fails loudly", "[ionization][pack]") {
    SECTION("not JSON at all") {
        const auto p = parseIonizationPack("{not json", "<broken>");
        CHECK(p.entries.empty());
        REQUIRE(p.errors.size() == 1);
        CHECK_THAT(p.errors[0], ContainsSubstring("not valid JSON"));
    }
    SECTION("an unknown schemaVersion is refused, not skipped") {
        const auto p = parseIonizationPack(R"({"schemaVersion":99,"note":"n"})", "<v99>");
        CHECK(p.entries.empty());
        REQUIRE(p.errors.size() == 1);
        CHECK_THAT(p.errors[0], ContainsSubstring("schemaVersion 99"));
    }
    SECTION("the inputs-not-predictions note is mandatory") {
        const auto p = parseIonizationPack(R"({"schemaVersion":1,"compounds":[]})", "<no-note>");
        REQUIRE(p.errors.size() == 1);
        CHECK_THAT(p.errors[0], ContainsSubstring("note"));
    }
    SECTION("an uncited pKa is rejected") {
        const auto p = parseIonizationPack(
            R"({"schemaVersion":1,"note":"n","compounds":[
                 {"id":"x","groups":[{"label":"carboxyl","pKa":4.2,"acidic":true}]}]})",
            "<uncited>");
        CHECK(p.entries.empty());
        REQUIRE(p.errors.size() == 1);
        CHECK_THAT(p.errors[0], ContainsSubstring("uncited pKa is not an input"));
    }
    SECTION("a group with no acidic/basic direction is rejected") {
        const auto p = parseIonizationPack(
            R"({"schemaVersion":1,"note":"n","compounds":[
                 {"id":"x","groups":[{"label":"carboxyl","pKa":4.2,"source":"a paper"}]}]})",
            "<no-direction>");
        CHECK(p.entries.empty());
        REQUIRE(p.errors.size() == 1);
        CHECK_THAT(p.errors[0], ContainsSubstring("acidic"));
    }
    SECTION("a melting point with no source is rejected") {
        const auto p = parseIonizationPack(
            R"({"schemaVersion":1,"note":"n","compounds":[{"id":"x","meltingPointC":150.0}]})",
            "<uncited-mp>");
        CHECK(p.entries.empty());
        REQUIRE(p.errors.size() == 1);
        CHECK_THAT(p.errors[0], ContainsSubstring("meltingPointSource"));
    }
    SECTION("a duplicate id is an error, not an arbitrary winner") {
        const auto p = parseIonizationPack(
            R"({"schemaVersion":1,"note":"n","compounds":[
                 {"id":"x","meltingPointC":150.0,"meltingPointSource":"a monograph"},
                 {"id":"x","meltingPointC":160.0,"meltingPointSource":"a monograph"}]})",
            "<dupe>");
        CHECK(p.entries.size() == 1);
        REQUIRE(p.errors.size() == 1);
        CHECK_THAT(p.errors[0], ContainsSubstring("duplicate id"));
    }
    SECTION("an entry that supplies no input at all is rejected") {
        const auto p = parseIonizationPack(
            R"({"schemaVersion":1,"note":"n","compounds":[{"id":"x"}]})", "<empty>");
        CHECK(p.entries.empty());
        REQUIRE(p.errors.size() == 1);
        CHECK_THAT(p.errors[0], ContainsSubstring("no input at all"));
    }
    SECTION("an out-of-range pKa is rejected") {
        const auto p = parseIonizationPack(
            R"({"schemaVersion":1,"note":"n","compounds":[
                 {"id":"x","groups":[{"label":"carboxyl","pKa":42.0,"acidic":true,
                                      "source":"a paper"}]}]})",
            "<absurd>");
        CHECK(p.entries.empty());
        REQUIRE(p.errors.size() == 1);
        CHECK_THAT(p.errors[0], ContainsSubstring("outside the range"));
    }
    SECTION("a missing pack file is an error, never an empty success") {
        const auto p = loadIonizationPack("/nonexistent/ionization.json");
        CHECK(p.entries.empty());
        REQUIRE(p.errors.size() == 1);
        CHECK_THAT(p.errors[0], ContainsSubstring("not found"));
    }
}

TEST_CASE("Mass and isotope envelope never depend on the pack", "[ionization]") {
    useInTreePack();
    // Constructed with an EMPTY pack on purpose: if these numbers survive that,
    // they genuinely do not depend on a cited pKa.
    const RealIonization emptyPack{IonizationPack{}};

    for (const Molecule& m : {ibuprofen(), omeprazole(), absentCompound()}) {
        const IonizationReport r = emptyPack.analyze(m);
        INFO(m.id);
        CHECK(r.mass.monoisotopic.provenance == Provenance::Measured);
        CHECK(r.mass.average.provenance == Provenance::Measured);
        CHECK(!r.mass.elements.empty());
        CHECK(r.mass.average.value == Catch::Approx(m.molWeight).margin(0.05));
        CHECK(r.envelope.peaks.size() >= 2);
        CHECK(!r.envelope.source.empty());
        // But everything the pack would have supplied is absent, by name.
        CHECK(r.speciation.points.empty());
        CHECK(r.speciation.isoelectricPoint.provenance == Provenance::NotComputed);
        CHECK_THAT(r.speciation.isoelectricPoint.source, ContainsSubstring("pKa"));
    }
}

TEST_CASE("A compound absent from the pack is NotComputed by name, never defaulted",
          "[ionization]") {
    useInTreePack();
    const RealIonization ion{inTreePack()};
    const IonizationReport r = ion.analyze(absentCompound());

    // Composition is still fully available - that is the invariant.
    REQUIRE(r.mass.monoisotopic.provenance == Provenance::Measured);
    CHECK(r.mass.formula == "C17H17Cl2N");
    REQUIRE(r.envelope.peaks.size() >= 3);
    // Two chlorines make M+2 large; a wrong isotope table would not produce this.
    bool bigM2 = false;
    for (const auto& p : r.envelope.peaks)
        if (p.nominalShift == 2 && p.intensity > 0.5) bigM2 = true;
    CHECK(bigM2);

    // Every pack-dependent Quantity is NotComputed AND names its prerequisite.
    const Quantity dependent[] = {
        r.speciation.isoelectricPoint, r.speciation.logDAtPh74,
        r.solubility.intrinsic,        r.solubility.pHmax,
        r.solubility.solubilityAtPh74, r.solubility.doseNumber,
        r.solubility.dissolutionNumber, r.solubility.absorptionNumber,
        r.buffer.betaAtPh74,           r.buffer.maxCapacity,
        r.buffer.maxCapacityPh,        r.dissolution.timeTo85Pct};
    for (const Quantity& q : dependent) {
        CHECK(q.provenance == Provenance::NotComputed);
        CHECK(!q.source.empty());
    }
    CHECK_THAT(r.solubility.intrinsic.source, ContainsSubstring("melting point"));
    // No curve is fabricated from a plausible default.
    CHECK(r.speciation.points.empty());
    CHECK(r.solubility.curve.empty());
    CHECK(r.buffer.curve.empty());
    CHECK(r.dissolution.points.empty());

    // logP is NOT pack-dependent: Wildman-Crippen computes it from structure, so
    // it must still be there, tiered Predicted and sourced to the method.
    CHECK(r.speciation.logP.provenance == Provenance::Predicted);
    CHECK_THAT(r.speciation.logP.source, ContainsSubstring("Crippen"));
}

TEST_CASE("A pack-measured logP beats the Wildman-Crippen estimate", "[ionization]") {
    useInTreePack();
    const RealIonization ion{inTreePack()};

    const Quantity measured = ion.logPFor(ibuprofen());
    CHECK(measured.provenance == Provenance::Measured);
    CHECK_THAT(measured.value, WithinAbs(3.97, 1e-12));

    // The same molecule with no pack entry falls back to Crippen, which must
    // disagree - otherwise this test would pass for the wrong reason.
    const RealIonization noPack{IonizationPack{}};
    const Quantity estimated = noPack.logPFor(ibuprofen());
    REQUIRE(estimated.provenance == Provenance::Predicted);
    CHECK_THAT(estimated.source, ContainsSubstring("Crippen"));
    // The method's own published RMS error travels with the number.
    CHECK(estimated.error > 0.0);
    CHECK(estimated.value != Catch::Approx(measured.value));
}

TEST_CASE("An ampholyte in the pack gets a real isoelectric point", "[ionization]") {
    useInTreePack();
    const RealIonization ion{inTreePack()};
    const IonizationReport r = ion.analyze(omeprazole());

    // Two independent groups -> four microstates, and a genuine pI because the net
    // charge runs +1 -> 0 -> -1.
    REQUIRE(r.speciation.labels.size() == 4);
    REQUIRE(r.speciation.points.size() > 100);
    REQUIRE(r.speciation.isoelectricPoint.provenance != Provenance::NotComputed);
    CHECK_THAT(r.speciation.isoelectricPoint.value, WithinAbs(0.5 * (4.06 + 8.80), 0.06));

    // Every pH point's microspecies fractions must sum to one: the ladder is a
    // distribution, and a distribution that does not close is not one.
    for (const auto& p : r.speciation.points) {
        double sum = 0.0;
        for (double f : p.microspeciesFractions) sum += f;
        INFO("pH " << p.pH);
        CHECK_THAT(sum, WithinAbs(1.0, 1e-10));
    }

    // The GSE ran because the pack carries a melting point.
    CHECK(r.solubility.intrinsic.provenance != Provenance::NotComputed);
    CHECK(!r.solubility.curve.empty());
    CHECK(!r.buffer.curve.empty());
    // A polyprotic compound must SAY that its monoprotic solubility profile used
    // one group, rather than presenting the curve as the whole story.
    bool saidSo = false;
    for (const auto& w : r.solubility.warnings)
        if (w.find("monoprotic model") != std::string::npos) saidSo = true;
    CHECK(saidSo);
}

TEST_CASE("A monoprotic acid has no isoelectric point and says so", "[ionization]") {
    useInTreePack();
    const RealIonization ion{inTreePack()};
    const IonizationReport r = ion.analyze(ibuprofen());

    // Net charge runs 0 -> -1 and never crosses zero, so NotComputed here is
    // correct chemistry rather than a missing input - and the reason must say
    // that, not name a prerequisite the user could go and supply.
    CHECK(r.speciation.isoelectricPoint.provenance == Provenance::NotComputed);
    CHECK(!r.speciation.isoelectricPoint.source.empty());
    // logD, which needs only a pKa and a logP, IS computed.
    REQUIRE(r.speciation.logDAtPh74.provenance != Provenance::NotComputed);
    // Ibuprofen at pH 7.4 is essentially fully ionized, so logD sits far below
    // logP 3.97.
    CHECK(r.speciation.logDAtPh74.value < 2.0);
    CHECK(r.solubility.intrinsic.provenance != Provenance::NotComputed);
}

TEST_CASE("Dissolution needs formulation inputs the pack cannot supply", "[ionization]") {
    useInTreePack();
    const RealIonization ion{inTreePack()};
    const IonizationReport r = ion.analyze(ibuprofen());

    // A dose and a particle size distribution are properties of a formulation, not
    // of a molecule, so even a fully-covered compound gets NotComputed here - and
    // the reason names all four rather than assuming any of them.
    REQUIRE(r.dissolution.timeTo85Pct.provenance == Provenance::NotComputed);
    CHECK_THAT(r.dissolution.timeTo85Pct.source, ContainsSubstring("dose"));
    CHECK_THAT(r.dissolution.timeTo85Pct.source, ContainsSubstring("particle radius"));
    CHECK(r.dissolution.points.empty());
    CHECK(!r.dissolution.assumptions.empty());
}

TEST_CASE("The pack never carries a Predicted pKa", "[ionization][pack]") {
    // The structural guarantee: there is no code path that writes a Predicted or
    // Heuristic pKa, so no shipped or hand-written pack can introduce one. A
    // Predicted pKa would require a named model and a benchmark error, and this
    // build contains no pKa model at all.
    const auto p = parseIonizationPack(
        R"({"schemaVersion":1,"note":"INPUTS","compounds":[
             {"id":"x","groups":[{"label":"carboxyl","pKa":4.2,"acidic":true,
                                  "source":"a potentiometric determination in the literature"}]}]})",
        "<one>");
    REQUIRE(p.errors.empty());
    REQUIRE(p.entries.count("x") == 1);
    const auto& g = p.entries.at("x").groups.at(0);
    CHECK(g.pKa.provenance == Provenance::Measured);
    CHECK(g.pKa.unit.empty());   // a pKa is a log10 ratio and carries no unit
}

TEST_CASE("Pass-through entry points reach the real chemistry", "[ionization]") {
    useInTreePack();
    const RealIonization ion{inTreePack()};

    SECTION("formula parses and refuses") {
        const auto fm = ion.formula("C9H13N");
        REQUIRE(fm.has_value());
        CHECK(fm->formula == "C9H13N");
        CHECK_THAT(fm->monoisotopic.value, WithinAbs(135.104799, 1e-4));
        // A neutral has no m/z, and that is stated rather than left blank.
        CHECK(fm->mz.provenance == Provenance::NotComputed);
        CHECK(!ion.formula("this is not a formula").has_value());
    }
    SECTION("an unparsable envelope says why instead of returning an empty plot") {
        const IsotopeEnvelope env = ion.envelope("!!!", 1e-4);
        CHECK(env.peaks.empty());
        CHECK_THAT(env.source, ContainsSubstring("not computed"));
    }
    SECTION("balance reaches the integer null space") {
        const BalancedEquation b = ion.balance({"C3H8", "O2"}, {"CO2", "H2O"}, {});
        REQUIRE(b.balanced);
        REQUIRE(b.reactantCoefficients.size() == 2);
        REQUIRE(b.productCoefficients.size() == 2);
        CHECK(b.reactantCoefficients[0] == 1);
        CHECK(b.reactantCoefficients[1] == 5);
        CHECK(b.productCoefficients[0] == 3);
        CHECK(b.productCoefficients[1] == 4);
    }
    SECTION("solvePh reaches the tableau solver and converges") {
        // 0.1 M acetic acid, the phase's headline fixture: pH 2.881, residuals
        // below 1e-10.
        SpeciationProblem p;
        p.components = {"H", "A"};
        p.totals = {0.1, 0.1};
        p.species = {"HA", "A-"};
        p.stoichiometry = {{1, 1}, {0, 1}};
        p.logK = {4.756, 0.0};
        p.charges = {0.0, -1.0};
        const SpeciationResult r = ion.solvePh(p);
        REQUIRE(r.converged);
        CHECK(r.massBalanceResidual < 1e-10);
        CHECK_THAT(r.pH, WithinAbs(2.881, 0.002));
    }
    SECTION("buffer reaches Van Slyke") {
        // 0.1 M buffer at pH == pKa: beta = 0.0576 M/pH.
        const BufferReport b = ion.buffer({BufferComponent{"acetate", 4.756, 0.1}});
        REQUIRE(!b.curve.empty());
        double atPka = 0.0;
        for (const auto& p : b.curve)
            if (std::abs(p.pH - 4.756) < 0.011) atPka = std::max(atPka, p.beta);
        CHECK_THAT(atPka, WithinAbs(0.0576, 5e-4));
    }
    SECTION("solubility refuses without a melting point rather than defaulting one") {
        const Quantity logP =
            makeQuantity(3.97, "", 0.0, Provenance::Measured, "a shake-flask determination");
        const std::vector<IonizableGroup> groups = {
            {"carboxyl", makeQuantity(4.91, "", 0.0, Provenance::Measured, "a pKa table"), true}};
        const SolubilityReport withMp = ion.solubility(ibuprofen(), groups, logP, 76.0);
        CHECK(withMp.intrinsic.provenance != Provenance::NotComputed);
        // meltingPointC <= 0 is the contract's "not supplied".
        const SolubilityReport noMp = ion.solubility(ibuprofen(), groups, logP, 0.0);
        CHECK(noMp.intrinsic.provenance == Provenance::NotComputed);
        CHECK(!noMp.intrinsic.source.empty());
    }
}
