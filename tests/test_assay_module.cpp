// Tests for the assay module adapter and the design simulator. Every case here
// defends a behaviour a user can see: the module dispatches a model to the engine
// that owns it, it refuses an experiment it does not have the metadata for, and the
// design simulator's confidence intervals cover the truth at the rate they claim.
//
// The coverage case is the important one. It is the only test in the suite that
// checks a reported uncertainty against reality rather than checking a number
// against another number, and it is the reason the design simulator reuses the
// production import/QC/fit path instead of a private one.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "assay/Design.h"
#include "contracts/Services.h"
#include "modules/AssayModule.h"

using namespace biocad;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

// A noiseless descending 4PL: signal 100 at zero concentration, plateau 0,
// midpoint 1e-7 mol/L, unit slope. y = 100 / (1 + x/1e-7).
std::string noiselessPlateCsv() {
    std::string s =
        "plate_id,well,role,sample_id,series_id,concentration,conc_unit,replicate,readout,"
        "readout_unit\n";
    const char* rows = "ABCDEFGH";
    for (int r = 0; r < 4; ++r) {
        s += std::string("P1,") + rows[r] + "1,negative,ctrl,,0,M," + std::to_string(r) +
             ",0,RFU\n";
        s += std::string("P1,") + rows[r] + "2,positive,ctrl,,0,M," + std::to_string(r) +
             ",100,RFU\n";
    }
    for (int i = 0; i < 10; ++i) {
        const double conc = 1.0e-5 / std::pow(10.0, i * 0.5);
        const double y = 100.0 / (1.0 + conc / 1.0e-7);
        char buf[192];
        std::snprintf(buf, sizeof buf, "P1,A%d,sample,cmpd1,s1,%.10g,M,0,%.10g,RFU\n", i + 3,
                      conc, y);
        s += buf;
    }
    return s;
}

AssayDesignSpec baseSpec() {
    AssayDesignSpec d;
    d.truthModel = AssayModel::FourParameterLogistic;
    // 5PL letter order: A (signal at zero), B (slope), C (midpoint), D (plateau).
    d.truthParameters = {100.0, 1.0, 1.0e-7, 0.0};
    for (int i = 0; i < 10; ++i) d.concentrations.push_back(1.0e-5 / std::pow(10.0, i * 0.5));
    d.replicates = 3;
    d.rows = 8;
    d.columns = 12;
    d.additiveNoiseSd = 2.0;
    d.seed = 20260813;
    d.replicateRuns = 200;
    return d;
}

}  // namespace

TEST_CASE("A long CSV round-trips import -> qc -> fit", "[assay][module]") {
    RealAssay assay;
    const auto ds = assay.import(noiselessPlateCsv());
    REQUIRE(ds);
    REQUIRE(ds->plates.size() == 1);
    REQUIRE(ds->detectedLayout == "long");

    const QcReport qc = assay.qc(ds->plates.front());
    // Zero-variance controls with a 100-unit window: perfect separation.
    CHECK_THAT(qc.zPrime.value, WithinAbs(1.0, 1e-12));
    CHECK(qc.interpretation.find("excellent") != std::string::npos);

    std::vector<Well> series;
    for (const auto& w : ds->plates.front().wells)
        if (w.role == WellRole::Sample) series.push_back(w);
    REQUIRE(series.size() == 10);

    const FitResult fit = assay.fit(series, AssayModel::FourParameterLogistic, false);
    REQUIRE(fit.converged);
    // The midpoint is recovered exactly from noiseless data, and it is a MODEL
    // number even though every well feeding it was measured.
    CHECK_THAT(std::log10(fit.derivedEc50.value), WithinAbs(-7.0, 1e-6));
    CHECK(fit.derivedEc50.provenance == Provenance::Model);
    CHECK_FALSE(fit.extrapolated);
    CHECK(fit.observations == 10);
}

TEST_CASE("Services::valid() requires the assay module", "[assay][module]") {
    RealAssay assay;
    Services s;
    // valid() only tests for null, so distinct non-null addresses exercise the exact
    // startup predicate without dragging in the Windows-only modules.
    char sentinel = 0;
    auto* p = &sentinel;
    s.library = reinterpret_cast<ILibrary*>(p);
    s.stability = reinterpret_cast<IStabilityModule*>(p);
    s.admet = reinterpret_cast<IAdmetModule*>(p);
    s.absorption = reinterpret_cast<IAbsorptionModule*>(p);
    s.similarity = reinterpret_cast<ISimilarityModule*>(p);
    s.legal = reinterpret_cast<ILegalModule*>(p);
    s.docking = reinterpret_cast<IDockingModule*>(p);
    s.runs = reinterpret_cast<IRunStore*>(p);
    s.pharmacodynamics = reinterpret_cast<IPharmacodynamicsModule*>(p);
    s.sequence = reinterpret_cast<ISequenceModule*>(p);
    s.structure = reinterpret_cast<IStructureModule*>(p);
    s.metabolismFacts = reinterpret_cast<IMetabolismFactsModule*>(p);
    s.alerts = reinterpret_cast<IAlertsModule*>(p);
    s.ionization = reinterpret_cast<IIonizationModule*>(p);
    s.nucleicAcid = reinterpret_cast<INucleicAcidModule*>(p);
    s.populationPk = reinterpret_cast<IPopulationPkModule*>(p);
    CHECK_FALSE(s.valid());
    s.assay = &assay;
    CHECK(s.valid());
}

TEST_CASE("An ITC fit from a bare well list is refused by name", "[assay][module]") {
    RealAssay assay;
    const auto ds = assay.import(noiselessPlateCsv());
    REQUIRE(ds);
    std::vector<Well> series;
    for (const auto& w : ds->plates.front().wells)
        if (w.role == WellRole::Sample) series.push_back(w);

    const FitResult f = assay.fit(series, AssayModel::WisemanIsotherm, false);
    CHECK_FALSE(f.converged);
    // The refusal must NAME the missing metadata. "Failed" would be indistinguishable
    // from a numerical problem, and inventing a cell volume would produce a
    // thermodynamic parameter set from a number nobody measured.
    CHECK(f.note.find("cell volume") != std::string::npos);
    CHECK(f.note.find("blank") != std::string::npos);
}

TEST_CASE("The same design seed reproduces the same run vector", "[assay][design]") {
    RealAssay assay;
    AssayDesignSpec spec = baseSpec();
    spec.replicateRuns = 40;
    const AssayDesignReport a = assay.simulate(spec);
    const AssayDesignReport b = assay.simulate(spec);
    REQUIRE(a.recoveredEc50.size() == b.recoveredEc50.size());
    REQUIRE_FALSE(a.recoveredEc50.empty());
    // Byte comparison, not tolerance: PCG64-DXSM plus a Box-Muller transform owned by
    // BioCAD exists precisely so this is exact.
    CHECK(std::memcmp(a.recoveredEc50.data(), b.recoveredEc50.data(),
                      a.recoveredEc50.size() * sizeof(double)) == 0);

    AssayDesignSpec other = spec;
    other.seed += 1;
    const AssayDesignReport c = assay.simulate(other);
    CHECK(c.recoveredEc50 != a.recoveredEc50);
}

TEST_CASE("A design's 95% intervals cover the truth about 95% of the time",
          "[assay][design]") {
    RealAssay assay;
    AssayDesignSpec spec = baseSpec();
    spec.replicateRuns = 1000;
    const AssayDesignReport r = assay.simulate(spec);
    REQUIRE(r.empiricalCoveragePct.provenance == Provenance::Model);
    REQUIRE(r.convergenceRatePct.value > 95.0);
    // Measured at 94.4% for this spec. The band is deliberately wide (+/- 5 points):
    // this asserts that the reported interval MEANS something, not that a Monte Carlo
    // estimate hit a particular value.
    CHECK(r.empiricalCoveragePct.value > 90.0);
    CHECK(r.empiricalCoveragePct.value < 100.0);
    // The recovered midpoints must be centred on the truth, or the interval is
    // covering the wrong thing for the right fraction of runs.
    CHECK_THAT(std::log10(r.medianEc50.value), WithinAbs(-7.0, 0.02));
}

TEST_CASE("The D-optimal ladder narrows the interval it was optimised for",
          "[assay][design]") {
    RealAssay assay;
    AssayDesignSpec spec = baseSpec();
    spec.replicateRuns = 300;

    const assay::DOptimalLadder opt = assay::dOptimalLadder(spec);
    REQUIRE(opt.concentrations.size() == spec.concentrations.size());
    CHECK(opt.dEfficiencyGain > 1.0);
    // Every chosen point must be a point the bench can actually reach.
    for (double c : opt.concentrations) {
        bool onLadder = false;
        for (double cand : opt.candidates)
            if (std::fabs(std::log10(c / cand)) < 1e-9) onLadder = true;
        CHECK(onLadder);
    }

    AssayDesignSpec dopt = spec;
    dopt.concentrations = opt.concentrations;
    const AssayDesignReport uniform = assay.simulate(spec);
    const AssayDesignReport optimal = assay.simulate(dopt);
    CHECK(optimal.ec50CiWidthLog10.value < uniform.ec50CiWidthLog10.value);
}

TEST_CASE("A plate whose controls overlap is reported unusable", "[assay][module]") {
    RealAssay assay;
    std::string csv =
        "plate_id,well,role,sample_id,series_id,concentration,conc_unit,replicate,readout,"
        "readout_unit\n";
    const char* rows = "ABCDEFGH";
    for (int r = 0; r < 6; ++r) {
        // Controls 15 units apart with a spread of 20: Z-prime is negative.
        csv += std::string("P9,") + rows[r] + "1,negative,ctrl,,0,M," + std::to_string(r) + "," +
               std::to_string(40 + 20 * (r % 2)) + ",RFU\n";
        csv += std::string("P9,") + rows[r] + "2,positive,ctrl,,0,M," + std::to_string(r) + "," +
               std::to_string(55 + 20 * (r % 2)) + ",RFU\n";
    }
    const auto ds = assay.import(csv);
    REQUIRE(ds);
    const QcReport qc = assay.qc(ds->plates.front());
    CHECK(qc.zPrime.value <= 0.0);
    CHECK(qc.interpretation.find("unusable") != std::string::npos);
}

TEST_CASE("The dilution calculator is exact and warns below the transfer floor",
          "[assay][design]") {
    // 250.3 g/mol, 10 mM, 1 mL: 250.3 * 0.01 * 0.001 = 2.503 mg.
    CHECK_THAT(1000.0 * assay::massForStock(250.3, 0.01, 0.001), WithinRel(2.503, 1e-12));

    const assay::DilutionPlan plan = assay::serialDilution(0.01, 1.0e-5, 10.0, 5, 100.0);
    REQUIRE(plan.steps.size() == 5);
    // The first well comes from the stock: 1000-fold, so 0.1 uL into 99.9 uL.
    CHECK_THAT(plan.steps.front().transferVolumeUl, WithinRel(0.1, 1e-12));
    CHECK_THAT(plan.steps.front().diluentVolumeUl, WithinRel(99.9, 1e-12));
    CHECK_THAT(plan.steps.back().concentration, WithinRel(1.0e-9, 1e-12));
    // 0.1 uL is below any hand pipette's honest floor, and the assumed CV no longer
    // holds there: the plan must say so rather than printing the volume.
    CHECK_FALSE(plan.warnings.empty());
    CHECK_THAT(plan.compoundedCvAtLastStep, WithinRel(std::sqrt(5.0), 1e-12));

    const assay::DilutionPlan bad = assay::serialDilution(1.0e-6, 1.0e-5, 3.0, 4, 100.0);
    CHECK(bad.steps.empty());
    CHECK_FALSE(bad.warnings.empty());
}

TEST_CASE("A trace-based truth model is refused, not simulated", "[assay][design]") {
    RealAssay assay;
    AssayDesignSpec spec = baseSpec();
    spec.truthModel = AssayModel::WisemanIsotherm;
    const AssayDesignReport r = assay.simulate(spec);
    CHECK(r.empiricalCoveragePct.provenance == Provenance::NotComputed);
    CHECK_FALSE(r.warnings.empty());
}
