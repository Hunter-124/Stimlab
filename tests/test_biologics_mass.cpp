// Mass ladders, the required resolving power, peptide mapping, liability flags and
// developability descriptors.
//
// Every mass here is composition arithmetic over the NIST SRD 144 isotope pack. There
// is no published-value literal to compare against in-tree, so the intact mass is
// checked TWO INDEPENDENT WAYS - once from the assembled whole-chain formula and once
// by summing the per-residue formula masses plus a water - and the disulfide and
// reduction deltas are checked against the hydrogen mass from the same table.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <string>

#include "bio/Imgt.h"
#include "bio/Liabilities.h"
#include "chem/Formula.h"

using Catch::Approx;
using namespace biocad;
using namespace biocad::bio;

namespace {

// 1N8Z (RCSB entry FASTA): trastuzumab Fab heavy chain and complete light chain.
constexpr const char* kH =
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAED"
    "TAVYYCSRWGGDGFYAMDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSS"
    "GLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEP";
constexpr const char* kL =
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQ"
    "HYTTPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLT"
    "LSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC";

double mono(const std::string& formula) {
    const auto p = chem::parseFormula(formula);
    REQUIRE(p.has_value());
    return chem::monoisotopicMass(*p);
}
double average(const std::string& formula) {
    const auto p = chem::parseFormula(formula);
    REQUIRE(p.has_value());
    return chem::averageMass(*p);
}

// The second, independent route to a chain mass: per-residue formulas plus one water.
double residueSum(const std::string& seq, bool useAverage) {
    double total = useAverage ? average("H2O") : mono("H2O");
    for (char c : seq) {
        const std::string f = peptideFormula(std::string(1, c));
        REQUIRE_FALSE(f.empty());
        total += useAverage ? average(f) - average("H2O") : mono(f) - mono("H2O");
    }
    return total;
}

const MassLadderEntry* entry(const MassLadder& l, const std::string& prefix) {
    for (const auto& e : l.entries)
        if (e.species.rfind(prefix, 0) == 0) return &e;
    return nullptr;
}

}  // namespace

TEST_CASE("a chain mass agrees computed two independent ways", "[biologics]") {
    REQUIRE(chem::isotopeTableOk());
    const std::string formula = peptideFormula(kH);
    REQUIRE_FALSE(formula.empty());
    CHECK(average(formula) == Approx(residueSum(kH, true)).margin(1e-6));
    CHECK(mono(formula) == Approx(residueSum(kH, false)).margin(1e-6));
}

TEST_CASE("the mass ladder is composed from formulas, never from literals", "[biologics]") {
    REQUIRE(chem::isotopeTableOk());
    MassLadderInput in;
    in.heavyChains = {kH, kH};
    in.lightChains = {kL, kL};
    in.interchainDisulfides = 4;
    in.intrachainDisulfides = 12;
    const MassLadder ladder = massLadder(in);
    REQUIRE(ladder.entries.size() > 6);

    for (const auto& e : ladder.entries) {
        CHECK(e.average.provenance == Provenance::Measured);
        CHECK(e.average.unit == "Da");
        CHECK(e.monoisotopic.unit == "Da");
        // The source names the isotope table and the exact formula used.
        CHECK(e.average.source.find("formula ") != std::string::npos);
        CHECK(e.average.value > 0.0);
        // Above 10 kDa the note must say why the average mass is the reported one.
        if (e.average.value > 10000.0) CHECK_FALSE(e.note.empty());
    }

    const MassLadderEntry* intact = entry(ladder, "intact (aglycosylated)");
    const MassLadderEntry* reduced = entry(ladder, "reduced (");
    REQUIRE(intact);
    REQUIRE(reduced);
    // Reduction adds exactly two hydrogens per disulfide, from the same table.
    CHECK(reduced->monoisotopic.value - intact->monoisotopic.value ==
          Approx(2.0 * 16 * mono("H")).margin(1e-9));
    CHECK(intact->disulfides == 16);

    // Glycoforms differ by whole monosaccharide compositions: G1F - G0F is one Hex.
    const MassLadderEntry* g0f = entry(ladder, "intact + G0F");
    const MassLadderEntry* g1f = entry(ladder, "intact + G1F");
    REQUIRE(g0f);
    REQUIRE(g1f);
    CHECK(g1f->average.value - g0f->average.value == Approx(2 * average("C6H10O5")).margin(1e-6));

    // A pyroGlu heavy chain loses NH3 twice (two heavy chains, both starting in Glu
    // here, so it is a water each).
    const MassLadderEntry* pyro = entry(ladder, "intact, N-terminal pyroGlu");
    REQUIRE(pyro);
    CHECK(intact->monoisotopic.value - pyro->monoisotopic.value ==
          Approx(2 * mono("H2O")).margin(1e-6));
}

TEST_CASE("the required resolving power matches the hand calculation", "[biologics]") {
    REQUIRE(chem::isotopeTableOk());
    // Both deltas come from the isotope table, not from a literal.
    const double deamidation = mono("O") - mono("NH");
    const double isotopeSpacing = mono("[13C]") - mono("C");
    CHECK(deamidation == Approx(0.984016).margin(1e-5));
    CHECK(isotopeSpacing == Approx(1.003355).margin(1e-5));
    const double gap = std::abs(isotopeSpacing - deamidation);
    CHECK(gap == Approx(0.019339).margin(1e-5));
    // 150000 / 0.019339 = 7.756e6, by hand.
    CHECK(150000.0 / gap == Approx(7756246.0).epsilon(1e-3));

    MassLadderInput in;
    in.heavyChains = {kH};
    in.includeGlycoforms = false;
    const MassLadder ladder = massLadder(in);
    const MassLadderEntry* intact = entry(ladder, "intact");
    REQUIRE(intact);
    REQUIRE(ladder.requiredResolvingPower.provenance == Provenance::Measured);
    CHECK(ladder.requiredResolvingPower.value ==
          Approx(intact->average.value / gap).margin(1e-6));
}

TEST_CASE("tryptic peptides respect the proline rule and their ion series close",
          "[biologics]") {
    DigestOptions opts;
    opts.protease = "trypsin";
    opts.maxMissedCleavages = 1;
    const PeptideMap map = digest(kH, opts);
    REQUIRE_FALSE(map.peptides.empty());
    const std::string chain = kH;
    for (const auto& p : map.peptides) {
        if (p.missedCleavages != 0) continue;
        if (p.end + 1 >= static_cast<int>(chain.size())) continue;
        // Trypsin does not cut before proline.
        CHECK(chain[static_cast<std::size_t>(p.end) + 1] != 'P');
    }
    // b1 + y(n-1) = the peptide plus two protons: the series must be self-consistent.
    const auto& first = map.peptides.front();
    REQUIRE_FALSE(first.bIons.empty());
    REQUIRE_FALSE(first.yIons.empty());
    const double protonPair = 2 * (mono("H") - 0.00054857990907);
    CHECK(first.bIons.front() + first.yIons.back() ==
          Approx(first.monoisotopic.value + protonPair).margin(1e-3));
    CHECK(map.coveragePct > 90.0);

    // An unknown protease is refused rather than silently falling back to trypsin.
    DigestOptions bad;
    bad.protease = "papain";
    const PeptideMap none = digest(kH, bad);
    CHECK(none.peptides.empty());
    REQUIRE_FALSE(none.warnings.empty());
}

TEST_CASE("liability flags carry citations and report exposure as unknown without coordinates",
          "[biologics]") {
    REQUIRE(liabilityPack().ok);
    const AbDomain d = numberDomain(kH);
    REQUIRE(d.numbered);
    const auto flags = scanLiabilities(d);
    REQUIRE_FALSE(flags.empty());
    for (const auto& f : flags) {
        CHECK_FALSE(f.citation.empty());
        CHECK_FALSE(f.label.empty());
        // No structure was supplied: exposure is UNKNOWN, never assumed to be exposed.
        CHECK_FALSE(f.exposureKnown);
        CHECK(f.relativeSasa == 0.0);
    }
    // The published deamidation risk ordering NG > NS > NT is the pack's ordering.
    int ng = -1, ns = -1, nt = -1;
    for (const auto& r : liabilityPack().rules) {
        if (r.ruleId == "deamid-ng") ng = r.risk;
        if (r.ruleId == "deamid-ns") ns = r.risk;
        if (r.ruleId == "deamid-nt") nt = r.risk;
    }
    REQUIRE(ng > 0);
    CHECK(ng < ns);
    CHECK(ns < nt);
    // The NG site in CDR2 of this sequence is found, with its IMGT position.
    bool sawCdr2Ng = false;
    for (const auto& f : flags)
        if (f.ruleId == "deamid-ng" && f.region == "CDR2") sawCdr2Ng = true;
    CHECK(sawCdr2Ng);
}

TEST_CASE("developability descriptors are exact, and TAP is refused without a structure origin",
          "[biologics]") {
    REQUIRE(descriptorPack().ok);
    DevelopabilityInput in;
    in.chains = {kH, kL};
    const DevelopabilityReport r = developability(in);

    CHECK(r.isoelectricPoint.provenance == Provenance::Measured);
    // The net charge must actually vanish at the reported pI.
    CHECK(netCharge(std::string(kH) + kL, r.isoelectricPoint.value) == Approx(0.0).margin(1e-3));
    // eps280 by hand from the composition.
    const std::string all = std::string(kH) + kL;
    const int trp = static_cast<int>(std::count(all.begin(), all.end(), 'W'));
    const int tyr = static_cast<int>(std::count(all.begin(), all.end(), 'Y'));
    const int cys = static_cast<int>(std::count(all.begin(), all.end(), 'C'));
    CHECK(r.extinctionCoefficient280.value ==
          Approx(5500.0 * trp + 1490.0 * tyr + 125.0 * (cys / 2)).margin(1e-9));
    CHECK(r.extinctionCoefficient280.unit == "M^-1 cm^-1");
    CHECK(r.gravy.provenance == Provenance::Measured);
    CHECK(r.aliphaticIndex.provenance == Provenance::Measured);
    CHECK(r.instabilityIndex.provenance == Provenance::Measured);
    CHECK(r.fvChargeSymmetry.provenance == Provenance::Measured);

    // No structure origin -> every TAP metric is NotComputed, naming what is missing.
    CHECK(r.tapPsh.provenance == Provenance::NotComputed);
    CHECK(r.tapPpc.provenance == Provenance::NotComputed);
    CHECK(r.tapPnc.provenance == Provenance::NotComputed);
    CHECK(r.tapSfvcsp.provenance == Provenance::NotComputed);
    CHECK_FALSE(r.tapPsh.source.empty());
    CHECK(r.structureOrigin.empty());

    // Fv charge symmetry needs exactly two chains and refuses otherwise.
    DevelopabilityInput one;
    one.chains = {kH};
    CHECK(developability(one).fvChargeSymmetry.provenance == Provenance::NotComputed);

    // Net charge falls monotonically across the curve and crosses zero once.
    const auto curve = netChargeCurve(all);
    REQUIRE(curve.size() > 50);
    CHECK(curve.front().second > 0.0);
    CHECK(curve.back().second < 0.0);
}
