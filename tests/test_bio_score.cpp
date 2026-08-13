// Tests for the two exact geometry scores. Each case is one that would actually catch a
// bug: lDDT's defining identity and its rigid-motion invariance, the locality of its
// per-residue breakdown, SASA against a closed-form single sphere, and the reproducibility
// contract carried in every source string.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <string>

#include "bio/Score.h"

using namespace biocad;
using namespace biocad::bio;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;

namespace {

constexpr double kPi = 3.14159265358979323846;

Atom mkAtom(const char* name, const char* element, double x, double y, double z) {
    Atom a;
    a.name = name;
    a.element = element;
    a.x = x;
    a.y = y;
    a.z = z;
    return a;
}

// A crude extended chain: consecutive CA 3.8 A apart, plus two more atoms per residue so
// there are inter-residue distances to score.
Structure makeChain(int residues) {
    Structure s;
    s.id = "TEST";
    Model m;
    Chain c;
    c.id = "A";
    for (int i = 0; i < residues; ++i) {
        Residue r;
        r.name = (i % 2) ? "ALA" : "SER";
        r.authSeqId = i + 1;
        const double x = 3.8 * i;
        const double y = 1.2 * std::sin(0.6 * i);
        const double z = 1.2 * std::cos(0.6 * i);
        r.atoms.push_back(mkAtom(" CA ", "C", x, y, z));
        r.atoms.push_back(mkAtom(" CB ", "C", x + 1.5, y + 1.4, z));
        r.atoms.push_back(mkAtom(" O  ", "O", x + 0.6, y - 1.2, z + 0.4));
        c.residues.push_back(r);
    }
    m.chains.push_back(c);
    s.models.push_back(m);
    return s;
}

// A genuine 3D rigid motion: rotate about z, then about x, then translate.
void rigidMove(Structure& s, double angle, double tx, double ty, double tz) {
    const double ca = std::cos(angle);
    const double sa = std::sin(angle);
    for (auto& m : s.models) {
        for (auto& ch : m.chains) {
            for (auto& r : ch.residues) {
                for (auto& a : r.atoms) {
                    const double x1 = ca * a.x - sa * a.y;
                    const double y1 = sa * a.x + ca * a.y;
                    const double z1 = a.z;
                    a.x = x1 + tx;
                    a.y = ca * y1 - sa * z1 + ty;
                    a.z = sa * y1 + ca * z1 + tz;
                }
            }
        }
    }
}

Structure singleAtom(const char* element) {
    Structure s;
    Model m;
    Chain c;
    c.id = "A";
    Residue r;
    r.name = "LIG";
    r.authSeqId = 1;
    r.atoms.push_back(mkAtom(" C1 ", element, 0.0, 0.0, 0.0));
    c.residues.push_back(r);
    m.chains.push_back(c);
    s.models.push_back(m);
    return s;
}

}  // namespace

TEST_CASE("lDDT of a structure against itself is exactly 1.0", "[bio]") {
    const Structure s = makeChain(20);
    const LddtResult r = lddt(s, s);

    // Exact equality is the point: every distance difference is identically zero, so all
    // four tolerance tests pass and the mean of four ones is one, never 0.999...
    REQUIRE(r.global.value == 1.0);
    REQUIRE(r.global.provenance == Provenance::Measured);
    REQUIRE(r.global.unit.empty());
    REQUIRE(r.consideredPairs > 0);
    REQUIRE(r.perResidue.size() == 20);
    REQUIRE(r.unmatchedModelResidues == 0);
    REQUIRE(r.unmatchedReferenceResidues == 0);
}

TEST_CASE("lDDT is invariant under a rigid rotation and translation", "[bio]") {
    const Structure ref = makeChain(20);

    Structure moved = ref;
    rigidMove(moved, 0.6457718232379019 /* 37 degrees */, 13.0, -4.0, 2.5);
    // Superposition-free means a rigid motion cannot change the score at all.
    REQUIRE(lddt(moved, ref).global.value == 1.0);

    Structure perturbed = ref;
    perturbed.models[0].chains[0].residues[10].atoms[0].x += 0.6;
    const double before = lddt(perturbed, ref).global.value;
    rigidMove(perturbed, 0.6457718232379019, 13.0, -4.0, 2.5);
    REQUIRE_THAT(lddt(perturbed, ref).global.value, WithinAbs(before, 1e-12));
}

TEST_CASE("lDDT degrades locally with the size of a single-atom perturbation", "[bio]") {
    const Structure ref = makeChain(20);

    // 0.3 A moves every distance by at most 0.3 A, below all four tolerances (0.5/1/2/4),
    // so the score is still exactly 1. That is arithmetic, not insensitivity.
    Structure tiny = ref;
    tiny.models[0].chains[0].residues[10].atoms[0].x += 0.3;
    REQUIRE(lddt(tiny, ref).global.value == 1.0);

    Structure small = ref;
    small.models[0].chains[0].residues[10].atoms[0].x += 0.6;
    const LddtResult ls = lddt(small, ref);
    REQUIRE(ls.global.value < 1.0);
    REQUIRE(ls.global.value > 0.95);

    Structure big = ref;
    big.models[0].chains[0].residues[10].atoms[0].x += 10.0;
    const LddtResult lb = lddt(big, ref);
    REQUIRE(lb.global.value < ls.global.value);
    REQUIRE(lb.perResidue[10].score < ls.perResidue[10].score);
    // The damage is attributed to residue 11 specifically; a distant residue is untouched.
    REQUIRE(lb.perResidue[0].score == 1.0);
    REQUIRE(lb.perResidue[10].score < 0.8);
    REQUIRE(lb.perResidue[10].key.authSeqId == 11);
}

TEST_CASE("lDDT refuses a comparison it cannot make", "[bio]") {
    const Structure ref = makeChain(20);

    const Structure empty;
    REQUIRE(lddt(empty, ref).global.provenance == Provenance::NotComputed);
    REQUIRE_THAT(lddt(empty, ref).global.source, ContainsSubstring("atoms"));

    const Structure shorter = makeChain(19);
    const LddtResult mismatch = lddt(shorter, ref);
    REQUIRE(mismatch.global.provenance == Provenance::NotComputed);
    REQUIRE_THAT(mismatch.global.source, ContainsSubstring("57"));
    REQUIRE_THAT(mismatch.global.source, ContainsSubstring("60"));
}

TEST_CASE("pairResidues reports what it could not match", "[bio]") {
    const Structure ref = makeChain(20);
    const Structure shorter = makeChain(19);

    const ResiduePairing p = pairResidues(*shorter.model(), *ref.model());
    REQUIRE(p.pairs.size() == 19);
    REQUIRE(p.unmatchedReference == 1);
    REQUIRE(p.unmatchedModel == 0);
    REQUIRE(p.pairs.front().key.label() == "A:1");

    // Insertion codes are identity: 100 and 100A must not pair with each other.
    Structure a = makeChain(1);
    Structure b = makeChain(1);
    b.models[0].chains[0].residues[0].insertionCode = 'A';
    REQUIRE(pairResidues(*a.model(), *b.model()).pairs.empty());
}

TEST_CASE("SASA of an isolated carbon matches the closed-form sphere", "[bio]") {
    const SasaResult r = sasa(singleAtom("C"));
    const double analytic = 4.0 * kPi * (1.70 + 1.40) * (1.70 + 1.40);

    // Every one of the 92 points is accessible, so the quadrature is exact here: the
    // measured error against 4*pi*(R+probe)^2 = 120.762822 A^2 is 0.
    REQUIRE_THAT(r.total.value, WithinAbs(analytic, 1e-9));
    REQUIRE(r.total.unit == "A^2");
    REQUIRE(r.total.provenance == Provenance::Measured);
}

TEST_CASE("SASA of overlapping atoms is less than the sum of the isolated values", "[bio]") {
    Structure two = singleAtom("C");
    two.models[0].chains[0].residues[0].atoms.push_back(mkAtom(" C2 ", "C", 0.5, 0.0, 0.0));

    const double one = sasa(singleAtom("C")).total.value;
    const double both = sasa(two).total.value;
    REQUIRE(both < 2.0 * one);
    REQUIRE(both > one);   // two spheres 0.5 A apart still expose more than one alone
}

TEST_CASE("SASA states everything needed to reproduce it", "[bio]") {
    const SasaResult r = sasa(singleAtom("C"));
    // A bare SASA number is not reproducible; two tools disagree by ~10% on probe radius,
    // point count, radii set and hydrogen policy alone, so all four are in the source.
    REQUIRE_THAT(r.total.source, ContainsSubstring("Shrake-Rupley"));
    REQUIRE_THAT(r.total.source, ContainsSubstring("1.40"));
    REQUIRE_THAT(r.total.source, ContainsSubstring("92"));
    REQUIRE_THAT(r.total.source, ContainsSubstring("Bondi"));
    REQUIRE_THAT(r.total.source, ContainsSubstring("C 1.70"));
    REQUIRE_THAT(r.total.source, ContainsSubstring("united-atom"));
    REQUIRE(r.method == r.total.source);
}

TEST_CASE("SASA per-residue areas sum to the total and carry Tien 2013 relatives", "[bio]") {
    const Structure s = makeChain(20);
    const SasaResult r = sasa(s);

    double sum = 0.0;
    for (const ResidueSasa& rs : r.perResidue) {
        sum += rs.absolute.value;
    }
    REQUIRE_THAT(sum, WithinAbs(r.total.value, 1e-9));
    REQUIRE(r.perResidue.front().relative.provenance == Provenance::Measured);
    REQUIRE(r.perResidue.front().relative.unit.empty());
    REQUIRE_THAT(r.perResidue.front().relative.source, ContainsSubstring("Tien"));

    // A ligand has no published maximum accessibility, so there is no relative value to
    // report - an invented denominator would be worse than none.
    const SasaResult lig = sasa(singleAtom("C"));
    REQUIRE(lig.perResidue.front().relative.provenance == Provenance::NotComputed);
}

TEST_CASE("SASA hydrogen policy and the radii fallback are explicit", "[bio]") {
    Structure withH = singleAtom("C");
    withH.models[0].chains[0].residues[0].atoms.push_back(mkAtom(" H1 ", "H", 1.0, 0.0, 0.0));

    const double unitedAtom = sasa(withH).total.value;
    REQUIRE_THAT(unitedAtom, WithinAbs(sasa(singleAtom("C")).total.value, 1e-12));

    SasaOptions allAtom;
    allAtom.includeHydrogens = true;
    REQUIRE(sasa(withH, allAtom).total.value != unitedAtom);
    REQUIRE_THAT(sasa(withH, allAtom).total.source, ContainsSubstring("all-atom"));

    REQUIRE(vdwRadius("SE") == 1.90);
    REQUIRE(vdwRadius(" c ") == 1.70);
    REQUIRE(vdwRadius("UNOBTANIUM") == 1.80);   // documented fallback, never zero
    REQUIRE(maxAccessibility("TRP") == 285.0);
    REQUIRE(maxAccessibility("HOH") == 0.0);

    const Structure none;
    REQUIRE(sasa(none).total.provenance == Provenance::NotComputed);
}
