// Cartoon / ribbon geometry. Every claim here is checkable without a GPU, so it is checked:
// the spline interpolates its guide points exactly, the frames stay orthonormal, and THE FLIP
// CHECK IS MEASURED - the total twist of the ribbon normal about its tangent is computed with
// the check on and off and the two numbers are compared, rather than asserting a flag.
//
// The fixture is ideal_helix_strand.pdb: 12 residues of ideal alpha helix, a 3-residue coil and
// 10 residues of ideal extended strand, built from standard bond lengths and angles.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <filesystem>
#include <map>

#include "bio/Cartoon.h"
#include "bio/PdbReader.h"

using namespace biocad::bio;

namespace {

std::filesystem::path fixture(const char* name) {
    return std::filesystem::path(BIOCAD_TEST_FIXTURES) / name;
}

double length(const Vec3d& v) { return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z); }
double dot(const Vec3d& a, const Vec3d& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

struct Loaded {
    Structure   structure;
    Annotations annotations;
};

Loaded load() {
    Loaded l;
    l.structure = readPdbFile(fixture("ideal_helix_strand.pdb"), &l.annotations);
    return l;
}

// A sub-chain of `chain` with a forced secondary structure, so the helix and strand paths can
// be exercised in isolation.
Chain slice(const Chain& chain, std::size_t first, std::size_t last) {
    Chain out;
    out.id = chain.id;
    for (std::size_t i = first; i < last; ++i) out.residues.push_back(chain.residues[i]);
    return out;
}

}  // namespace

TEST_CASE("secondary structure comes from the HELIX and SHEET records", "[bio][cartoon]") {
    const Loaded l = load();
    REQUIRE(l.structure.warnings.empty());
    const Model* m = l.structure.model(1);
    REQUIRE(m != nullptr);
    REQUIRE(l.annotations.helices.size() == 1);
    REQUIRE(l.annotations.strands.size() == 1);

    const std::vector<SsType> ss = assignFromAnnotations(m->chains.at(0), l.annotations);
    REQUIRE(ss.size() == 25);
    int helix = 0, strand = 0, coil = 0;
    for (SsType t : ss) {
        if (t == SsType::Helix) ++helix;
        else if (t == SsType::Strand) ++strand;
        else ++coil;
    }
    REQUIRE(helix == 12);
    REQUIRE(strand == 10);
    REQUIRE(coil == 3);
}

TEST_CASE("the spline passes exactly through its guide points", "[bio][cartoon]") {
    const Loaded l = load();
    const Chain& chain = l.structure.model(1)->chains.at(0);
    const std::vector<SsType> ss = assignFromAnnotations(chain, l.annotations);
    const CartoonOptions opts;
    const std::vector<GuidePoint> guides = guidePoints(chain, ss, opts);
    REQUIRE(guides.size() == 25);

    const std::vector<FrameSample> samples = sampleSpline(guides, opts);
    // linearSegments samples per segment, plus the closing knot.
    REQUIRE(samples.size() == (guides.size() - 1) * 8 + 1);

    std::size_t knots = 0;
    for (const FrameSample& s : samples) {
        if (!s.knot) continue;
        ++knots;
        const std::size_t gi =
            static_cast<std::size_t>(s.u == 0.0 ? s.guideIndex : s.guideIndex + 1);
        const GuidePoint& g = guides[gi];
        const Vec3d d{s.position.x - g.position.x, s.position.y - g.position.y,
                      s.position.z - g.position.z};
        REQUIRE(length(d) < 1e-12);
    }
    REQUIRE(knots == guides.size());
}

TEST_CASE("every spline frame is orthonormal", "[bio][cartoon]") {
    const Loaded l = load();
    const Chain& chain = l.structure.model(1)->chains.at(0);
    const std::vector<SsType> ss = assignFromAnnotations(chain, l.annotations);
    const CartoonOptions opts;
    const std::vector<FrameSample> samples = sampleSpline(guidePoints(chain, ss, opts), opts);
    REQUIRE(samples.size() > 100);

    for (const FrameSample& s : samples) {
        REQUIRE_THAT(length(s.tangent), Catch::Matchers::WithinAbs(1.0, 1e-12));
        REQUIRE_THAT(length(s.normal), Catch::Matchers::WithinAbs(1.0, 1e-12));
        REQUIRE_THAT(length(s.binormal), Catch::Matchers::WithinAbs(1.0, 1e-12));
        REQUIRE(std::fabs(dot(s.tangent, s.normal)) < 1e-12);
        REQUIRE(std::fabs(dot(s.tangent, s.binormal)) < 1e-12);
        REQUIRE(std::fabs(dot(s.normal, s.binormal)) < 1e-12);
    }
}

TEST_CASE("the 180-degree flip check is measured, not asserted", "[bio][cartoon]") {
    const Loaded l = load();
    const Chain& chain = l.structure.model(1)->chains.at(0);

    // The strand is where successive peptide planes are most nearly anti-parallel, so it is
    // where the corkscrew is unambiguous: without the flip check the ribbon normal rotates
    // about the tangent by roughly 180 degrees PER RESIDUE.
    const Chain strand = slice(chain, 15, 25);
    const std::vector<SsType> ss(strand.residues.size(), SsType::Strand);

    CartoonOptions on;
    CartoonOptions off;
    off.flipCheck = false;

    const std::vector<GuidePoint> gOn = guidePoints(strand, ss, on);
    const std::vector<GuidePoint> gOff = guidePoints(strand, ss, off);
    REQUIRE(gOn.size() == 10);
    REQUIRE(gOff.size() == 10);
    // Roughly every other normal needs negating; with the check off, none are touched.
    int flipped = 0;
    for (const GuidePoint& g : gOn) flipped += g.flipped ? 1 : 0;
    REQUIRE(flipped >= 4);
    for (const GuidePoint& g : gOff) REQUIRE_FALSE(g.flipped);

    const double twistOn = totalTwistDegrees(sampleSpline(gOn, on));
    const double twistOff = totalTwistDegrees(sampleSpline(gOff, off));
    const double perResidueOn = twistOn / static_cast<double>(gOn.size() - 1);
    const double perResidueOff = twistOff / static_cast<double>(gOff.size() - 1);

    REQUIRE(perResidueOff > 160.0);   // measured 184.5 deg/residue on this fixture
    REQUIRE(perResidueOn < 20.0);     // measured 6.3 deg/residue
    REQUIRE(twistOff > 10.0 * twistOn);
}

TEST_CASE("mesh counts are 256 triangles per residue to the expected order", "[bio][cartoon]") {
    const Loaded l = load();
    const Model* m = l.structure.model(1);
    const Chain& chain = m->chains.at(0);
    const CartoonOptions opts;

    const Chain helix = slice(chain, 0, 12);
    const Chain strand = slice(chain, 15, 25);
    const CartoonMesh mh =
        buildCartoon(helix, std::vector<SsType>(helix.residues.size(), SsType::Helix), opts);
    const CartoonMesh msd =
        buildCartoon(strand, std::vector<SsType>(strand.residues.size(), SsType::Strand), opts);

    // 8 samples/residue x 16 radial segments x 2 triangles = 256 per residue interval, plus a
    // 16-triangle cap at each end.
    REQUIRE(mh.residues == 12);
    REQUIRE(mh.triangleCount() == (12 - 1) * 8 * 16 * 2 + 32);
    REQUIRE(msd.residues == 10);
    REQUIRE(msd.triangleCount() == (10 - 1) * 8 * 16 * 2 + 32);
    // 2848 / 12 = 237 and 2336 / 10 = 234: the 256 figure is per residue INTERVAL, so a short
    // segment reads slightly under it and a 300-residue chain lands at 76576 triangles.
    REQUIRE(mh.triangleCount() == 2848);
    REQUIRE(msd.triangleCount() == 2336);

    for (const CartoonMesh* mesh : {&mh, &msd}) {
        REQUIRE(mesh->indices.size() % 3 == 0);
        REQUIRE(mesh->palette.size() == 3);
        for (std::uint32_t i : mesh->indices) REQUIRE(i < mesh->vertices.size());
        for (const CartoonVertex& v : mesh->vertices) {
            const double n = std::sqrt(static_cast<double>(v.nx) * v.nx +
                                       static_cast<double>(v.ny) * v.ny +
                                       static_cast<double>(v.nz) * v.nz);
            REQUIRE_THAT(n, Catch::Matchers::WithinAbs(1.0, 1e-5));
        }
    }

    // Colour indices follow the secondary structure, so all three palette entries appear on a
    // chain that has all three states.
    const std::vector<SsType> ss = assignFromAnnotations(chain, l.annotations);
    const CartoonMesh whole = buildCartoon(chain, ss, opts);
    std::map<std::uint32_t, int> used;
    for (const CartoonVertex& v : whole.vertices) ++used[v.colorIndex];
    REQUIRE(used.size() == 3);
    REQUIRE(whole.triangleCount() == (25 - 1) * 8 * 16 * 2 + 32);

    // The whole-model overload merges chains with the index offsets applied.
    const CartoonMesh merged = buildCartoon(*m, l.annotations, opts);
    REQUIRE(merged.residues == whole.residues);
    REQUIRE(merged.triangleCount() == whole.triangleCount());
    for (std::uint32_t i : merged.indices) REQUIRE(i < merged.vertices.size());
}

TEST_CASE("a strand carries an arrowhead that tapers to a point", "[bio][cartoon]") {
    const Loaded l = load();
    const Chain& chain = l.structure.model(1)->chains.at(0);
    const std::vector<SsType> ss = assignFromAnnotations(chain, l.annotations);
    const CartoonOptions opts;
    const std::vector<FrameSample> samples = sampleSpline(guidePoints(chain, ss, opts), opts);

    double widest = 0.0, narrowest = 1.0;
    for (const FrameSample& s : samples) {
        if (s.ss != SsType::Strand) continue;
        widest = std::max(widest, s.widthScale);
        narrowest = std::min(narrowest, s.widthScale);
    }
    REQUIRE_THAT(widest, Catch::Matchers::WithinAbs(opts.arrowFactor, 1e-12));
    REQUIRE(narrowest < 0.2);
}

TEST_CASE("a chain with fewer than two CA atoms reports instead of drawing", "[bio][cartoon]") {
    Chain one;
    one.id = "Z";
    Residue r;
    r.name = "ALA";
    r.authSeqId = 1;
    Atom ca;
    ca.name = " CA ";
    ca.element = "C";
    r.atoms.push_back(ca);
    one.residues.push_back(r);

    const CartoonMesh mesh = buildCartoon(one, {SsType::Coil});
    REQUIRE(mesh.empty());
    REQUIRE(mesh.warnings.size() == 1);
    REQUIRE(mesh.warnings[0].find("fewer than two") != std::string::npos);
}
