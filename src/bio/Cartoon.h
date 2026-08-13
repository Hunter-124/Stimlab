// bio/Cartoon.h - cartoon / ribbon geometry for a protein chain, as a plain indexed mesh.
//
// The pipeline is four stages, each separately testable without a GPU:
//
//   1. GUIDE POINTS. One per residue, from CA, plus a ribbon normal from the peptide plane
//      (the carbonyl direction O-C). This stage owns the FLIP CHECK described below.
//   2. SPLINE. A cardinal (Catmull-Rom family) spline through the guide points with a
//      per-sample orthonormal frame (tangent, normal, binormal). The spline INTERPOLATES its
//      guide points exactly - a sample at a knot is the CA-derived guide point itself.
//   3. CROSS-SECTIONS. A closed ring of points in the frame's normal/binormal plane whose
//      shape depends on the secondary structure: ellipse for helix, flat rectangle for strand
//      (with an arrowhead over the last residue of the strand), small circle for coil.
//   4. EXTRUSION. Consecutive rings are stitched into triangles; the chain ends get caps.
//
// THE FLIP CHECK, AND WHY IT IS NOT OPTIONAL. Successive peptide planes along a chain point in
// alternating directions: in an extended strand the carbonyl of residue i points one way and
// residue i+1's points nearly the opposite way (~180 degrees apart). Interpolating those raw
// normals rotates the ribbon by half a turn per residue, so a flat strand renders as a
// corkscrew. The fix is to compare each normal with its predecessor and negate it when their
// dot product is negative. CartoonOptions::flipCheck exists so the failure can be MEASURED
// rather than asserted: bio::totalTwistDegrees() over the frames of the same chain, with the
// check on and off, is the number the test compares.
//
// SECONDARY STRUCTURE IS NOT ASSIGNED HERE. It comes in as a per-residue vector, built by
// assignFromAnnotations() out of the HELIX / SHEET (or _struct_conf / _struct_sheet_range)
// records the readers already parse. This file deliberately contains no third assignment
// algorithm: a cartoon that disagreed with the panel's own secondary structure would be two
// answers to one question.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "bio/Annotations.h"
#include "bio/Structure.h"

namespace biocad::bio {

enum class SsType { Coil, Helix, Strand };

// Mol*'s cartoon defaults, with the role each one plays HERE.
struct CartoonOptions {
    // Cardinal-spline tangent scale. 0.5 is the classic uniform Catmull-Rom tangent
    // (P[i+1]-P[i-1])/2 and is used for coil and strand.
    double standardTension = 0.5;
    // The same scale for helices. A larger value lengthens the Hermite tangents, which keeps
    // the ribbon hugging the helical path instead of cutting the corner at every turn.
    double helixTension = 0.9;
    // Guide-point smoothing: a helix CA sits off the helix axis, so the guide point is moved
    // this fraction of the way towards the midpoint of its two neighbours. 0 = raw CA.
    double standardShift = 0.5;
    // Length of the phantom control point beyond each chain terminus, in units of half a
    // residue spacing, so the first and last real segments have a defined tangent. 2 =
    // reflect the neighbour exactly one residue spacing past the terminus.
    double overhangFactor = 2.0;
    // Base half-thickness of every cross-section, in Angstrom.
    double sizeFactor = 0.2;
    // Width / thickness of the flat cross-sections: a strand ribbon is 5x wider than thick.
    double aspectRatio = 5.0;
    // Width multiplier at the base of a strand's arrowhead.
    double arrowFactor = 1.5;
    // Spline samples per residue. 8 samples x 16 radial x 2 triangles = 256 triangles/residue.
    int linearSegments = 8;
    // Points around one cross-section ring.
    int radialSegments = 16;
    // See the header comment. Off is a rendering bug, kept reachable only so it can be measured.
    bool flipCheck = true;
    // Flat caps on the first and last ring so a ribbon is not an open tube.
    bool caps = true;
};

struct Vec3d {
    double x = 0, y = 0, z = 0;
};

struct GuidePoint {
    Vec3d  position;      // CA, optionally shifted towards the neighbour midpoint
    Vec3d  ca;            // the raw CA, kept so the spline can be checked against the input
    Vec3d  normal;        // ribbon normal (peptide plane), flip-corrected when enabled
    SsType ss = SsType::Coil;
    int    residueIndex = 0;   // index into Chain::residues
    bool   normalFromGeometry = false;  // true when O/C were absent and a fallback was used
    bool   flipped = false;             // this normal was negated by the flip check
};

struct FrameSample {
    Vec3d  position;
    Vec3d  tangent;
    Vec3d  normal;
    Vec3d  binormal;
    SsType ss = SsType::Coil;
    int    guideIndex = 0;   // the guide point this sample's segment starts at
    double u = 0.0;          // 0..1 within that segment
    bool   knot = false;     // u == 0: the sample sits exactly on the guide point
    double widthScale = 1.0; // 1, or the arrowhead taper over a strand's last residue
};

struct CartoonVertex {
    float px = 0, py = 0, pz = 0;
    float nx = 0, ny = 0, nz = 0;
    std::uint32_t colorIndex = 0;   // index into CartoonMesh::palette
};

struct CartoonMesh {
    std::vector<CartoonVertex>  vertices;
    std::vector<std::uint32_t>  indices;      // triangle list
    std::vector<std::uint32_t>  palette;      // 0xAABBGGRR, indexed by CartoonVertex::colorIndex
    std::size_t                 residues = 0; // residues that contributed a guide point
    std::vector<std::string>    warnings;

    [[nodiscard]] std::size_t triangleCount() const { return indices.size() / 3; }
    [[nodiscard]] bool empty() const { return indices.empty(); }
};

// Per-residue secondary structure from the HELIX / SHEET records the readers already parsed.
// Residues outside every record are Coil - that is the honest default, because a PDB entry
// without HELIX records has not stated that its helices are absent, only that it did not list
// them, and the panel says so.
std::vector<SsType> assignFromAnnotations(const Chain& chain, const Annotations& annotations);

// Guide points for one chain. Residues without a CA (waters, ions, ligands) are skipped, which
// is why GuidePoint carries its residue index. Fewer than two guide points yields an empty
// vector: a single point defines no curve.
std::vector<GuidePoint> guidePoints(const Chain& chain, const std::vector<SsType>& ss,
                                    const CartoonOptions& options = {});

// Spline samples with orthonormal frames. Samples are ordered along the chain; every segment
// contributes options.linearSegments samples plus one final sample at the very end.
std::vector<FrameSample> sampleSpline(const std::vector<GuidePoint>& guides,
                                      const CartoonOptions& options = {});

// Sum over consecutive samples of the absolute rotation of the ribbon normal about the local
// tangent, in degrees. This is the corkscrew metric: with the flip check on it stays small,
// with it off it accumulates roughly 180 degrees per residue.
double totalTwistDegrees(const std::vector<FrameSample>& samples);

// The whole pipeline for one chain.
CartoonMesh buildCartoon(const Chain& chain, const std::vector<SsType>& ss,
                         const CartoonOptions& options = {});

// Every chain of a model, merged into one mesh (index offsets applied). Chains with fewer than
// two guide points are reported in CartoonMesh::warnings, not silently dropped.
CartoonMesh buildCartoon(const Model& model, const Annotations& annotations,
                         const CartoonOptions& options = {});

}  // namespace biocad::bio
