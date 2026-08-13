#include "bio/Cartoon.h"

#include <algorithm>
#include <cmath>

namespace biocad::bio {
namespace {

constexpr double kPi = 3.14159265358979323846;

Vec3d add(const Vec3d& a, const Vec3d& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
Vec3d sub(const Vec3d& a, const Vec3d& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
Vec3d scale(const Vec3d& a, double s) { return {a.x * s, a.y * s, a.z * s}; }
double dot(const Vec3d& a, const Vec3d& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
Vec3d cross(const Vec3d& a, const Vec3d& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double norm(const Vec3d& a) { return std::sqrt(dot(a, a)); }
Vec3d normalize(const Vec3d& a) {
    const double n = norm(a);
    return n > 1e-12 ? scale(a, 1.0 / n) : Vec3d{0, 0, 0};
}

// Any unit vector perpendicular to `t`. Used only when the peptide plane is unavailable AND the
// geometric fallback degenerates (three collinear CAs), which is rare and always reported.
Vec3d anyPerpendicular(const Vec3d& t) {
    const Vec3d ref = std::fabs(t.x) < 0.9 ? Vec3d{1, 0, 0} : Vec3d{0, 1, 0};
    return normalize(cross(ref, t));
}

const Atom* pick(const Residue& r, const char* name) {
    // Residue::atom() already compares on the trimmed name and returns the first match, which
    // for a disordered residue is the A conformer. The cartoon follows one conformer; drawing
    // both would give the chain two backbones.
    return r.atom(name);
}

Vec3d position(const Atom& a) { return {a.x, a.y, a.z}; }

bool insideRange(int seqId, char ins, int startSeq, char startIns, int endSeq, char endIns) {
    // Insertion codes make residue numbers a lexicographic pair, not an integer: 100 < 100A.
    const auto key = [](int s, char i) { return std::pair<int, char>{s, i == '\0' ? ' ' : i}; };
    const auto k = key(seqId, ins);
    return k >= key(startSeq, startIns) && k <= key(endSeq, endIns);
}

struct RingPoint {
    Vec3d offset;   // in the (normal, binormal) plane, Angstrom
    Vec3d normal;   // outward surface normal, same plane
};

// One closed cross-section ring. The WIDE axis is the ribbon normal (the peptide-plane
// direction) and the THIN axis is the binormal, so a strand's flat face contains the peptide
// planes - which is what makes a sheet look like a sheet.
std::vector<RingPoint> ring(SsType ss, int radialSegments, double sizeFactor, double aspectRatio,
                            double widthScale) {
    std::vector<RingPoint> out;
    out.reserve(static_cast<std::size_t>(radialSegments));
    const double thin = sizeFactor;
    const double wide = sizeFactor * aspectRatio;

    if (ss == SsType::Strand) {
        // Flat rectangle: four straight sides with per-side normals, so the edges stay crisp
        // instead of being smoothed into a lozenge. Corner vertices take their side's normal,
        // which leaves a one-vertex shading seam at each corner - correct for a flat ribbon.
        const double w = wide * widthScale;
        const double h = thin;
        const int perSide = std::max(1, radialSegments / 4);
        const int total = perSide * 4;
        struct Side { Vec3d from, to, normal; };
        const Side sides[4] = {
            {{ w,  h, 0}, {-w,  h, 0}, {0, 1, 0}},
            {{-w,  h, 0}, {-w, -h, 0}, {-1, 0, 0}},
            {{-w, -h, 0}, { w, -h, 0}, {0, -1, 0}},
            {{ w, -h, 0}, { w,  h, 0}, {1, 0, 0}},
        };
        for (int k = 0; k < total; ++k) {
            const Side& s = sides[k / perSide];
            const double f = static_cast<double>(k % perSide) / perSide;
            RingPoint p;
            p.offset = add(s.from, scale(sub(s.to, s.from), f));
            p.normal = s.normal;
            out.push_back(p);
        }
        return out;
    }

    // Helix: ellipse (wide x thin). Coil: circle of the base thickness - a tube, because a
    // flat loop would imply a plane the residue does not have.
    const double a = ss == SsType::Helix ? wide * widthScale : thin;
    const double b = thin;
    for (int k = 0; k < radialSegments; ++k) {
        const double th = 2.0 * kPi * static_cast<double>(k) / radialSegments;
        const double c = std::cos(th), s = std::sin(th);
        RingPoint p;
        p.offset = {a * c, b * s, 0};
        // Ellipse outward normal at (a cos, b sin) is proportional to (b cos, a sin).
        p.normal = normalize(Vec3d{b * c, a * s, 0});
        out.push_back(p);
    }
    return out;
}

std::uint32_t colorIndexFor(SsType ss) {
    switch (ss) {
        case SsType::Helix: return 1;
        case SsType::Strand: return 2;
        case SsType::Coil: break;
    }
    return 0;
}

}  // namespace

std::vector<SsType> assignFromAnnotations(const Chain& chain, const Annotations& annotations) {
    std::vector<SsType> ss(chain.residues.size(), SsType::Coil);
    for (std::size_t i = 0; i < chain.residues.size(); ++i) {
        const Residue& r = chain.residues[i];
        for (const HelixRecord& h : annotations.helices) {
            if (h.chainId != chain.id) continue;
            if (insideRange(r.authSeqId, r.insertionCode, h.startSeqId, h.startInsertionCode,
                            h.endSeqId, h.endInsertionCode)) {
                ss[i] = SsType::Helix;
                break;
            }
        }
        if (ss[i] != SsType::Coil) continue;
        for (const StrandRecord& s : annotations.strands) {
            if (s.chainId != chain.id) continue;
            if (insideRange(r.authSeqId, r.insertionCode, s.startSeqId, s.startInsertionCode,
                            s.endSeqId, s.endInsertionCode)) {
                ss[i] = SsType::Strand;
                break;
            }
        }
    }
    return ss;
}

std::vector<GuidePoint> guidePoints(const Chain& chain, const std::vector<SsType>& ss,
                                    const CartoonOptions& options) {
    std::vector<GuidePoint> guides;
    guides.reserve(chain.residues.size());
    for (std::size_t i = 0; i < chain.residues.size(); ++i) {
        const Residue& r = chain.residues[i];
        const Atom* ca = pick(r, "CA");
        // A residue with no CA is not a polymer position for cartoon purposes: waters, ions and
        // ligands are drawn by the stick path, not by the ribbon.
        if (!ca || r.oneLetter() == 'X') continue;
        GuidePoint g;
        g.ca = position(*ca);
        g.position = g.ca;
        g.residueIndex = static_cast<int>(i);
        g.ss = i < ss.size() ? ss[i] : SsType::Coil;
        const Atom* c = pick(r, "C");
        const Atom* o = pick(r, "O");
        if (c && o) {
            g.normal = normalize(sub(position(*o), position(*c)));
        } else {
            g.normalFromGeometry = true;   // filled in the neighbour pass below
        }
        guides.push_back(g);
    }
    if (guides.size() < 2) return {};

    // Geometric fallback for residues whose carbonyl is absent: the normal of the plane through
    // three consecutive CAs. Reported through GuidePoint::normalFromGeometry because it is a
    // different quantity from the peptide-plane normal, not a silent substitute.
    for (std::size_t i = 0; i < guides.size(); ++i) {
        if (!guides[i].normalFromGeometry && norm(guides[i].normal) > 1e-12) continue;
        guides[i].normalFromGeometry = true;
        Vec3d n{0, 0, 0};
        if (i > 0 && i + 1 < guides.size())
            n = cross(sub(guides[i - 1].ca, guides[i].ca), sub(guides[i + 1].ca, guides[i].ca));
        if (norm(n) <= 1e-12) {
            const std::size_t j = i + 1 < guides.size() ? i + 1 : i - 1;
            n = anyPerpendicular(normalize(sub(guides[j].ca, guides[i].ca)));
        }
        guides[i].normal = normalize(n);
    }

    // THE FLIP CHECK. Successive peptide planes alternate direction; interpolating the raw
    // normals turns a flat strand into a corkscrew (~180 degrees of twist per residue). Each
    // normal is compared with its ALREADY-CORRECTED predecessor, so one flip does not undo the
    // next. Disabling this is a rendering bug and exists only to be measured.
    if (options.flipCheck) {
        for (std::size_t i = 1; i < guides.size(); ++i) {
            if (dot(guides[i].normal, guides[i - 1].normal) < 0.0) {
                guides[i].normal = scale(guides[i].normal, -1.0);
                guides[i].flipped = true;
            }
        }
    }

    // Guide-point smoothing for helices: a helical CA sits ~2.3 A off the axis, so the raw CA
    // path is a coil of the same pitch as the helix. Moving the point standardShift of the way
    // to the midpoint of its neighbours pulls it towards the axis.
    if (options.standardShift > 0.0) {
        std::vector<Vec3d> shifted(guides.size());
        for (std::size_t i = 0; i < guides.size(); ++i) shifted[i] = guides[i].ca;
        for (std::size_t i = 1; i + 1 < guides.size(); ++i) {
            if (guides[i].ss != SsType::Helix) continue;
            const Vec3d mid = scale(add(guides[i - 1].ca, guides[i + 1].ca), 0.5);
            shifted[i] = add(scale(guides[i].ca, 1.0 - options.standardShift),
                             scale(mid, options.standardShift));
        }
        for (std::size_t i = 0; i < guides.size(); ++i) guides[i].position = shifted[i];
    }
    return guides;
}

std::vector<FrameSample> sampleSpline(const std::vector<GuidePoint>& guides,
                                      const CartoonOptions& options) {
    std::vector<FrameSample> out;
    const int n = static_cast<int>(guides.size());
    if (n < 2 || options.linearSegments < 1) return out;

    // Phantom control points past each terminus, so the first and last real segments have a
    // defined tangent. overhangFactor / 2 = 1 reflects the neighbour exactly one spacing out.
    const double over = options.overhangFactor * 0.5;
    const Vec3d pre = add(guides[0].position,
                          scale(sub(guides[0].position, guides[1].position), over));
    const Vec3d post = add(guides[static_cast<std::size_t>(n - 1)].position,
                           scale(sub(guides[static_cast<std::size_t>(n - 1)].position,
                                     guides[static_cast<std::size_t>(n - 2)].position), over));
    const auto P = [&](int i) -> Vec3d {
        if (i < 0) return pre;
        if (i >= n) return post;
        return guides[static_cast<std::size_t>(i)].position;
    };
    const auto N = [&](int i) -> Vec3d {
        const int c = std::clamp(i, 0, n - 1);
        return guides[static_cast<std::size_t>(c)].normal;
    };

    // A strand's arrowhead sits over the segment that ENDS the strand run. When the run ends at
    // the chain terminus there is no following segment, so the last segment carries the arrow.
    std::vector<bool> arrowSegment(static_cast<std::size_t>(std::max(0, n - 1)), false);
    for (int i = 0; i < n; ++i) {
        if (guides[static_cast<std::size_t>(i)].ss != SsType::Strand) continue;
        const bool runEnds =
            i + 1 >= n || guides[static_cast<std::size_t>(i + 1)].ss != SsType::Strand;
        if (!runEnds) continue;
        const int seg = std::min(i, n - 2);
        if (seg >= 0) arrowSegment[static_cast<std::size_t>(seg)] = true;
    }

    const auto emit = [&](int seg, double u, bool knot) {
        const Vec3d p0 = P(seg - 1), p1 = P(seg), p2 = P(seg + 1), p3 = P(seg + 2);
        const bool helix = guides[static_cast<std::size_t>(seg)].ss == SsType::Helix ||
                           guides[static_cast<std::size_t>(std::min(seg + 1, n - 1))].ss ==
                               SsType::Helix;
        const double tension = helix ? options.helixTension : options.standardTension;
        const Vec3d m1 = scale(sub(p2, p0), tension);
        const Vec3d m2 = scale(sub(p3, p1), tension);
        const double u2 = u * u, u3 = u2 * u;
        // Cubic Hermite. h00(0) = 1 and every other basis function vanishes at u = 0, so a
        // knot sample IS the guide point - exactly, not approximately.
        const double h00 = 2 * u3 - 3 * u2 + 1;
        const double h10 = u3 - 2 * u2 + u;
        const double h01 = -2 * u3 + 3 * u2;
        const double h11 = u3 - u2;
        FrameSample s;
        s.position = add(add(scale(p1, h00), scale(m1, h10)), add(scale(p2, h01), scale(m2, h11)));
        const double d00 = 6 * u2 - 6 * u;
        const double d10 = 3 * u2 - 4 * u + 1;
        const double d01 = -6 * u2 + 6 * u;
        const double d11 = 3 * u2 - 2 * u;
        Vec3d der = add(add(scale(p1, d00), scale(m1, d10)), add(scale(p2, d01), scale(m2, d11)));
        if (norm(der) <= 1e-12) der = sub(p2, p1);
        s.tangent = normalize(der);
        // The ribbon normal is interpolated between the two guide normals and then made
        // perpendicular to the tangent by Gram-Schmidt, which is what keeps the frame
        // orthonormal for every sample instead of only at the knots.
        Vec3d nn = add(scale(N(seg), 1.0 - u), scale(N(seg + 1), u));
        nn = sub(nn, scale(s.tangent, dot(nn, s.tangent)));
        if (norm(nn) <= 1e-9) nn = anyPerpendicular(s.tangent);
        s.normal = normalize(nn);
        s.binormal = normalize(cross(s.tangent, s.normal));
        s.ss = guides[static_cast<std::size_t>(seg)].ss;
        s.guideIndex = seg;
        s.u = u;
        s.knot = knot;
        s.widthScale = 1.0;
        if (s.ss == SsType::Strand && seg < static_cast<int>(arrowSegment.size()) &&
            arrowSegment[static_cast<std::size_t>(seg)]) {
            // The arrowhead: the ribbon widens to arrowFactor at the base and tapers to a point
            // at the strand's C-terminal end, which is what makes the direction readable.
            s.widthScale = options.arrowFactor * (1.0 - u);
        }
        out.push_back(s);
    };

    for (int seg = 0; seg + 1 < n; ++seg) {
        for (int j = 0; j < options.linearSegments; ++j) {
            const double u = static_cast<double>(j) / options.linearSegments;
            emit(seg, u, j == 0);
        }
    }
    emit(n - 2, 1.0, true);   // the final knot, so the ribbon reaches the last residue
    return out;
}

double totalTwistDegrees(const std::vector<FrameSample>& samples) {
    double total = 0.0;
    for (std::size_t i = 1; i < samples.size(); ++i) {
        // Rotation of the normal about the shared tangent between consecutive samples. Both
        // normals are projected into the plane perpendicular to the averaged tangent so the
        // angle measures TWIST and not the curve's own bending.
        const Vec3d t = normalize(add(samples[i - 1].tangent, samples[i].tangent));
        if (norm(t) <= 1e-12) continue;
        Vec3d a = sub(samples[i - 1].normal, scale(t, dot(samples[i - 1].normal, t)));
        Vec3d b = sub(samples[i].normal, scale(t, dot(samples[i].normal, t)));
        if (norm(a) <= 1e-12 || norm(b) <= 1e-12) continue;
        a = normalize(a);
        b = normalize(b);
        const double c = std::clamp(dot(a, b), -1.0, 1.0);
        total += std::acos(c) * 180.0 / kPi;
    }
    return total;
}

CartoonMesh buildCartoon(const Chain& chain, const std::vector<SsType>& ss,
                         const CartoonOptions& options) {
    CartoonMesh mesh;
    // 0xAABBGGRR, matching render::AtomInst::rgba so the mesh path shares the shader's packing.
    mesh.palette = {0xFFB4B4B4u /* coil, grey */, 0xFF4C4CE6u /* helix, red */,
                    0xFF32C8F0u /* strand, amber */};

    const std::vector<GuidePoint> guides = guidePoints(chain, ss, options);
    if (guides.size() < 2) {
        mesh.warnings.push_back("chain " + chain.id +
                                " has fewer than two residues with a CA atom: no ribbon");
        return mesh;
    }
    mesh.residues = guides.size();
    const std::vector<FrameSample> samples = sampleSpline(guides, options);
    if (samples.size() < 2) return mesh;

    // Rounded DOWN to a multiple of four so the rectangular strand cross-section (four sides,
    // equal point counts) emits exactly as many vertices per ring as the elliptical one. Rings
    // of different sizes would break the fixed-stride stitch below.
    const int R = std::max(4, (options.radialSegments / 4) * 4);
    const std::size_t ringCount = samples.size();
    mesh.vertices.reserve(ringCount * static_cast<std::size_t>(R) + 2);
    mesh.indices.reserve((ringCount - 1) * static_cast<std::size_t>(R) * 6 +
                         static_cast<std::size_t>(R) * 6);

    for (const FrameSample& s : samples) {
        const std::vector<RingPoint> rp =
            ring(s.ss, R, options.sizeFactor, options.aspectRatio, s.widthScale);
        const std::uint32_t ci = colorIndexFor(s.ss);
        for (const RingPoint& p : rp) {
            // The ring is authored in the frame's own 2D basis: x along the ribbon normal
            // (wide), y along the binormal (thin). Both offsets and normals map through it.
            const Vec3d world = add(s.position, add(scale(s.normal, p.offset.x),
                                                    scale(s.binormal, p.offset.y)));
            const Vec3d nrm = normalize(add(scale(s.normal, p.normal.x),
                                            scale(s.binormal, p.normal.y)));
            CartoonVertex v;
            v.px = static_cast<float>(world.x);
            v.py = static_cast<float>(world.y);
            v.pz = static_cast<float>(world.z);
            v.nx = static_cast<float>(nrm.x);
            v.ny = static_cast<float>(nrm.y);
            v.nz = static_cast<float>(nrm.z);
            v.colorIndex = ci;
            mesh.vertices.push_back(v);
        }
    }

    // Stitch consecutive rings. Every ring emits exactly R vertices (see the clamp above), so
    // the stitch is plain index arithmetic over a fixed stride.
    const std::size_t stride = mesh.vertices.size() / ringCount;
    for (std::size_t r = 0; r + 1 < ringCount; ++r) {
        const std::size_t a0 = r * stride;
        const std::size_t b0 = (r + 1) * stride;
        for (std::size_t k = 0; k < stride; ++k) {
            const std::size_t k1 = (k + 1) % stride;
            const auto A = static_cast<std::uint32_t>(a0 + k);
            const auto B = static_cast<std::uint32_t>(a0 + k1);
            const auto C = static_cast<std::uint32_t>(b0 + k);
            const auto D = static_cast<std::uint32_t>(b0 + k1);
            mesh.indices.insert(mesh.indices.end(), {A, C, B});
            mesh.indices.insert(mesh.indices.end(), {B, C, D});
        }
    }

    if (options.caps) {
        auto cap = [&](std::size_t ringIndex, const Vec3d& normal, bool reverse) {
            const std::size_t base = ringIndex * stride;
            Vec3d centre{0, 0, 0};
            for (std::size_t k = 0; k < stride; ++k) {
                const CartoonVertex& v = mesh.vertices[base + k];
                centre = add(centre, Vec3d{v.px, v.py, v.pz});
            }
            centre = scale(centre, 1.0 / static_cast<double>(stride));
            const std::uint32_t colour = mesh.vertices[base].colorIndex;
            const auto hub = static_cast<std::uint32_t>(mesh.vertices.size());
            CartoonVertex c;
            c.px = static_cast<float>(centre.x);
            c.py = static_cast<float>(centre.y);
            c.pz = static_cast<float>(centre.z);
            c.nx = static_cast<float>(normal.x);
            c.ny = static_cast<float>(normal.y);
            c.nz = static_cast<float>(normal.z);
            c.colorIndex = colour;
            mesh.vertices.push_back(c);
            // The rim vertices are reused rather than duplicated with the cap normal: the cap
            // is 0.4 A across at the default sizeFactor, so a shading seam there is invisible
            // and the alternative doubles the vertex count of every chain terminus.
            for (std::size_t k = 0; k < stride; ++k) {
                const auto A = static_cast<std::uint32_t>(base + k);
                const auto B = static_cast<std::uint32_t>(base + (k + 1) % stride);
                if (reverse)
                    mesh.indices.insert(mesh.indices.end(), {hub, B, A});
                else
                    mesh.indices.insert(mesh.indices.end(), {hub, A, B});
            }
        };
        cap(0, scale(samples.front().tangent, -1.0), true);
        cap(ringCount - 1, samples.back().tangent, false);
    }
    return mesh;
}

CartoonMesh buildCartoon(const Model& model, const Annotations& annotations,
                         const CartoonOptions& options) {
    CartoonMesh merged;
    merged.palette = {0xFFB4B4B4u, 0xFF4C4CE6u, 0xFF32C8F0u};
    for (const Chain& chain : model.chains) {
        const std::vector<SsType> ss = assignFromAnnotations(chain, annotations);
        CartoonMesh part = buildCartoon(chain, ss, options);
        for (const std::string& w : part.warnings) merged.warnings.push_back(w);
        if (part.empty()) continue;
        const auto offset = static_cast<std::uint32_t>(merged.vertices.size());
        merged.vertices.insert(merged.vertices.end(), part.vertices.begin(), part.vertices.end());
        for (std::uint32_t i : part.indices) merged.indices.push_back(i + offset);
        merged.residues += part.residues;
    }
    return merged;
}

}  // namespace biocad::bio
