#include "chem/Embed3D.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <vector>

#include <Eigen/Dense>

namespace biocad::chem {
namespace {

// --------------------------------------------------------------- element data
// Single-bond covalent radii (Cordero 2008), Angstrom. Carbon-like default.
double covRadius(int z) {
    switch (z) {
        case 1:  return 0.31;  case 5:  return 0.84;  case 6:  return 0.76;
        case 7:  return 0.71;  case 8:  return 0.66;  case 9:  return 0.57;
        case 11: return 1.66;  case 12: return 1.41;  case 13: return 1.21;
        case 14: return 1.11;  case 15: return 1.07;  case 16: return 1.05;
        case 17: return 1.02;  case 19: return 2.03;  case 20: return 1.76;
        case 26: return 1.32;  case 33: return 1.19;  case 34: return 1.20;
        case 35: return 1.20;  case 53: return 1.39;  case 3:  return 1.28;
        case 30: return 1.22;  case 78: return 1.36;  case 2:  return 0.28;
        case 10: return 0.58;  default: return 0.76;
    }
}

// Van-der-Waals radii (Bondi), Angstrom. ~1.7 default.
double vdwR(int z) {
    switch (z) {
        case 1:  return 1.20;  case 5:  return 1.92;  case 6:  return 1.70;
        case 7:  return 1.55;  case 8:  return 1.52;  case 9:  return 1.47;
        case 14: return 2.10;  case 15: return 1.80;  case 16: return 1.80;
        case 17: return 1.75;  case 35: return 1.85;  case 53: return 1.98;
        case 11: return 2.27;  case 12: return 1.73;  case 19: return 2.75;
        default: return 1.70;
    }
}

constexpr double kPi = 3.14159265358979323846;

// --------------------------------------------------------------- vec helpers
using V3 = std::array<double, 3>;
V3 vsub(const V3& a, const V3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
V3 vadd(const V3& a, const V3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
V3 vscale(const V3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
double vdot(const V3& a, const V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
double vlen(const V3& a) { return std::sqrt(std::max(0.0, vdot(a, a))); }
V3 vnorm(const V3& a) {
    const double l = vlen(a);
    return l > 1e-12 ? vscale(a, 1.0 / l) : V3{0, 0, 0};
}

// Deterministic LCG for symmetry-breaking jitter (no global RNG; reproducible).
struct Lcg {
    std::uint64_t s;
    explicit Lcg(unsigned seed) : s(seed ? seed : 0x9E3779B9u) {}
    double next() {  // (-1, 1)
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return (static_cast<double>(s >> 11) / 9007199254740992.0) * 2.0 - 1.0;
    }
};

// ----------------------------------------------------------- per-atom geometry
// Ideal VSEPR bond angle (radians) at a central atom from its hybridization.
double idealAngle(const Molecule& m, int j) {
    const Atom& a = m.atoms[j];
    bool hasTriple = false, hasDouble = false;
    for (int bi : a.bonds) {
        if (m.bonds[bi].order == 3.0) hasTriple = true;
        else if (m.bonds[bi].order == 2.0) hasDouble = true;
    }
    if (hasTriple) return kPi;                       // sp     180
    if (a.aromatic || hasDouble) return 2.0943951;   // sp2    120
    return 1.9106332;                                // sp3    109.47
}

// Ideal length (Angstrom) of one bond, from covalent radii + bond order.
double bondLength(const Molecule& m, const Bond& b) {
    const double base = covRadius(m.atoms[b.a].z) + covRadius(m.atoms[b.b].z);
    if (b.aromatic) return base * 0.92;
    if (b.order == 3.0) return base * 0.78;
    if (b.order == 2.0) return base * 0.87;
    return base;
}

// --------------------------------------------------------------- force field
struct Pair { int i, j; double d0; };

// A short steepest-descent relaxation: bonded (1-2) + angle (1-3) springs pull
// the geometry to ideal distances; soft van-der-Waals repulsion removes clashes.
void relax(std::vector<V3>& p, const std::vector<int>& zs,
           const std::vector<Pair>& bonds, const std::vector<Pair>& angles,
           const std::unordered_set<long long>& excluded) {
    const int n = static_cast<int>(p.size());
    const int kMaxIter = 600;
    const double kBond = 1.0, kAngle = 0.55, kRep = 0.6, kVdwScale = 0.62;

    auto springForce = [&](std::vector<V3>& f, const std::vector<Pair>& set, double k) {
        for (const Pair& s : set) {
            V3 d = vsub(p[s.i], p[s.j]);
            double r = vlen(d);
            if (r < 1e-6) { d = {1e-3, 0, 0}; r = 1e-3; }
            const V3 step = vscale(d, k * (r - s.d0) / r);  // toward d0
            f[s.i] = vsub(f[s.i], step);
            f[s.j] = vadd(f[s.j], step);
        }
    };

    for (int it = 0; it < kMaxIter; ++it) {
        std::vector<V3> f(n, V3{0, 0, 0});
        springForce(f, bonds, kBond);
        springForce(f, angles, kAngle);
        // Non-bonded clash repulsion (skip 1-2 / 1-3 pairs).
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) {
                if (excluded.count(static_cast<long long>(i) * n + j)) continue;
                V3 d = vsub(p[i], p[j]);
                double r = vlen(d);
                const double dmin = kVdwScale * (vdwR(zs[i]) + vdwR(zs[j]));
                if (r >= dmin) continue;
                if (r < 1e-6) { d = {1e-3, 0, 0}; r = 1e-3; }
                const V3 push = vscale(d, kRep * (dmin - r) / r);
                f[i] = vadd(f[i], push);
                f[j] = vsub(f[j], push);
            }
        const double lr = 0.10 * (1.0 - 0.5 * it / kMaxIter);
        for (int i = 0; i < n; ++i) {
            V3 mv = vscale(f[i], lr);
            const double ml = vlen(mv);
            if (ml > 0.15) mv = vscale(mv, 0.15 / ml);  // cap per-step move
            p[i] = vadd(p[i], mv);
        }
    }
}

// ----------------------------------------------------- hydrogen placement
// Place `nH` hydrogen directions around a center given the existing (fixed)
// heavy-bond directions, by relaxing points on a sphere away from each other
// (a Thomson-style spread that recovers tetrahedral / trigonal geometry).
std::vector<V3> hydrogenDirs(const std::vector<V3>& fixed, int nH, Lcg& rng) {
    std::vector<V3> nd;
    nd.reserve(nH);
    // Seed each new direction opposite the current occupied directions.
    for (int h = 0; h < nH; ++h) {
        V3 base{0, 0, 0};
        for (const V3& d : fixed) base = vadd(base, d);
        for (const V3& d : nd) base = vadd(base, d);
        base = vscale(base, -1.0);
        if (vlen(base) < 1e-3) {  // degenerate: pick a perpendicular to dir[0]
            const V3 ref = !fixed.empty() ? fixed[0] : V3{1, 0, 0};
            V3 perp = std::abs(ref[0]) < 0.9 ? V3{1, 0, 0} : V3{0, 1, 0};
            perp = vsub(perp, vscale(ref, vdot(perp, ref)));
            base = vlen(perp) > 1e-6 ? perp : V3{0, 0, 1};
        }
        base = vadd(vnorm(base), V3{rng.next() * 0.05, rng.next() * 0.05, rng.next() * 0.05});
        nd.push_back(vnorm(base));
    }
    // Relax the new directions away from fixed + each other, on the unit sphere.
    for (int it = 0; it < 80; ++it)
        for (int a = 0; a < nH; ++a) {
            V3 force{0, 0, 0};
            auto repel = [&](const V3& other) {
                V3 d = vsub(nd[a], other);
                const double r2 = std::max(1e-3, vdot(d, d));
                force = vadd(force, vscale(d, 1.0 / (r2 * std::sqrt(r2))));
            };
            for (const V3& d : fixed) repel(d);
            for (int b = 0; b < nH; ++b)
                if (b != a) repel(nd[b]);
            // Project the force onto the tangent plane, take a small step.
            force = vsub(force, vscale(nd[a], vdot(force, nd[a])));
            nd[a] = vnorm(vadd(nd[a], vscale(force, 0.10)));
        }
    return nd;
}

}  // namespace

double covalentRadius(int z) { return covRadius(z); }
double vdwRadius(int z) { return vdwR(z); }

Conformer embed3D(const Molecule& m, unsigned seed) {
    Conformer out;
    const int n = m.heavyAtomCount();
    if (n == 0) return out;

    Lcg rng(seed);
    std::vector<V3> p(n, V3{0, 0, 0});
    std::vector<int> zs(n);
    for (int i = 0; i < n; ++i) zs[i] = m.atoms[i].z;

    // --- distance matrix: bond lengths on edges, weighted shortest path else,
    //     with 1-3 entries overridden by the law-of-cosines through-space value.
    const double kInf = std::numeric_limits<double>::infinity();
    Eigen::MatrixXd D = Eigen::MatrixXd::Constant(n, n, kInf);
    for (int i = 0; i < n; ++i) D(i, i) = 0.0;
    for (const Bond& b : m.bonds) {
        const double L = bondLength(m, b);
        D(b.a, b.b) = std::min(D(b.a, b.b), L);
        D(b.b, b.a) = D(b.a, b.b);
    }
    for (int k = 0; k < n; ++k)  // Floyd-Warshall (n small)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (D(i, k) + D(k, j) < D(i, j)) D(i, j) = D(i, k) + D(k, j);

    // Replace any unreachable pairs (disconnected fragments) with a large finite.
    double maxFinite = 1.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (std::isfinite(D(i, j))) maxFinite = std::max(maxFinite, D(i, j));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (!std::isfinite(D(i, j))) D(i, j) = maxFinite * 1.6;

    // 1-3 through-space override (compact, ring-friendly start).
    std::vector<Pair> bonds, angles;
    std::unordered_set<long long> excluded;
    auto excl = [&](int i, int j) {
        const int a = std::min(i, j), b = std::max(i, j);
        excluded.insert(static_cast<long long>(a) * n + b);
    };
    for (const Bond& b : m.bonds) {
        bonds.push_back({b.a, b.b, bondLength(m, b)});
        excl(b.a, b.b);
    }
    for (int j = 0; j < n; ++j) {
        const auto& nb = m.atoms[j].nbr;
        const double theta = idealAngle(m, j);
        for (size_t x = 0; x < nb.size(); ++x)
            for (size_t y = x + 1; y < nb.size(); ++y) {
                const int i = nb[x], k = nb[y];
                const double l1 = D(i, j), l2 = D(j, k);
                const double d13 = std::sqrt(std::max(
                    0.01, l1 * l1 + l2 * l2 - 2.0 * l1 * l2 * std::cos(theta)));
                D(i, k) = d13;
                D(k, i) = d13;
                angles.push_back({i, k, d13});
                excl(i, k);
            }
    }

    // --- classical metric-matrix (MDS) embedding via Eigen ------------------
    if (n == 1) {
        p[0] = {0, 0, 0};
    } else if (n == 2) {
        p[0] = {0, 0, 0};
        p[1] = {D(0, 1), 0, 0};
    } else {
        Eigen::MatrixXd D2(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) D2(i, j) = D(i, j) * D(i, j);
        const Eigen::MatrixXd J =
            Eigen::MatrixXd::Identity(n, n) - Eigen::MatrixXd::Constant(n, n, 1.0 / n);
        const Eigen::MatrixXd B = -0.5 * J * D2 * J;  // Gram matrix
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(B);
        const Eigen::VectorXd ev = es.eigenvalues();        // ascending
        const Eigen::MatrixXd V = es.eigenvectors();
        for (int axis = 0; axis < 3; ++axis) {
            const int col = n - 1 - axis;                   // top-3 descending
            const double lam = std::max(0.0, ev(col));
            const double s = std::sqrt(lam);
            for (int i = 0; i < n; ++i) p[i][axis] = V(i, col) * s;
        }
    }

    // Symmetry-breaking jitter so degenerate/planar starts relax cleanly.
    for (int i = 0; i < n; ++i)
        for (int d = 0; d < 3; ++d) p[i][d] += rng.next() * 0.06;

    relax(p, zs, bonds, angles, excluded);

    // --- place explicit hydrogens from the relaxed heavy geometry -----------
    std::vector<V3> ap = p;                    // full atom system (heavy + H)
    std::vector<int> az = zs;
    std::vector<std::vector<int>> hOf(n);      // H indices attached to each heavy
    for (int i = 0; i < n; ++i) {
        const int h = m.atoms[i].totalH();
        if (h <= 0) continue;
        std::vector<V3> fixed;
        for (int nb : m.atoms[i].nbr) fixed.push_back(vnorm(vsub(p[nb], p[i])));
        const double chLen = covRadius(zs[i]) + covRadius(1);
        for (const V3& dir : hydrogenDirs(fixed, h, rng)) {
            hOf[i].push_back(static_cast<int>(ap.size()));
            ap.push_back(vadd(p[i], vscale(dir, chLen)));
            az.push_back(1);
        }
    }

    // --- relax the FULL system (heavy + H) so independently-placed hydrogens
    //     on neighbouring centres (e.g. tert-butyl methyls) cannot clash. -----
    const int N = static_cast<int>(ap.size());
    std::vector<Pair> fb, fa;
    std::unordered_set<long long> fex;
    auto fexcl = [&](int i, int j) {
        fex.insert(static_cast<long long>(std::min(i, j)) * N + std::max(i, j));
    };
    std::vector<std::vector<int>> adj(N);
    auto edgeLen = [&](int a, int b) -> double {
        if (a < n && b < n) {
            for (int bi : m.atoms[a].bonds)
                if (m.bonds[bi].other(a) == b) return bondLength(m, m.bonds[bi]);
            return covRadius(az[a]) + covRadius(az[b]);
        }
        const int heavy = a < n ? a : b;
        return covRadius(zs[heavy]) + covRadius(1);
    };
    auto addEdge = [&](int i, int j) {
        fb.push_back({i, j, edgeLen(i, j)});
        fexcl(i, j);
        adj[i].push_back(j);
        adj[j].push_back(i);
    };
    for (const Bond& b : m.bonds) addEdge(b.a, b.b);
    for (int i = 0; i < n; ++i)
        for (int hidx : hOf[i]) addEdge(i, hidx);
    for (int j = 0; j < n; ++j) {  // 1-3 angle springs centred on heavy atoms
        const double theta = idealAngle(m, j);
        const auto& nb = adj[j];
        for (size_t x = 0; x < nb.size(); ++x)
            for (size_t y = x + 1; y < nb.size(); ++y) {
                const int i = nb[x], k = nb[y];
                const double l1 = edgeLen(j, i), l2 = edgeLen(j, k);
                const double d13 = std::sqrt(std::max(
                    0.01, l1 * l1 + l2 * l2 - 2.0 * l1 * l2 * std::cos(theta)));
                fa.push_back({i, k, d13});
                fexcl(i, k);
            }
    }
    relax(ap, az, fb, fa, fex);

    // --- assemble final conformer ------------------------------------------
    out.heavyCount = n;
    out.pos.resize(N);
    out.z = az;
    for (int i = 0; i < N; ++i) out.pos[i] = {ap[i][0], ap[i][1], ap[i][2]};
    for (const Bond& b : m.bonds) out.bonds.emplace_back(b.a, b.b);
    for (int i = 0; i < n; ++i)
        for (int hidx : hOf[i]) out.bonds.emplace_back(i, hidx);

    // --- center at centroid + final NaN guard ------------------------------
    Vec3 c{0, 0, 0};
    for (const Vec3& v : out.pos) { c.x += v.x; c.y += v.y; c.z += v.z; }
    const double inv = 1.0 / out.size();
    c.x *= inv; c.y *= inv; c.z *= inv;
    for (int i = 0; i < out.size(); ++i) {
        Vec3& v = out.pos[i];
        v.x -= c.x; v.y -= c.y; v.z -= c.z;
        if (!std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z))
            v = {static_cast<double>(i) * 1.5, 0.0, 0.0};  // defensive fallback
    }
    return out;
}

}  // namespace biocad::chem
