#include "chem/Rings.h"

#include <algorithm>
#include <cstdint>
#include <queue>

namespace biocad::chem {

namespace {

// A GF(2) bond-incidence vector. Rings are compared and reduced as bitsets over
// bond indices, which is what makes the independence test cheap and exact.
using BitRow = std::vector<std::uint64_t>;

BitRow makeRow(std::size_t bondCount) {
    return BitRow((bondCount + 63) / 64, 0ULL);
}

void setBit(BitRow& r, int bit) { r[static_cast<std::size_t>(bit) / 64] |= 1ULL << (bit % 64); }
bool testBit(const BitRow& r, int bit) {
    return (r[static_cast<std::size_t>(bit) / 64] >> (bit % 64)) & 1ULL;
}
bool anyBit(const BitRow& r) {
    for (std::uint64_t w : r)
        if (w) return true;
    return false;
}
// Index of the lowest set bit; only called on non-empty rows.
int lowestBit(const BitRow& r) {
    for (std::size_t w = 0; w < r.size(); ++w) {
        if (!r[w]) continue;
        for (int b = 0; b < 64; ++b)
            if ((r[w] >> b) & 1ULL) return static_cast<int>(w) * 64 + b;
    }
    return -1;
}
void xorInto(BitRow& dst, const BitRow& src) {
    for (std::size_t i = 0; i < dst.size(); ++i) dst[i] ^= src[i];
}

struct Candidate {
    BitRow           bits;
    std::vector<int> atoms;  // walk order
    std::vector<int> bonds;  // walk order, bonds[i] connects atoms[i] -> atoms[i+1]
};

// Number of connected components of the atom graph. Needed for the Euler bound,
// and the reason a disconnected input like CC.CC correctly reports dim 0 rather
// than a negative number.
int componentCount(const Molecule& m) {
    const int n = static_cast<int>(m.atoms.size());
    std::vector<char> seen(static_cast<std::size_t>(n), 0);
    int comps = 0;
    std::vector<int> stack;
    for (int s = 0; s < n; ++s) {
        if (seen[static_cast<std::size_t>(s)]) continue;
        ++comps;
        stack.push_back(s);
        seen[static_cast<std::size_t>(s)] = 1;
        while (!stack.empty()) {
            const int u = stack.back();
            stack.pop_back();
            for (int v : m.atoms[static_cast<std::size_t>(u)].nbr) {
                if (seen[static_cast<std::size_t>(v)]) continue;
                seen[static_cast<std::size_t>(v)] = 1;
                stack.push_back(v);
            }
        }
    }
    return comps;
}

// Rotate/reflect a cycle to a canonical form: smallest atom index first, then
// the smaller of its two ring neighbours second. Two candidates describing the
// same ring therefore compare equal, which keeps the output deterministic.
void canonicalizeWalk(std::vector<int>& atoms, std::vector<int>& bonds) {
    const std::size_t n = atoms.size();
    std::size_t start = 0;
    for (std::size_t i = 1; i < n; ++i)
        if (atoms[i] < atoms[start]) start = i;
    const bool forward = atoms[(start + 1) % n] < atoms[(start + n - 1) % n];

    std::vector<int> a(n), b(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (forward) {
            a[i] = atoms[(start + i) % n];
            b[i] = bonds[(start + i) % n];
        } else {
            a[i] = atoms[(start + n - i) % n];
            // Walking backwards, the bond leaving a[i] is the one that entered
            // the corresponding atom in the forward walk.
            b[i] = bonds[(start + n - i - 1) % n];
        }
    }

    atoms = std::move(a);
    bonds = std::move(b);
}

}  // namespace

RingInfo perceiveRings(const Molecule& m) {
    RingInfo info;
    const int n = static_cast<int>(m.atoms.size());
    const int e = static_cast<int>(m.bonds.size());
    if (n == 0 || e == 0) return info;

    const int dim = e - n + componentCount(m);
    if (dim <= 0) return info;  // acyclic (possibly disconnected): nothing to do

    // --- BFS shortest-path trees from every atom (Horton's construction) ---
    const std::size_t un = static_cast<std::size_t>(n);
    std::vector<std::vector<int>> dist(un, std::vector<int>(un, -1));
    std::vector<std::vector<int>> parent(un, std::vector<int>(un, -1));
    std::vector<std::vector<int>> parentBond(un, std::vector<int>(un, -1));
    for (int r = 0; r < n; ++r) {
        std::queue<int> q;
        dist[static_cast<std::size_t>(r)][static_cast<std::size_t>(r)] = 0;
        q.push(r);
        while (!q.empty()) {
            const int u = q.front();
            q.pop();
            const Atom& au = m.atoms[static_cast<std::size_t>(u)];
            for (std::size_t k = 0; k < au.nbr.size(); ++k) {
                const int v = au.nbr[k];
                if (dist[static_cast<std::size_t>(r)][static_cast<std::size_t>(v)] >= 0) continue;
                dist[static_cast<std::size_t>(r)][static_cast<std::size_t>(v)] =
                    dist[static_cast<std::size_t>(r)][static_cast<std::size_t>(u)] + 1;
                parent[static_cast<std::size_t>(r)][static_cast<std::size_t>(v)] = u;
                parentBond[static_cast<std::size_t>(r)][static_cast<std::size_t>(v)] = au.bonds[k];
                q.push(v);
            }
        }
    }

    // Path from x back to root r, as (atoms, bonds) with atoms.front() == x.
    auto pathToRoot = [&](int r, int x, std::vector<int>& atoms, std::vector<int>& bonds) {
        atoms.clear();
        bonds.clear();
        int cur = x;
        while (cur != r) {
            atoms.push_back(cur);
            bonds.push_back(parentBond[static_cast<std::size_t>(r)][static_cast<std::size_t>(cur)]);
            cur = parent[static_cast<std::size_t>(r)][static_cast<std::size_t>(cur)];
        }
        atoms.push_back(r);
    };

    std::vector<Candidate> cands;
    std::vector<int> px, pxb, py, pyb;
    std::vector<char> onPath(un, 0);
    for (int r = 0; r < n; ++r) {
        for (int bi = 0; bi < e; ++bi) {
            const Bond& bd = m.bonds[static_cast<std::size_t>(bi)];
            const int x = bd.a, y = bd.b;
            if (dist[static_cast<std::size_t>(r)][static_cast<std::size_t>(x)] < 0 ||
                dist[static_cast<std::size_t>(r)][static_cast<std::size_t>(y)] < 0)
                continue;  // different component from r
            pathToRoot(r, x, px, pxb);
            pathToRoot(r, y, py, pyb);
            // The two tree paths must meet only at r, and neither may already
            // use bond bi, otherwise the "cycle" degenerates.
            bool bad = (std::find(pxb.begin(), pxb.end(), bi) != pxb.end()) ||
                       (std::find(pyb.begin(), pyb.end(), bi) != pyb.end());
            if (!bad) {
                for (int a : px) onPath[static_cast<std::size_t>(a)] = 1;
                for (int a : py) {
                    if (a == r) continue;
                    if (onPath[static_cast<std::size_t>(a)]) {
                        bad = true;
                        break;
                    }
                }
                for (int a : px) onPath[static_cast<std::size_t>(a)] = 0;
            }
            if (bad) continue;

            // Walk: x -> ... -> r -> ... -> y, then bond bi closes y -> x.
            Candidate c;
            c.atoms = px;                                              // x .. r
            for (std::size_t i = py.size(); i-- > 1;) c.atoms.push_back(py[i - 1]);
            c.bonds = pxb;                                             // x .. r bonds
            for (std::size_t i = pyb.size(); i-- > 0;) c.bonds.push_back(pyb[i]);
            c.bonds.push_back(bi);                                     // y -> x closure
            if (c.atoms.size() < 3 || c.bonds.size() != c.atoms.size()) continue;

            c.bits = makeRow(static_cast<std::size_t>(e));
            for (int b : c.bonds) setBit(c.bits, b);
            canonicalizeWalk(c.atoms, c.bonds);
            cands.push_back(std::move(c));
        }
    }

    // Smallest rings first; ties broken by the canonical atom walk so the chosen
    // SSSR is reproducible for fused systems where it is not unique.
    std::sort(cands.begin(), cands.end(), [](const Candidate& l, const Candidate& r) {
        if (l.atoms.size() != r.atoms.size()) return l.atoms.size() < r.atoms.size();
        return l.atoms < r.atoms;
    });
    cands.erase(std::unique(cands.begin(), cands.end(),
                            [](const Candidate& l, const Candidate& r) {
                                return l.bits == r.bits;
                            }),
                cands.end());

    // --- greedy GF(2) independence over bond-incidence vectors ---
    std::vector<BitRow> basis;   // reduced rows
    std::vector<int>    pivot;   // pivot bit of each reduced row
    for (const Candidate& c : cands) {
        if (static_cast<int>(info.atomRings.size()) == dim) break;
        BitRow row = c.bits;
        for (std::size_t i = 0; i < basis.size(); ++i)
            if (testBit(row, pivot[i])) xorInto(row, basis[i]);
        if (!anyBit(row)) continue;  // dependent on already-accepted rings
        const int p = lowestBit(row);
        basis.push_back(row);
        pivot.push_back(p);
        info.atomRings.push_back(c.atoms);
        info.bondRings.push_back(c.bonds);
    }
    return info;
}

void annotateRings(Molecule& m, const RingInfo& info) {
    for (Atom& a : m.atoms) a.inRing = false;
    for (Bond& b : m.bonds) b.inRing = false;
    for (const auto& ring : info.atomRings)
        for (int a : ring) m.atoms[static_cast<std::size_t>(a)].inRing = true;
    for (const auto& ring : info.bondRings)
        for (int b : ring) m.bonds[static_cast<std::size_t>(b)].inRing = true;
}

int ringSizeOf(const RingInfo& info, int atomIndex) {
    int best = 0;
    for (const auto& ring : info.atomRings) {
        if (std::find(ring.begin(), ring.end(), atomIndex) == ring.end()) continue;
        const int sz = static_cast<int>(ring.size());
        if (best == 0 || sz < best) best = sz;
    }
    return best;
}

bool inRingOfSize(const RingInfo& info, int atomIndex, int size) {
    for (const auto& ring : info.atomRings)
        if (static_cast<int>(ring.size()) == size &&
            std::find(ring.begin(), ring.end(), atomIndex) != ring.end())
            return true;
    return false;
}

int ringCountOf(const RingInfo& info, int atomIndex) {
    int c = 0;
    for (const auto& ring : info.atomRings)
        if (std::find(ring.begin(), ring.end(), atomIndex) != ring.end()) ++c;
    return c;
}

}  // namespace biocad::chem
