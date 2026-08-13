#include "chem/Canonical.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace biocad::chem {

namespace {

std::uint64_t mix(std::uint64_t h, std::uint64_t v) {
    h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return h;
}

// Bond flavours as small integers so refinement codes and written bond symbols
// agree on what "the same bond" means.
int bondCode(const Bond& b) {
    if (b.aromatic) return 1;
    if (b.order == 2.0) return 3;
    if (b.order == 3.0) return 4;
    if (b.order == 1.0) return 2;
    return 5;
}

// Packed initial invariant. Field widths are generous enough that no field can
// bleed into its neighbour for any real molecule, which is what makes the packing
// order-independent rather than merely usually-distinct.
std::uint64_t initialInvariant(const Atom& a) {
    std::uint64_t v = 0;
    v = (v << 8) | static_cast<std::uint64_t>(a.z & 0xff);
    v = (v << 4) | static_cast<std::uint64_t>(std::min(a.degree(), 15));
    v = (v << 4) | static_cast<std::uint64_t>(std::clamp(a.charge, -7, 7) + 8);
    v = (v << 4) | static_cast<std::uint64_t>(std::min(std::max(a.totalH(), 0), 15));
    v = (v << 1) | static_cast<std::uint64_t>(a.aromatic ? 1 : 0);
    v = (v << 1) | static_cast<std::uint64_t>(a.inRing ? 1 : 0);
    return v;
}

// Dense ranks (0..k-1) of `code`, ties sharing a rank. Ranks - not the raw hash
// values - are what feed the next refinement round, so the result depends only on
// the ORDER of the codes and never on the magic constants inside mix().
std::vector<int> denseRanks(const std::vector<std::uint64_t>& code) {
    std::vector<std::uint64_t> sorted = code;
    std::sort(sorted.begin(), sorted.end());
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
    std::vector<int> rank(code.size(), 0);
    for (std::size_t i = 0; i < code.size(); ++i) {
        rank[i] = static_cast<int>(
            std::lower_bound(sorted.begin(), sorted.end(), code[i]) - sorted.begin());
    }
    return rank;
}

// One Morgan refinement round: hash each atom's own rank together with the sorted
// multiset of (neighbour rank, bond flavour). Sorting is essential - it is what
// makes the code independent of the order neighbours happen to sit in nbr[].
std::vector<std::uint64_t> refineOnce(const Molecule& m, const std::vector<int>& rank) {
    const std::size_t n = m.atoms.size();
    std::vector<std::uint64_t> out(n, 0);
    std::vector<std::uint64_t> desc;
    for (std::size_t i = 0; i < n; ++i) {
        desc.clear();
        for (int bi : m.atoms[i].bonds) {
            const Bond& b = m.bonds[bi];
            const int v = b.other(static_cast<int>(i));
            desc.push_back(static_cast<std::uint64_t>(rank[v]) * 8u +
                           static_cast<std::uint64_t>(bondCode(b)));
        }
        std::sort(desc.begin(), desc.end());
        std::uint64_t h = mix(0x243f6a8885a308d3ULL, static_cast<std::uint64_t>(rank[i]));
        for (std::uint64_t d : desc) h = mix(h, d);
        out[i] = h;
    }
    return out;
}

// Refine to a stable partition (the rank vector stops moving).
std::vector<int> refineToFixpoint(const Molecule& m, std::vector<std::uint64_t> code) {
    std::vector<int> rank = denseRanks(code);
    for (std::size_t iter = 0; iter <= m.atoms.size(); ++iter) {
        code = refineOnce(m, rank);
        std::vector<int> next = denseRanks(code);
        if (next == rank) break;
        rank = std::move(next);
    }
    return rank;
}

bool isOrganicSubset(int z) {
    switch (z) {
        case 5: case 6: case 7: case 8: case 9:
        case 15: case 16: case 17: case 35: case 53: return true;
        default: return false;
    }
}

// Mirrors the implicit-H rule inside parseSmiles(). It is duplicated on purpose:
// the writer must predict EXACTLY what the parser will infer, otherwise it would
// drop an H that the round trip cannot recover. Keep the two in step.
std::vector<int> valenceList(int z, int charge) {
    switch (z) {
        case 1:  return {1};
        case 5:  return {3};
        case 6:  return {4};
        case 7:  return charge > 0 ? std::vector<int>{4}
                       : charge < 0 ? std::vector<int>{2}
                                    : std::vector<int>{3};
        case 8:  return charge < 0 ? std::vector<int>{1}
                       : charge > 0 ? std::vector<int>{3}
                                    : std::vector<int>{2};
        case 9:  return {1};
        case 14: return {4};
        case 15: return {3, 5};
        case 16: return {2, 4, 6};
        case 17: return {1};
        case 35: return {1};
        case 53: return {1};
        case 34: return {2, 4, 6};
        default: return {};
    }
}

int impliedH(const Molecule& m, int idx) {
    const Atom& a = m.atoms[idx];
    double sob = 0.0;
    for (int bi : a.bonds) sob += m.bonds[bi].order;
    const long sobR = std::lround(sob);
    const auto vals = valenceList(a.z, a.charge);
    if (vals.empty()) return 0;
    int target = vals.back();
    for (int v : vals)
        if (v >= sobR) { target = v; break; }
    return static_cast<int>(std::max(0L, static_cast<long>(target) - sobR));
}

std::string elementSymbol(const Atom& a) {
    const char* s = symbolByZ(a.z);
    std::string sym = s ? s : "*";
    if (a.aromatic && !sym.empty())
        sym[0] = static_cast<char>(std::tolower(static_cast<unsigned char>(sym[0])));
    return sym;
}

// Bracket only where the bare symbol would lose information: a charge, an atom
// outside the organic subset, or an H count the parser would not re-infer (the
// classic case being aromatic pyrrole-type [nH]).
std::string atomToken(const Molecule& m, int idx) {
    const Atom& a = m.atoms[idx];
    const std::string sym = elementSymbol(a);
    const bool plain = a.charge == 0 && isOrganicSubset(a.z) && a.totalH() == impliedH(m, idx);
    if (plain) return sym;

    std::string t = "[";
    t += sym;
    const int h = a.totalH();
    if (h == 1) t += "H";
    else if (h > 1) t += "H" + std::to_string(h);
    if (a.charge != 0) {
        t += (a.charge > 0) ? '+' : '-';
        const int mag = a.charge > 0 ? a.charge : -a.charge;
        if (mag > 1) t += std::to_string(mag);
    }
    t += "]";
    return t;
}

// '-' is suppressed except between two aromatic atoms, where a bare bond would be
// read back as part of an aromatic system (the biphenyl link).
std::string bondToken(const Molecule& m, int bi) {
    const Bond& b = m.bonds[bi];
    if (b.aromatic) {
        const bool bothArom = m.atoms[b.a].aromatic && m.atoms[b.b].aromatic;
        return bothArom ? "" : ":";
    }
    if (b.order == 2.0) return "=";
    if (b.order == 3.0) return "#";
    if (m.atoms[b.a].aromatic && m.atoms[b.b].aromatic) return "-";
    return "";
}

std::string ringDigit(int d) {
    if (d < 10) return std::string(1, static_cast<char>('0' + d));
    return "%" + std::string(1, static_cast<char>('0' + d / 10)) +
           std::string(1, static_cast<char>('0' + d % 10));
}

// Neighbours of `u` as (bond index, neighbour), ordered by canonical rank then by
// bond index so parallel edges are still deterministic.
std::vector<std::pair<int, int>> orderedNeighbors(const Molecule& m,
                                                  const std::vector<int>& rank, int u) {
    std::vector<std::pair<int, int>> out;
    out.reserve(m.atoms[u].bonds.size());
    for (int bi : m.atoms[u].bonds) out.emplace_back(bi, m.bonds[bi].other(u));
    std::sort(out.begin(), out.end(), [&](const auto& x, const auto& y) {
        if (rank[x.second] != rank[y.second]) return rank[x.second] < rank[y.second];
        return x.first < y.first;
    });
    return out;
}

// Writes one connected component, rooted at `root`.
std::string writeComponent(const Molecule& m, const std::vector<int>& rank, int root) {
    const std::size_t nb = m.bonds.size();
    std::vector<char> visited(m.atoms.size(), 0);
    std::vector<char> isRingClosure(nb, 0), classified(nb, 0);

    // Pass 1: classify bonds as spanning-tree edges or ring closures, walking in
    // exactly the order pass 2 will use, so the two passes agree.
    std::function<void(int, int)> classify = [&](int u, int parentBond) {
        visited[u] = 1;
        for (const auto& [bi, v] : orderedNeighbors(m, rank, u)) {
            if (bi == parentBond || classified[bi]) continue;
            if (!visited[v]) {
                classified[bi] = 1;
                classify(v, bi);
            } else {
                classified[bi] = 1;
                isRingClosure[bi] = 1;
            }
        }
    };
    classify(root, -1);

    // Pass 2: emit. Digits are allocated on first need and released on closure, so
    // a large molecule reuses low digits instead of overflowing past 9.
    std::fill(visited.begin(), visited.end(), 0);
    std::vector<int> digitOf(nb, -1);
    std::vector<char> digitInUse(100, 0);
    std::string out;

    std::function<void(int, int)> emit = [&](int u, int parentBond) {
        visited[u] = 1;
        out += atomToken(m, u);

        const auto nbrs = orderedNeighbors(m, rank, u);
        for (const auto& [bi, v] : nbrs) {
            if (!isRingClosure[bi]) continue;
            if (digitOf[bi] < 0) {
                int d = 1;
                while (d < 100 && digitInUse[d]) ++d;
                digitOf[bi] = d;
                digitInUse[d] = 1;
                out += bondToken(m, bi) + ringDigit(d);
            } else {
                out += bondToken(m, bi) + ringDigit(digitOf[bi]);
                digitInUse[digitOf[bi]] = 0;
            }
        }

        std::vector<std::pair<int, int>> children;
        for (const auto& [bi, v] : nbrs) {
            if (bi == parentBond || isRingClosure[bi] || visited[v]) continue;
            children.emplace_back(bi, v);
        }
        for (std::size_t k = 0; k < children.size(); ++k) {
            const auto [bi, v] = children[k];
            if (visited[v]) continue;  // reached through an earlier branch
            const bool last = (k + 1 == children.size());
            if (!last) out += "(";
            out += bondToken(m, bi);
            emit(v, bi);
            if (!last) out += ")";
        }
    };
    emit(root, -1);
    return out;
}

}  // namespace

std::vector<int> canonicalRanks(const Molecule& m) {
    const std::size_t n = m.atoms.size();
    if (n == 0) return {};

    std::vector<std::uint64_t> code(n);
    for (std::size_t i = 0; i < n; ++i) code[i] = initialInvariant(m.atoms[i]);
    std::vector<int> rank = refineToFixpoint(m, std::move(code));

    // Break ties one at a time: lowest tied class first, lowest atom index inside
    // it (atoms in a class already share their invariant tuple, so the tuple
    // comparison is already satisfied), then refine again. Refining after EACH
    // individual break is what turns a stable partition into a total order rather
    // than merely a finer partition.
    for (std::size_t guard = 0; guard <= n; ++guard) {
        std::map<int, std::vector<int>> classes;
        for (std::size_t i = 0; i < n; ++i) classes[rank[i]].push_back(static_cast<int>(i));

        int chosen = -1;
        for (const auto& [r, members] : classes) {
            if (members.size() > 1) { chosen = members.front(); break; }
        }
        if (chosen < 0) break;  // already a total order

        std::vector<std::uint64_t> seeded(n);
        for (std::size_t i = 0; i < n; ++i)
            seeded[i] = static_cast<std::uint64_t>(rank[i]) * 2u + 1u;
        seeded[static_cast<std::size_t>(chosen)] -= 1u;  // strictly ahead of its class
        rank = refineToFixpoint(m, std::move(seeded));
    }
    return rank;
}

std::string canonicalSmiles(const Molecule& m) {
    if (m.empty()) return {};
    const std::vector<int> rank = canonicalRanks(m);

    // Connected components; each is rooted at its own lowest-ranked atom so the
    // walk never depends on input atom order.
    const int n = static_cast<int>(m.atoms.size());
    std::vector<int> comp(n, -1);
    int ncomp = 0;
    for (int i = 0; i < n; ++i) {
        if (comp[i] >= 0) continue;
        std::vector<int> stack{i};
        comp[i] = ncomp;
        while (!stack.empty()) {
            const int u = stack.back();
            stack.pop_back();
            for (int v : m.atoms[u].nbr)
                if (comp[v] < 0) { comp[v] = ncomp; stack.push_back(v); }
        }
        ++ncomp;
    }

    std::vector<int> root(static_cast<std::size_t>(ncomp), -1);
    for (int i = 0; i < n; ++i) {
        int& r = root[static_cast<std::size_t>(comp[i])];
        if (r < 0 || rank[i] < rank[r]) r = i;
    }

    std::vector<std::string> parts;
    parts.reserve(static_cast<std::size_t>(ncomp));
    for (int c = 0; c < ncomp; ++c) parts.push_back(writeComponent(m, rank, root[c]));
    // Ordering components by their own canonical string is what makes "CC.c1ccccc1"
    // and "c1ccccc1.CC" the same molecule to this writer.
    std::sort(parts.begin(), parts.end());

    std::string out;
    for (std::size_t i = 0; i < parts.size(); ++i) {
        if (i) out += ".";
        out += parts[i];
    }
    return out;
}

std::uint64_t graphHash(const Molecule& m) {
    std::uint64_t h = mix(0xcbf29ce484222325ULL, m.atoms.size());
    h = mix(h, m.bonds.size());
    const std::string s = canonicalSmiles(m);
    for (unsigned char c : s) h = mix(h, c);
    return h;
}

}  // namespace biocad::chem
