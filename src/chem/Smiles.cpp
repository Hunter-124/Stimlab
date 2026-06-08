#include "chem/Smiles.h"

#include <cctype>
#include <cmath>
#include <functional>
#include <map>
#include <vector>

namespace stimlab::chem {

namespace {

// Standard valences used to fill implicit hydrogens (organic subset).
std::vector<int> valenceList(int z, int charge) {
    switch (z) {
        case 1:  return {1};
        case 5:  return {3};                 // B
        case 6:  return {4};                 // C
        case 7:  return charge > 0 ? std::vector<int>{4}
                       : charge < 0 ? std::vector<int>{2}
                                    : std::vector<int>{3};  // N
        case 8:  return charge < 0 ? std::vector<int>{1}
                       : charge > 0 ? std::vector<int>{3}
                                    : std::vector<int>{2};  // O
        case 9:  return {1};                 // F
        case 14: return {4};                 // Si
        case 15: return {3, 5};              // P
        case 16: return {2, 4, 6};           // S
        case 17: return {1};                 // Cl
        case 35: return {1};                 // Br
        case 53: return {1};                 // I
        case 34: return {2, 4, 6};           // Se
        default: return {};
    }
}

void perceiveRings(Molecule& m) {
    // A bond is in a ring iff it is NOT a bridge. Find bridges via DFS low-link.
    const int n = static_cast<int>(m.atoms.size());
    std::vector<int> disc(n, -1), low(n, 0);
    int timer = 0;

    std::function<void(int, int)> dfs = [&](int u, int parentBond) {
        disc[u] = low[u] = timer++;
        for (int bi : m.atoms[u].bonds) {
            if (bi == parentBond) continue;
            const int v = m.bonds[bi].other(u);
            if (disc[v] == -1) {
                dfs(v, bi);
                low[u] = std::min(low[u], low[v]);
                if (low[v] > disc[u]) {
                    // bridge -> not a ring bond
                } else {
                    m.bonds[bi].inRing = true;
                }
            } else {
                low[u] = std::min(low[u], disc[v]);
                m.bonds[bi].inRing = true;  // back edge closes a ring
            }
        }
    };
    for (int i = 0; i < n; ++i)
        if (disc[i] == -1) dfs(i, -1);

    for (int i = 0; i < n; ++i) {
        for (int bi : m.atoms[i].bonds)
            if (m.bonds[bi].inRing) { m.atoms[i].inRing = true; break; }
    }
}

void computeImplicitH(Molecule& m) {
    for (auto& a : m.atoms) {
        if (a.bracket) continue;  // explicit H already known
        double sob = 0.0;
        for (int bi : a.bonds) sob += m.bonds[bi].order;
        const long sobR = std::lround(sob);
        const auto vals = valenceList(a.z, a.charge);
        if (vals.empty()) { a.implicitH = 0; continue; }
        int target = vals.back();
        for (int v : vals)
            if (v >= sobR) { target = v; break; }
        a.implicitH = static_cast<int>(std::max(0L, target - sobR));
    }
}

}  // namespace

std::optional<Molecule> parseSmiles(std::string_view s) {
    Molecule m;
    std::vector<int> branch;
    int prev = -1;
    double pendingOrder = 0.0;  // 0 = unspecified
    bool pendingArom = false;

    struct RingHalf { int atom; double order; bool arom; bool open; };
    std::map<int, RingHalf> rings;

    auto addAtom = [&](int z, bool arom, bool bracket, int charge, int h) {
        Atom a;
        a.z = z; a.aromatic = arom; a.bracket = bracket; a.charge = charge; a.explicitH = h;
        m.atoms.push_back(a);
        return static_cast<int>(m.atoms.size()) - 1;
    };
    auto connect = [&](int a, int b, double order, bool arom) {
        Bond bd; bd.a = a; bd.b = b;
        if (order == 0.0) {
            if (m.atoms[a].aromatic && m.atoms[b].aromatic) { bd.order = 1.5; bd.aromatic = true; }
            else bd.order = 1.0;
        } else {
            bd.order = order; bd.aromatic = arom;
        }
        const int idx = static_cast<int>(m.bonds.size());
        m.bonds.push_back(bd);
        m.atoms[a].nbr.push_back(b);  m.atoms[a].bonds.push_back(idx);
        m.atoms[b].nbr.push_back(a);  m.atoms[b].bonds.push_back(idx);
    };

    size_t i = 0;
    auto linkPrev = [&](int idx) {
        if (prev >= 0) connect(prev, idx, pendingOrder, pendingArom);
        prev = idx; pendingOrder = 0.0; pendingArom = false;
    };

    while (i < s.size()) {
        const char c = s[i];
        if (c == '(') { branch.push_back(prev); ++i; continue; }
        if (c == ')') { if (!branch.empty()) { prev = branch.back(); branch.pop_back(); } ++i; continue; }
        if (c == '-') { pendingOrder = 1.0; ++i; continue; }
        if (c == '=') { pendingOrder = 2.0; ++i; continue; }
        if (c == '#') { pendingOrder = 3.0; ++i; continue; }
        if (c == ':') { pendingOrder = 1.5; pendingArom = true; ++i; continue; }
        if (c == '/' || c == '\\') { pendingOrder = 1.0; ++i; continue; }  // ignore E/Z
        if (c == '.') { prev = -1; pendingOrder = 0.0; pendingArom = false; ++i; continue; }

        if (std::isdigit(static_cast<unsigned char>(c)) || c == '%') {
            int rn = 0;
            if (c == '%') {
                if (i + 2 >= s.size()) return std::nullopt;
                rn = (s[i + 1] - '0') * 10 + (s[i + 2] - '0'); i += 3;
            } else { rn = c - '0'; ++i; }
            auto it = rings.find(rn);
            if (it == rings.end() || !it->second.open) {
                rings[rn] = RingHalf{prev, pendingOrder, pendingArom, true};
            } else {
                const double order = pendingOrder != 0.0 ? pendingOrder : it->second.order;
                connect(it->second.atom, prev, order, pendingArom || it->second.arom);
                it->second.open = false;
            }
            pendingOrder = 0.0; pendingArom = false;
            continue;
        }

        if (c == '[') {
            const size_t j = s.find(']', i);
            if (j == std::string_view::npos) return std::nullopt;
            const std::string_view tok = s.substr(i + 1, j - i - 1);
            size_t k = 0;
            while (k < tok.size() && std::isdigit(static_cast<unsigned char>(tok[k]))) ++k;  // isotope
            if (k >= tok.size()) return std::nullopt;

            bool arom = false;
            int z = 0;
            // Element symbol: try two-letter, then one-letter (upper or aromatic lower).
            std::string sym1(1, tok[k]);
            std::string sym2 = (k + 1 < tok.size()) ? std::string{tok[k], tok[k + 1]} : "";
            if (!sym2.empty() && findElementBySymbol(sym2) &&
                std::isupper(static_cast<unsigned char>(tok[k + 1])) == 0 &&
                std::isupper(static_cast<unsigned char>(tok[k]))) {
                // two-letter like Cl, Br, Si, Na, Se
                z = findElementBySymbol(sym2)->z; k += 2;
            } else if (const ElementInfo* e1 = findElementBySymbol(sym1)) {
                z = e1->z; ++k;
            } else if (std::islower(static_cast<unsigned char>(tok[k]))) {
                std::string up(1, static_cast<char>(std::toupper(tok[k])));
                const ElementInfo* ea = findElementBySymbol(up);
                if (!ea) return std::nullopt;
                z = ea->z; arom = true; ++k;
            } else {
                return std::nullopt;
            }

            while (k < tok.size() && tok[k] == '@') ++k;  // chirality (ignored)

            int hcount = 0;
            if (k < tok.size() && tok[k] == 'H') {
                ++k; hcount = 1;
                if (k < tok.size() && std::isdigit(static_cast<unsigned char>(tok[k]))) {
                    hcount = tok[k] - '0'; ++k;
                }
            }
            int charge = 0;
            while (k < tok.size() && (tok[k] == '+' || tok[k] == '-')) {
                const int sign = (tok[k] == '+') ? 1 : -1; ++k;
                if (k < tok.size() && std::isdigit(static_cast<unsigned char>(tok[k]))) {
                    charge += sign * (tok[k] - '0'); ++k;
                } else {
                    charge += sign;
                }
            }
            const int idx = addAtom(z, arom, true, charge, hcount);
            linkPrev(idx);
            i = j + 1;
            continue;
        }

        // Organic-subset atom (outside brackets).
        int z = 0;
        bool arom = false;
        if (c == 'C' && i + 1 < s.size() && s[i + 1] == 'l') { z = 17; i += 2; }
        else if (c == 'B' && i + 1 < s.size() && s[i + 1] == 'r') { z = 35; i += 2; }
        else {
            switch (c) {
                case 'B': z = 5; break;  case 'C': z = 6; break;  case 'N': z = 7; break;
                case 'O': z = 8; break;  case 'P': z = 15; break; case 'S': z = 16; break;
                case 'F': z = 9; break;  case 'I': z = 53; break;
                case 'b': z = 5; arom = true; break;  case 'c': z = 6; arom = true; break;
                case 'n': z = 7; arom = true; break;  case 'o': z = 8; arom = true; break;
                case 'p': z = 15; arom = true; break; case 's': z = 16; arom = true; break;
                default: return std::nullopt;
            }
            ++i;
        }
        const int idx = addAtom(z, arom, false, 0, 0);
        linkPrev(idx);
    }

    if (m.atoms.empty()) return std::nullopt;
    perceiveRings(m);
    computeImplicitH(m);
    return m;
}

}  // namespace stimlab::chem
