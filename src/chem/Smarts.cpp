#include "chem/Smarts.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <set>
#include <unordered_map>

#include "chem/Aromaticity.h"

namespace biocad::chem {

namespace {

// Every recursive subpattern gets a process-unique id so the memo key
// (id, atom) can never collide between two different $() subpatterns.
std::int64_t nextPatternId() {
    static std::atomic<std::int64_t> counter{1};
    return counter.fetch_add(1);
}

// ---- parser --------------------------------------------------------------

class Parser {
public:
    Parser(std::string_view s) : s_(s) {}

    bool run(SmartsPattern& out) {
        pat_.source.assign(s_);
        pat_.id = nextPatternId();
        if (!parseChain(-1)) return false;
        if (i_ < s_.size()) return fail("unexpected character");
        if (!openRings_.empty()) return fail("unclosed ring bond");
        if (pat_.atoms.empty()) return fail("empty pattern");
        out = std::move(pat_);
        return true;
    }

    const std::string& error() const { return err_; }

private:
    // ---- low-level helpers
    bool eof() const { return i_ >= s_.size(); }
    char peek(std::size_t off = 0) const { return i_ + off < s_.size() ? s_[i_ + off] : '\0'; }

    bool fail(const char* what) {
        if (err_.empty())
            err_ = std::string(what) + " at offset " + std::to_string(i_) + " in \"" +
                   std::string(s_) + "\"";
        return false;
    }

    int addAtomNode(AtomNode n) {
        pat_.atomNodes.push_back(n);
        return static_cast<int>(pat_.atomNodes.size()) - 1;
    }
    int addBondNode(BondNode n) {
        pat_.bondNodes.push_back(n);
        return static_cast<int>(pat_.bondNodes.size()) - 1;
    }

    // ---- chain / branch structure
    // `prevAtom` is the atom the next atom bonds back to (-1 at a component start).
    bool parseChain(int prevAtom) {
        int prev = prevAtom;
        int pendingBond = -1;   // parsed bond expression awaiting its right-hand atom
        bool pendingDot = false;

        while (!eof()) {
            const char c = peek();
            if (c == ')') break;

            if (c == '.') {
                ++i_;
                pendingDot = true;
                prev = -1;
                continue;
            }
            if (c == '(') {
                if (prev < 0) return fail("branch with no preceding atom");
                ++i_;
                if (!parseChain(prev)) return false;
                if (peek() != ')') return fail("unclosed branch");
                ++i_;
                continue;
            }
            if (isBondStart(c)) {
                if (pendingBond >= 0) return fail("two consecutive bond expressions");
                pendingBond = parseBondExpr();
                if (pendingBond < 0) return false;
                continue;
            }
            if (std::isdigit(static_cast<unsigned char>(c)) || c == '%') {
                if (prev < 0) return fail("ring closure with no preceding atom");
                if (!parseRingClosure(prev, pendingBond)) return false;
                pendingBond = -1;
                continue;
            }

            const int atom = parseAtom();
            if (atom < 0) return false;
            if (prev >= 0) {
                const int expr = pendingBond >= 0
                                     ? pendingBond
                                     : addBondNode(BondNode{ExprOp::Prim, BondPrimKind::Default,
                                                            -1, -1});
                addBond(prev, atom, expr);
            } else if (pendingBond >= 0 && !pendingDot) {
                return fail("bond expression with no preceding atom");
            }
            pendingBond = -1;
            pendingDot = false;
            prev = atom;
        }
        if (pendingBond >= 0) return fail("dangling bond expression");
        return true;
    }

    void addBond(int a, int b, int expr) {
        const int bi = static_cast<int>(pat_.bonds.size());
        pat_.bonds.push_back(PatternBond{a, b, expr});
        pat_.atoms[static_cast<std::size_t>(a)].bonds.push_back(bi);
        pat_.atoms[static_cast<std::size_t>(b)].bonds.push_back(bi);
    }

    bool parseRingClosure(int atom, int bondExpr) {
        int label = 0;
        if (peek() == '%') {
            ++i_;
            if (!std::isdigit(static_cast<unsigned char>(peek())) ||
                !std::isdigit(static_cast<unsigned char>(peek(1))))
                return fail("%% ring label needs two digits");
            label = (peek() - '0') * 10 + (peek(1) - '0');
            i_ += 2;
        } else {
            label = peek() - '0';
            ++i_;
        }
        auto it = openRings_.find(label);
        if (it == openRings_.end()) {
            openRings_.emplace(label, OpenRing{atom, bondExpr});
            return true;
        }
        const OpenRing open = it->second;
        openRings_.erase(it);
        if (open.atom == atom) return fail("ring bond to the same atom");
        int expr = open.bondExpr >= 0 ? open.bondExpr : bondExpr;
        if (expr < 0)
            expr = addBondNode(BondNode{ExprOp::Prim, BondPrimKind::Default, -1, -1});
        addBond(open.atom, atom, expr);
        return true;
    }

    // ---- bonds
    static bool isBondStart(char c) {
        return c == '-' || c == '=' || c == '#' || c == ':' || c == '~' || c == '@' || c == '/' ||
               c == '\\' || c == '!';
    }

    // Same precedence ladder as atoms: ! > & > , > ;
    int parseBondExpr() { return parseBondLow(); }

    int parseBondLow() {
        int lhs = parseBondOr();
        if (lhs < 0) return -1;
        while (peek() == ';') {
            ++i_;
            const int rhs = parseBondOr();
            if (rhs < 0) return -1;
            lhs = addBondNode(BondNode{ExprOp::And, BondPrimKind::Default, lhs, rhs});
        }
        return lhs;
    }
    int parseBondOr() {
        int lhs = parseBondAnd();
        if (lhs < 0) return -1;
        while (peek() == ',') {
            ++i_;
            const int rhs = parseBondAnd();
            if (rhs < 0) return -1;
            lhs = addBondNode(BondNode{ExprOp::Or, BondPrimKind::Default, lhs, rhs});
        }
        return lhs;
    }
    int parseBondAnd() {
        int lhs = parseBondUnary();
        if (lhs < 0) return -1;
        for (;;) {
            if (peek() == '&') {
                ++i_;
            } else if (!isBondStart(peek())) {
                break;  // implicit high-precedence and
            }
            const int rhs = parseBondUnary();
            if (rhs < 0) return -1;
            lhs = addBondNode(BondNode{ExprOp::And, BondPrimKind::Default, lhs, rhs});
        }
        return lhs;
    }
    int parseBondUnary() {
        if (peek() == '!') {
            ++i_;
            const int sub = parseBondUnary();
            if (sub < 0) return -1;
            return addBondNode(BondNode{ExprOp::Not, BondPrimKind::Default, sub, -1});
        }
        BondPrimKind k;
        switch (peek()) {
            case '-': k = BondPrimKind::Single; break;
            case '=': k = BondPrimKind::Double; break;
            case '#': k = BondPrimKind::Triple; break;
            case ':': k = BondPrimKind::Aromatic; break;
            case '~': k = BondPrimKind::Any; break;
            case '@': k = BondPrimKind::RingBond; break;
            // Directional bonds carry stereo the graph does not store; they are
            // plain single bonds here, which is what every toolkit does.
            case '/':
            case '\\': k = BondPrimKind::Single; break;
            default: fail("expected a bond primitive"); return -1;
        }
        ++i_;
        return addBondNode(BondNode{ExprOp::Prim, k, -1, -1});
    }

    // ---- atoms
    int parseAtom() {
        PatternAtom pa;
        if (peek() == '[') {
            ++i_;
            const int expr = parseAtomLow();
            if (expr < 0) return -1;
            if (peek() == ':') {  // atom map, retained for a later SMIRKS phase
                ++i_;
                if (!std::isdigit(static_cast<unsigned char>(peek()))) {
                    fail("atom map needs a number");
                    return -1;
                }
                int v = 0;
                while (std::isdigit(static_cast<unsigned char>(peek()))) v = v * 10 + (s_[i_++] - '0');
                pa.map = v;
            }
            if (peek() != ']') {
                fail("unclosed bracket atom");
                return -1;
            }
            ++i_;
            pa.expr = expr;
        } else {
            const int expr = parseOrganicAtom();
            if (expr < 0) return -1;
            pa.expr = expr;
        }
        pat_.atoms.push_back(std::move(pa));
        return static_cast<int>(pat_.atoms.size()) - 1;
    }

    int parseOrganicAtom() {
        if (peek() == '*') {
            ++i_;
            return addAtomNode(AtomNode{ExprOp::Prim, AtomPrimKind::Any, 0, -1, -1, -1, -1});
        }
        // a and A are legal unbracketed SMARTS primitives, not just [a] / [A].
        if (peek() == 'a') {
            ++i_;
            return addAtomNode(
                AtomNode{ExprOp::Prim, AtomPrimKind::AromaticAny, 0, -1, -1, -1, -1});
        }
        if (peek() == 'A') {
            ++i_;
            return addAtomNode(
                AtomNode{ExprOp::Prim, AtomPrimKind::AliphaticAny, 0, -1, -1, -1, -1});
        }
        if (peek() == '$') {
            fail("recursive SMARTS $() must be written inside brackets");
            return -1;
        }
        // Organic subset only, outside brackets, as SMARTS requires.
        static const char* kTwo[] = {"Cl", "Br"};
        for (const char* t : kTwo) {
            if (peek() == t[0] && peek(1) == t[1]) {
                i_ += 2;
                const ElementInfo* e = findElementBySymbol(t);
                return addAtomNode(
                    AtomNode{ExprOp::Prim, AtomPrimKind::Element, e->z, 0, -1, -1, -1});
            }
        }
        const char c = peek();
        std::string sym(1, c);
        int arom = std::islower(static_cast<unsigned char>(c)) ? 1 : 0;
        if (arom == 1) sym[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        static const std::string kAliphatic = "BCNOPSFI";
        static const std::string kAromatic = "bcnops";
        if ((arom == 0 && kAliphatic.find(c) == std::string::npos) ||
            (arom == 1 && kAromatic.find(c) == std::string::npos)) {
            fail("only the organic subset may appear outside brackets");
            return -1;
        }
        const ElementInfo* e = findElementBySymbol(sym);
        if (e == nullptr) {
            fail("unknown element");
            return -1;
        }
        ++i_;
        return addAtomNode(AtomNode{ExprOp::Prim, AtomPrimKind::Element, e->z, arom, -1, -1, -1});
    }

    int parseAtomLow() {
        int lhs = parseAtomOr();
        if (lhs < 0) return -1;
        while (peek() == ';') {
            ++i_;
            const int rhs = parseAtomOr();
            if (rhs < 0) return -1;
            lhs = addAtomNode(AtomNode{ExprOp::And, AtomPrimKind::Any, 0, -1, -1, lhs, rhs});
        }
        return lhs;
    }
    int parseAtomOr() {
        int lhs = parseAtomAnd();
        if (lhs < 0) return -1;
        while (peek() == ',') {
            ++i_;
            const int rhs = parseAtomAnd();
            if (rhs < 0) return -1;
            lhs = addAtomNode(AtomNode{ExprOp::Or, AtomPrimKind::Any, 0, -1, -1, lhs, rhs});
        }
        return lhs;
    }
    int parseAtomAnd() {
        int lhs = parseAtomUnary();
        if (lhs < 0) return -1;
        for (;;) {
            if (peek() == '&') {
                ++i_;
            } else if (!startsAtomPrim(peek())) {
                break;
            }
            const int rhs = parseAtomUnary();
            if (rhs < 0) return -1;
            lhs = addAtomNode(AtomNode{ExprOp::And, AtomPrimKind::Any, 0, -1, -1, lhs, rhs});
        }
        return lhs;
    }
    static bool startsAtomPrim(char c) {
        if (c == '\0' || c == ']' || c == ';' || c == ',' || c == '&' || c == ':') return false;
        return c == '!' || c == '*' || c == '#' || c == '$' || c == '+' || c == '-' || c == '@' ||
               c == '^' || std::isalnum(static_cast<unsigned char>(c));
    }
    int parseAtomUnary() {
        if (peek() == '!') {
            ++i_;
            const int sub = parseAtomUnary();
            if (sub < 0) return -1;
            return addAtomNode(AtomNode{ExprOp::Not, AtomPrimKind::Any, 0, -1, -1, sub, -1});
        }
        return parseAtomPrim();
    }

    int readNumber(int fallback) {
        if (!std::isdigit(static_cast<unsigned char>(peek()))) return fallback;
        int v = 0;
        while (std::isdigit(static_cast<unsigned char>(peek()))) v = v * 10 + (s_[i_++] - '0');
        return v;
    }

    int prim(AtomPrimKind k, int value, int arom = -1, int sub = -1) {
        return addAtomNode(AtomNode{ExprOp::Prim, k, value, arom, sub, -1, -1});
    }

    int parseAtomPrim() {
        const char c = peek();
        if (c == '\0' || c == ']') {
            fail("expected an atom primitive");
            return -1;
        }
        if (std::isdigit(static_cast<unsigned char>(c))) {
            fail("isotope specifications are not supported (the graph stores no isotopes)");
            return -1;
        }
        if (c == '^') {
            fail("hybridisation primitive ^n is not supported");
            return -1;
        }
        if (c == '*') {
            ++i_;
            return prim(AtomPrimKind::Any, 0);
        }
        if (c == '#') {
            ++i_;
            if (!std::isdigit(static_cast<unsigned char>(peek()))) {
                fail("#n needs an atomic number");
                return -1;
            }
            return prim(AtomPrimKind::Element, readNumber(0), -1);
        }
        if (c == '$') {
            ++i_;
            if (peek() != '(') {
                fail("recursive SMARTS needs $(");
                return -1;
            }
            ++i_;
            const std::size_t start = i_;
            int depth = 1;
            while (!eof() && depth > 0) {
                if (peek() == '(') ++depth;
                if (peek() == ')') --depth;
                ++i_;
            }
            if (depth != 0) {
                fail("unclosed recursive SMARTS");
                return -1;
            }
            const std::string_view inner = s_.substr(start, i_ - start - 1);
            Parser sub(inner);
            auto child = std::make_shared<SmartsPattern>();
            if (!sub.run(*child)) {
                err_ = "inside $(): " + sub.error();
                return -1;
            }
            pat_.recursives.push_back(std::move(child));
            return prim(AtomPrimKind::Recursive, 0, -1,
                        static_cast<int>(pat_.recursives.size()) - 1);
        }
        if (c == '+' || c == '-') {
            ++i_;
            int sign = (c == '+') ? 1 : -1;
            int n = 1;
            while (peek() == c) {  // ++ / -- / +++
                ++i_;
                ++n;
            }
            if (n == 1) n = readNumber(1);
            return prim(AtomPrimKind::Charge, sign * n);
        }
        if (c == '@') {  // chirality: accepted and ignored, see header
            ++i_;
            if (peek() == '@') ++i_;
            return prim(AtomPrimKind::Any, 0);
        }

        // Two-letter element symbols win over single-letter primitives, so
        // [Br] is bromine and [Hg] is mercury rather than B&r or H&g.
        if (std::isalpha(static_cast<unsigned char>(c)) &&
            std::islower(static_cast<unsigned char>(peek(1)))) {
            std::string sym{c, peek(1)};
            int arom = 0;
            if (std::islower(static_cast<unsigned char>(c))) {
                sym[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                arom = 1;  // aromatic se / as
            }
            const ElementInfo* e = findElementBySymbol(sym);
            if (e != nullptr) {
                i_ += 2;
                return prim(AtomPrimKind::Element, e->z, arom);
            }
        }
        switch (c) {
            case 'a': ++i_; return prim(AtomPrimKind::AromaticAny, 0);
            case 'A': ++i_; return prim(AtomPrimKind::AliphaticAny, 0);
            case 'H': ++i_; return prim(AtomPrimKind::TotalH, readNumber(1));
            case 'X': ++i_; return prim(AtomPrimKind::Connectivity, readNumber(1));
            case 'D': ++i_; return prim(AtomPrimKind::Degree, readNumber(1));
            // R with no number means "in at least one ring": value < 0.
            case 'R': ++i_; return prim(AtomPrimKind::RingCount, readNumber(-1));
            case 'r': {
                ++i_;
                const int n = readNumber(-1);
                if (n < 0) return prim(AtomPrimKind::InAnyRing, 0);
                if (n < 3) {
                    fail("ring size rn must be at least 3");
                    return -1;
                }
                return prim(AtomPrimKind::RingSize, n);
            }
            default: break;
        }
        if (std::isalpha(static_cast<unsigned char>(c))) {
            std::string sym(1, c);
            int arom = 0;
            if (std::islower(static_cast<unsigned char>(c))) {
                sym[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                arom = 1;
            }
            const ElementInfo* e = findElementBySymbol(sym);
            if (e != nullptr) {
                ++i_;
                return prim(AtomPrimKind::Element, e->z, arom);
            }
        }
        fail("unrecognised atom primitive");
        return -1;
    }

    struct OpenRing {
        int atom;
        int bondExpr;
    };

    std::string_view s_;
    std::size_t i_ = 0;
    std::string err_;
    SmartsPattern pat_;
    std::unordered_map<int, OpenRing> openRings_;
};

// ---- matcher -------------------------------------------------------------

class Matcher {
public:
    Matcher(const Molecule& mol, const RingInfo& rings) : mol_(mol), rings_(rings) {
        const std::size_t n = mol.atoms.size();
        // Cached once per molecule: ring membership counts are asked for on
        // every candidate test, and recomputing them per atom is quadratic.
        ringCnt_.resize(n);
        for (std::size_t i = 0; i < n; ++i) ringCnt_[i] = ringCountOf(rings, static_cast<int>(i));
        adj_.assign(n, {});
        for (std::size_t i = 0; i < n; ++i) {
            const Atom& a = mol.atoms[i];
            const std::size_t k = std::min(a.nbr.size(), a.bonds.size());
            for (std::size_t j = 0; j < k; ++j) adj_[i].push_back({a.nbr[j], a.bonds[j]});
        }
    }

    void run(const SmartsPattern& pat, std::size_t limit, std::vector<Match>& out) {
        search(pat, -1, -1, limit, &out);
    }

    bool any(const SmartsPattern& pat) { return search(pat, -1, -1, 1, nullptr); }

private:
    // ---- primitive evaluation
    bool evalAtom(const SmartsPattern& p, int node, int ai) {
        const AtomNode& n = p.atomNodes[static_cast<std::size_t>(node)];
        switch (n.op) {
            case ExprOp::Not: return !evalAtom(p, n.lhs, ai);
            case ExprOp::And: return evalAtom(p, n.lhs, ai) && evalAtom(p, n.rhs, ai);
            case ExprOp::Or: return evalAtom(p, n.lhs, ai) || evalAtom(p, n.rhs, ai);
            case ExprOp::Prim: break;
        }
        const Atom& a = mol_.atoms[static_cast<std::size_t>(ai)];
        switch (n.prim) {
            case AtomPrimKind::Any: return true;
            case AtomPrimKind::Element:
                if (a.z != n.value) return false;
                if (n.arom < 0) return true;
                return a.aromatic == (n.arom == 1);
            case AtomPrimKind::AromaticAny: return a.aromatic;
            case AtomPrimKind::AliphaticAny: return !a.aromatic;
            case AtomPrimKind::TotalH: return a.totalH() == n.value;
            case AtomPrimKind::Connectivity: return a.degree() + a.totalH() == n.value;
            case AtomPrimKind::Degree: return a.degree() == n.value;
            case AtomPrimKind::Charge: return a.charge == n.value;
            case AtomPrimKind::RingCount:
                return n.value < 0 ? ringCnt_[static_cast<std::size_t>(ai)] > 0
                                   : ringCnt_[static_cast<std::size_t>(ai)] == n.value;
            case AtomPrimKind::InAnyRing: return ringCnt_[static_cast<std::size_t>(ai)] > 0;
            case AtomPrimKind::RingSize: return inRingOfSize(rings_, ai, n.value);
            case AtomPrimKind::Recursive: return recursive(*p.recursives[static_cast<std::size_t>(n.sub)], ai);
        }
        return false;
    }

    // Memoised, because the same $() is re-tested against the same atom many
    // times: once per candidate filter pass and again on every extension.
    bool recursive(const SmartsPattern& sub, int ai) {
        const std::int64_t key = sub.id * 1000003LL + ai;
        auto it = recCache_.find(key);
        if (it != recCache_.end()) return it->second;
        const bool ok = search(sub, 0, ai, 1, nullptr);
        recCache_.emplace(key, ok);
        return ok;
    }

    bool evalBond(const SmartsPattern& p, int node, const Bond& b) {
        const BondNode& n = p.bondNodes[static_cast<std::size_t>(node)];
        switch (n.op) {
            case ExprOp::Not: return !evalBond(p, n.lhs, b);
            case ExprOp::And: return evalBond(p, n.lhs, b) && evalBond(p, n.rhs, b);
            case ExprOp::Or: return evalBond(p, n.lhs, b) || evalBond(p, n.rhs, b);
            case ExprOp::Prim: break;
        }
        switch (n.prim) {
            case BondPrimKind::Any: return true;
            case BondPrimKind::Single: return !b.aromatic && b.order == 1.0;
            case BondPrimKind::Double: return !b.aromatic && b.order == 2.0;
            case BondPrimKind::Triple: return !b.aromatic && b.order == 3.0;
            case BondPrimKind::Aromatic: return b.aromatic;
            case BondPrimKind::RingBond: return b.inRing;
            case BondPrimKind::Default: return b.aromatic || b.order == 1.0;
        }
        return false;
    }

    // ---- search
    // Returns true as soon as one match exists; when `out` is non-null it
    // collects distinct matches (by matched-atom SET, see header) up to limit.
    bool search(const SmartsPattern& pat, int anchorPat, int anchorMol, std::size_t limit,
                std::vector<Match>* out) {
        const std::size_t np = pat.atoms.size();
        const std::size_t nm = mol_.atoms.size();
        if (np == 0 || np > nm || limit == 0) return false;

        // Prefilter: evaluate each pattern atom's whole expression against every
        // molecule atom once. This is the cheap win - element/degree/aromaticity/
        // charge/H-count all collapse into one candidate list per pattern atom,
        // so the VF2 search never even considers most atoms.
        std::vector<std::vector<int>> cand(np);
        for (std::size_t pi = 0; pi < np; ++pi) {
            for (std::size_t mi = 0; mi < nm; ++mi) {
                if (anchorPat == static_cast<int>(pi) && static_cast<int>(mi) != anchorMol) continue;
                if (evalAtom(pat, pat.atoms[pi].expr, static_cast<int>(mi)))
                    cand[pi].push_back(static_cast<int>(mi));
            }
            if (cand[pi].empty()) return false;
            // A pattern atom needs at least as many incident bonds in the
            // molecule as in the pattern (monomorphism allows extra, never fewer).
            const std::size_t need = pat.atoms[pi].bonds.size();
            auto& c = cand[pi];
            c.erase(std::remove_if(c.begin(), c.end(),
                                   [&](int mi) {
                                       return adj_[static_cast<std::size_t>(mi)].size() < need;
                                   }),
                    c.end());
            if (c.empty()) return false;
        }

        const std::vector<int> order = buildOrder(pat, cand, anchorPat);
        std::vector<int> mapping(np, -1);
        std::vector<char> used(nm, 0);
        std::set<std::vector<int>> seen;
        bool found = false;
        extend(pat, order, cand, 0, mapping, used, seen, limit, out, found);
        return found;
    }

    // Most-constrained-first: start from the pattern atom with the fewest
    // candidates, then always grow along an existing pattern bond so every
    // extension is bond-checked immediately instead of at the end.
    static std::vector<int> buildOrder(const SmartsPattern& pat,
                                       const std::vector<std::vector<int>>& cand, int anchorPat) {
        const int np = static_cast<int>(pat.atoms.size());
        std::vector<char> picked(static_cast<std::size_t>(np), 0);
        std::vector<int> order;
        order.reserve(static_cast<std::size_t>(np));

        auto leastCandidates = [&]() {
            int best = -1;
            for (int i = 0; i < np; ++i) {
                if (picked[static_cast<std::size_t>(i)]) continue;
                if (best < 0 || cand[static_cast<std::size_t>(i)].size() <
                                    cand[static_cast<std::size_t>(best)].size())
                    best = i;
            }
            return best;
        };

        int first = anchorPat >= 0 ? anchorPat : leastCandidates();
        picked[static_cast<std::size_t>(first)] = 1;
        order.push_back(first);

        while (static_cast<int>(order.size()) < np) {
            int best = -1, bestTies = -1;
            for (int i = 0; i < np; ++i) {
                if (picked[static_cast<std::size_t>(i)]) continue;
                int ties = 0;
                for (int bi : pat.atoms[static_cast<std::size_t>(i)].bonds) {
                    const PatternBond& b = pat.bonds[static_cast<std::size_t>(bi)];
                    const int o = (b.a == i) ? b.b : b.a;
                    if (picked[static_cast<std::size_t>(o)]) ++ties;
                }
                if (ties == 0) continue;
                if (best < 0 || ties > bestTies ||
                    (ties == bestTies && cand[static_cast<std::size_t>(i)].size() <
                                             cand[static_cast<std::size_t>(best)].size())) {
                    best = i;
                    bestTies = ties;
                }
            }
            if (best < 0) best = leastCandidates();  // separate '.' component
            picked[static_cast<std::size_t>(best)] = 1;
            order.push_back(best);
        }
        return order;
    }

    void extend(const SmartsPattern& pat, const std::vector<int>& order,
                const std::vector<std::vector<int>>& cand, std::size_t depth,
                std::vector<int>& mapping, std::vector<char>& used,
                std::set<std::vector<int>>& seen, std::size_t limit, std::vector<Match>* out,
                bool& found) {
        if (out != nullptr && out->size() >= limit) return;
        if (found && out == nullptr) return;

        if (depth == order.size()) {
            found = true;
            if (out == nullptr) return;
            std::vector<int> key = mapping;
            std::sort(key.begin(), key.end());
            if (seen.insert(std::move(key)).second) out->push_back(Match{mapping});
            return;
        }

        const int pi = order[depth];
        for (int mi : cand[static_cast<std::size_t>(pi)]) {
            if (used[static_cast<std::size_t>(mi)]) continue;
            if (!bondsOk(pat, pi, mi, mapping)) continue;
            mapping[static_cast<std::size_t>(pi)] = mi;
            used[static_cast<std::size_t>(mi)] = 1;
            extend(pat, order, cand, depth + 1, mapping, used, seen, limit, out, found);
            used[static_cast<std::size_t>(mi)] = 0;
            mapping[static_cast<std::size_t>(pi)] = -1;
            if (found && out == nullptr) return;
            if (out != nullptr && out->size() >= limit) return;
        }
    }

    // Every pattern bond from `pi` to an already-mapped pattern atom must exist
    // in the molecule and satisfy its bond expression. Extra molecule bonds are
    // fine: this is monomorphism.
    bool bondsOk(const SmartsPattern& pat, int pi, int mi, const std::vector<int>& mapping) {
        for (int bi : pat.atoms[static_cast<std::size_t>(pi)].bonds) {
            const PatternBond& pb = pat.bonds[static_cast<std::size_t>(bi)];
            const int po = (pb.a == pi) ? pb.b : pb.a;
            const int mo = mapping[static_cast<std::size_t>(po)];
            if (mo < 0) continue;
            int molBond = -1;
            for (const auto& e : adj_[static_cast<std::size_t>(mi)])
                if (e.first == mo) {
                    molBond = e.second;
                    break;
                }
            if (molBond < 0) return false;
            if (!evalBond(pat, pb.expr, mol_.bonds[static_cast<std::size_t>(molBond)])) return false;
        }
        return true;
    }

    const Molecule& mol_;
    const RingInfo& rings_;
    std::vector<int> ringCnt_;
    std::vector<std::vector<std::pair<int, int>>> adj_;
    std::unordered_map<std::int64_t, bool> recCache_;
};

thread_local std::string tlsError;

}  // namespace

std::optional<SmartsPattern> parseSmarts(std::string_view smarts, std::string* error) {
    if (error != nullptr) error->clear();
    Parser p(smarts);
    SmartsPattern out;
    if (!p.run(out)) {
        if (error != nullptr) *error = p.error();
        return std::nullopt;
    }
    return out;
}

std::optional<SmartsPattern> parseSmarts(std::string_view smarts) {
    return parseSmarts(smarts, &tlsError);
}

const std::string& smartsError() { return tlsError; }

std::vector<Match> matchAll(const SmartsPattern& pat, const Molecule& mol, const RingInfo& rings,
                            std::size_t limit) {
    std::vector<Match> out;
    if (pat.empty() || mol.empty()) return out;
    Matcher m(mol, rings);
    m.run(pat, limit, out);
    return out;
}

std::vector<Match> matchAll(const SmartsPattern& pat, const Molecule& mol, std::size_t limit) {
    const RingInfo rings = perceiveRings(mol);
    return matchAll(pat, mol, rings, limit);
}

bool matches(const SmartsPattern& pat, const Molecule& mol, const RingInfo& rings) {
    if (pat.empty() || mol.empty()) return false;
    Matcher m(mol, rings);
    return m.any(pat);
}

bool matches(const SmartsPattern& pat, const Molecule& mol) {
    const RingInfo rings = perceiveRings(mol);
    return matches(pat, mol, rings);
}

PreparedMolecule prepareMolecule(const Molecule& mol) {
    PreparedMolecule out;
    out.mol = mol;
    out.rings = perceiveRings(out.mol);
    annotateRings(out.mol, out.rings);
    perceiveAromaticity(out.mol, out.rings);
    return out;
}

}  // namespace biocad::chem
