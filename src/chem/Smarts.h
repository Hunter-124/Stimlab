// chem/Smarts.h - SMARTS substructure query language: parser + VF2 matcher.
//
// WHY this exists: every published biotransformation rule set (GLORYx, SyGMa),
// every structural-alert set (PAINS, Brenk, Glaxo) and every site-of-metabolism
// model (SMARTCyp) is expressed in SMARTS. Without a matcher those rule packs
// cannot be data at all, and the alternative - hand-coding each rule in C++ -
// is exactly the fabricated, uncitable logic this project is removing.
//
// SCOPE AND HONESTY: a primitive that is parsed but not enforced would silently
// match the wrong atoms, which is worse than refusing the pattern. So every
// primitive listed below is enforced, and everything else (isotopes, ^n
// hybridisation) is a hard parse error naming itself. The one deliberate
// exception is tetrahedral chirality (@ / @@), which is accepted and ignored
// because chem::Molecule stores no stereochemistry at all: there is nothing to
// match against, and the resulting match set is a documented superset.
//
// Supported atom primitives inside [ ]:
//   *  any                      a / A   aromatic / aliphatic
//   C  aliphatic element        c       aromatic element        #n  by number
//   Hn total-H count (H == H1)  Xn      total connectivity (heavy degree + H)
//   Dn explicit (heavy) degree  +/-/+n/-n/++/-- formal charge
//   R  in >=1 SSSR ring         Rn      in exactly n SSSR rings (R0 == none)
//   r  in a ring                rn      in a ring of size n
//   $(...) recursive SMARTS anchored on this atom
//   :n atom map (retained for a later SMIRKS phase, never matched on)
// Logical operators with SMARTS precedence !  >  &  >  ,  >  ;
//
// Supported bond primitives: -  =  #  :  ~  @  and the same logical operators.
// / and \ are read as single bonds (no stereo in the graph). An omitted bond is
// "single or aromatic", as SMARTS requires. '.' separates components.
//
// MATCH SEMANTICS: subgraph *monomorphism*, not isomorphism - the molecule may
// carry extra bonds between mapped atoms. That is what SMARTS means, and
// isomorphism would make ring and fused-ring patterns fail.
//
// UNIQUENESS CONVENTION (fixed, and consistent everywhere): two matches are
// distinct iff their *sets* of matched molecule atoms differ. Order is ignored,
// so the automorphisms of a symmetric pattern do not multiply the count:
// c1ccccc1 matches benzene once and naphthalene twice (its two six-cycles).
// Match::atoms holds, in pattern-atom order, the first mapping found for that
// set, so a caller can still tell which pattern atom hit which molecule atom.
#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "chem/Molecule.h"
#include "chem/Rings.h"

namespace biocad::chem {

// ---- atom expression AST -------------------------------------------------

enum class AtomPrimKind {
    Any,           // *
    Element,       // C, c, #6      value = Z, arom = -1 any / 0 aliphatic / 1 aromatic
    AromaticAny,   // a
    AliphaticAny,  // A
    TotalH,        // Hn
    Connectivity,  // Xn
    Degree,        // Dn
    Charge,        // +n / -n
    RingCount,     // R (value < 0 means ">= 1"), Rn
    InAnyRing,     // r
    RingSize,      // rn
    Recursive,     // $(...)
};

enum class ExprOp { Prim, Not, And, Or };

struct AtomNode {
    ExprOp       op = ExprOp::Prim;
    AtomPrimKind prim = AtomPrimKind::Any;
    int          value = 0;
    int          arom = -1;  // Element only
    int          sub = -1;   // Recursive only: index into SmartsPattern::recursives
    int          lhs = -1;   // Not/And/Or operands: indices into SmartsPattern::atomNodes
    int          rhs = -1;
};

// ---- bond expression AST -------------------------------------------------

enum class BondPrimKind {
    Any,       // ~
    Single,    // -
    Double,    // =
    Triple,    // #
    Aromatic,  // :
    RingBond,  // @
    Default,   // omitted: single or aromatic
};

struct BondNode {
    ExprOp       op = ExprOp::Prim;
    BondPrimKind prim = BondPrimKind::Default;
    int          lhs = -1;
    int          rhs = -1;
};

// ---- pattern -------------------------------------------------------------

struct PatternAtom {
    int expr = -1;                    // root index into SmartsPattern::atomNodes
    int map = 0;                      // atom map label from :n, 0 when absent
    std::vector<int> bonds;           // incident indices into SmartsPattern::bonds
};

struct PatternBond {
    int a = 0;
    int b = 0;
    int expr = -1;                    // root index into SmartsPattern::bondNodes
};

struct SmartsPattern {
    std::string                                 source;
    std::vector<AtomNode>                       atomNodes;
    std::vector<BondNode>                       bondNodes;
    std::vector<PatternAtom>                    atoms;
    std::vector<PatternBond>                    bonds;
    // Recursive subpatterns. Each carries a process-unique id used as the
    // memoisation key, because the same $() is re-tested against the same atom
    // many times during one search.
    std::vector<std::shared_ptr<SmartsPattern>> recursives;
    std::int64_t                                id = 0;

    [[nodiscard]] bool empty() const { return atoms.empty(); }
};

// Parse SMARTS. Returns nullopt on malformed input and, when `error` is given,
// fills it with a message naming the problem and its offset so a bad rule pack
// can name its bad rule. A malformed pattern is never a pattern that matches
// nothing.
std::optional<SmartsPattern> parseSmarts(std::string_view smarts, std::string* error);
std::optional<SmartsPattern> parseSmarts(std::string_view smarts);

// Last parse error from the single-argument overload (thread-local).
const std::string& smartsError();

struct Match {
    std::vector<int> atoms;  // molecule atom indices, in pattern-atom order
};

// Ring-dependent primitives (R, Rn, r, rn, @) need ring perception, and
// aromaticity must already be perceived on `mol` (see chem/Aromaticity.h) -
// these functions never mutate the molecule they are asked about. The RingInfo
// overload exists so a rule-pack sweep perceives rings once, not per pattern.
std::vector<Match> matchAll(const SmartsPattern& pat, const Molecule& mol, const RingInfo& rings,
                            std::size_t limit = 1000);
std::vector<Match> matchAll(const SmartsPattern& pat, const Molecule& mol,
                            std::size_t limit = 1000);
bool matches(const SmartsPattern& pat, const Molecule& mol, const RingInfo& rings);
bool matches(const SmartsPattern& pat, const Molecule& mol);

// Convenience for callers holding a raw parsed SMILES: copy, perceive rings and
// aromaticity from the graph, and hand back both so uppercase-SMILES input is
// matched correctly by aromatic patterns.
struct PreparedMolecule {
    Molecule mol;
    RingInfo rings;
};
PreparedMolecule prepareMolecule(const Molecule& mol);

}  // namespace biocad::chem
