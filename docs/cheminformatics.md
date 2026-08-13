# The in-house cheminformatics engine

This page is written for a reader deciding **how much of `src/chem` to trust**. It describes
what the engine computes, how, and - at equal length - where it is wrong. Every claim below
names the file it came from, and every number was produced by running the code in this tree on
this machine (see [Measured evidence](#measured-evidence) for the exact figures and how they
were obtained).

## 1. Why there is no RDKit

`vcpkg.json:29` records the reason verbatim:

> NOTE: there is intentionally no rdkit - BioCAD has its own in-house cheminformatics engine in
> src/chem (rdkit is not a port in this vcpkg registry).

That is a supply-chain fact, not a preference: the dependency simply cannot be resolved for the
Windows triplets this app builds for, so there was no version of this project that could link
RDKit and be built by its own presets.

**State the consequence plainly.** Everything in `src/chem` - the SMILES parser
(`src/chem/Smiles.cpp`), ring perception (`src/chem/Rings.cpp`), aromaticity
(`src/chem/Aromaticity.cpp`), the canonical writer (`src/chem/Canonical.cpp`), the SMARTS
engine (`src/chem/Smarts.cpp`), the descriptors (`src/chem/Descriptors.cpp`), the
Wildman-Crippen logP / molar-refractivity implementation (`src/chem/Crippen.cpp`), the
bioactivation alert screen (`src/chem/Alerts.cpp`), the functional-group perception
(`src/chem/Analysis.cpp`) and the 3D embedder (`src/chem/Embed3D.cpp`) - is a
**reimplementation of a published algorithm**. A reimplementation is a place a bug can hide that
a widely-used library would not have, because a library with tens of thousands of users has had
its edge cases found by those users. This engine has had its edge cases found by a Catch2 suite
(`tests/test_chem.cpp`, `test_chem_rings.cpp`, `test_chem_canonical.cpp`, `test_chem_smarts.cpp`,
`test_chem_alerts.cpp`, `test_embed3d.cpp`) and by the probes quoted on this page. Those are not
the same thing, and this document exists so nobody has to guess which one they are relying on.

The practical rule that follows: **BioCAD's descriptors are BioCAD's descriptors.** Use them to
rank compounds inside one BioCAD run. Do not quote them as AlogP or as Ertl TPSA values and do
not diff them digit-for-digit against a toolkit that computes the same-named quantity.

## 2. Ring perception - SSSR over a GF(2) cycle basis

`src/chem/Rings.h` / `src/chem/Rings.cpp`. Before this existed the only ring information in the
tree came from the parser's bridge detection (`src/chem/Smiles.cpp:37-60`, a DFS low-link scan),
which answers exactly one question - "is this bond in some cycle?" - and says nothing about ring
size or membership count.

The algorithm, as documented in `Rings.h` and implemented in `Rings.cpp`:

1. **The answer's size is known in advance.** Euler's formula fixes the dimension of the cycle
   space: `dim = |bonds| - |atoms| + |components|`. An acyclic (or disconnected acyclic) input
   has `dim = 0`. The same expression is what `ringCount()` in `src/chem/Descriptors.cpp:96-118`
   returns, computed independently by a flood fill over components.
2. **Candidates come from Horton's set**: for every atom `r` and every bond `(x, y)`, the cycle
   `SP(r,x) + (x,y) + SP(y,r)` when the two shortest paths are vertex-disjoint. That set is
   proven to contain a minimum-weight cycle basis.
3. **Acceptance is a GF(2) independence test.** Each candidate is a bond-incidence bit vector
   (`BitRow`, a `std::vector<std::uint64_t>` bitset - `Rings.cpp:12-40`); candidates are sorted
   smallest-first and accepted greedily only when linearly independent of the already-accepted
   rings under Gaussian elimination over GF(2).
4. **Termination is structural**: the loop stops the moment the accepted count reaches `dim`, so
   the routine can neither over-report nor spin.

`RingInfo` (`Rings.h`) carries `atomRings[i]` and `bondRings[i]` for the same ring, in matching
walk order, smallest ring first. `annotateRings()` writes `Atom::inRing` / `Bond::inRing` back
onto the graph. The three query helpers are what the SMARTS ring primitives consume:
`ringSizeOf()` (smallest ring containing an atom, `rN`), `inRingOfSize()` (any ring of exactly
that size - a fusion atom of indole answers true for both 5 and 6) and `ringCountOf()` (`Rn`).

### What SSSR does not give you

**An SSSR is not unique**, and the header says so before this document does. Cubane is the
standard counterexample: six four-membered faces, but a cycle space of dimension
`12 - 8 + 1 = 5`, so any SSSR must discard one face and *which* face is discarded is arbitrary.
`perceiveRings()` returns **one valid SSSR**, chosen deterministically (ordered by ring size,
then by canonicalised atom sequence, so the same molecule always yields the same rings). That is
what ring-size queries and Hueckel counting need. Code that needs *every* small ring wants the
relevant-cycles set, which this deliberately does not compute - so do not read a cage-system
ring assignment as "the" ring system.

## 3. Aromaticity - perceived from the graph, not from the spelling

`src/chem/Aromaticity.h` / `.cpp`. The bug this fixes was real and load-bearing: before
perception existed, lowercase `c1ccccc1` was aromatic and the exactly equivalent Kekule form
`C1=CC=CC=C1` was not, so every aromatic SMARTS primitive, every structural alert and every
biotransformation rule would fire on one spelling of benzene and silently miss the other.

The model is per-ring Hueckel counting over the SSSR: every ring atom must be sp2-capable (one
sp3 centre disqualifies the ring outright), the pi-electron contributions are summed, and the
ring is aromatic when the sum is `4n + 2`. Rings are then re-evaluated **to a fixed point** so
fused systems resolve - naphthalene's second ring counts only 4 electrons in isolation from
Kekule input, and is upgraded to 6 once the first ring is known aromatic.

### Pi-electron contribution table

Reproduced from `src/chem/Aromaticity.h`, which is the authoritative copy:

| Atom in ring | e- | Why |
|---|---|---|
| C with a ring double bond (order 2 or 1.5) | 1 | one electron of the shared C=C pi bond lies in this ring |
| C with an exocyclic double bond to O/N/S (carbonyl / thiocarbonyl / imine carbon) | 0 | the pi pair is pulled out of the ring by the electronegative atom; the carbon stays sp2 |
| C with an exocyclic double bond to a carbon that is itself perceived aromatic (fusion carbon) | 1 | cross-ring conjugation at a fused bond; this is the fixed-point rule |
| C with an exocyclic double bond to a non-aromatic carbon | 0 | fulvene-type exocyclic alkene donates nothing |
| C(-) carbanion, no ring double bond | 2 | filled p orbital (cyclopentadienyl anion) |
| C(+) carbocation, no ring double bond | 0 | empty p orbital (tropylium) |
| N/P with a ring double bond and no ring H | 1 | pyridine-type: the lone pair is in the sigma plane, not the pi system |
| N/P with an H, or three sigma connections, and only single/aromatic ring bonds | 2 | pyrrole-type: the lone pair completes the sextet |
| N(-) amide / azolate anion | 2 | filled p orbital |
| O/S/Se with two ring single (or aromatic) bonds | 2 | furan/thiophene: one lone pair enters the pi system |
| O(+)/S(+) with a ring double bond | 1 | pyrylium / thiopyrylium |
| B (three connections) | 0 | empty p orbital (borole, borazine-type rings) |
| anything else | - | not sp2: the ring is rejected |

### What this model gets wrong

Also from `Aromaticity.h`, because a heuristic ring model has failure modes and hiding them is
the only unacceptable option:

- **2-pyridone, cytosine, and caffeine's pyrimidinedione ring are reported AROMATIC.** The
  exocyclic C=O contributes 0 and the amide N contributes 2, summing to 6. This matches RDKit's
  default model and ChEMBL-style conventions, but a chemist may reasonably call these amides.
  Measured: `O=c1cccc[nH]1` reports 6 aromatic atoms.
- **Aromaticity delocalised over a perimeter rather than over an SSSR ring is missed outright.**
  Azulene is the verified example: its SSSR is a 5-ring and a 7-ring counting 5 and 7 electrons,
  while the real 10-electron system runs around the bicyclic perimeter, which is not an SSSR
  ring. Measured: `c1ccc2cccc2cc1` reports **0** aromatic atoms. Macrocyclic annulenes such as
  [18]annulene fail the same way. Fixing this needs perimeter / electron-flow perception, which
  this deliberately is not.
- **No anti-aromaticity reasoning.** A 4n ring is simply "not aromatic", never flagged as
  destabilised (cyclobutadiene).
- **No charge-separated / mesomeric-betaine rings, no N-oxides, no metal-coordinated rings.**
- **The SSSR itself is not unique for cage systems** (section 2), so a cage ring assignment is
  one valid choice rather than the only one.

`perceiveAromaticity()` sets **flags only** (`Atom::aromatic`, `Bond::aromatic`) and never
produces a `data::Quantity`: a structural perception is not a measured or predicted number, so
no provenance tier is involved. A bond is flagged aromatic only if it lies inside a ring that
was perceived aromatic, so biphenyl's inter-ring single bond stays plain.
`normalizeAromaticBondOrders()` is a deliberately separate opt-in step that sets aromatic bond
orders to 1.5; doing it unconditionally would destroy the Kekule structure that valence, implicit
H counts and carbonyl detection rely on.

### Perception is not automatic - a finding, now closed at the source

`perceiveAromaticity()` / `perceiveRingsAndAromaticity()` must be **called**, and nothing calls
it implicitly, so for a while the descriptor path silently depended on how a pack author had typed
the SMILES. Every function whose answer depends on aromaticity now perceives on its **own private
copy**, which means the caller cannot forget and the caller's molecule is never mutated:
`chem::crippen()` (`src/chem/Crippen.cpp:208`), `chem::detectGroups()`
(`src/chem/Analysis.cpp:166`), `chem::screenAlerts()` (`src/chem/Alerts.cpp:141`) and
`chem::tpsa()` / `chem::tpsaIncludingSulfurAndPhosphorus()` (`src/chem/Descriptors.cpp:209-215`,
whose bodies moved into a static `tpsaPerceived()` so the public entry points can call
`prepareMolecule()` first).

That fixes all ~15 descriptor call sites at once, including `src/modules/RealBackend.cpp:41-53`,
which still hands over the raw output of `chem::parseSmiles`. It matters because
`assets/packs/cns-monoamine.json:10` stores caffeine in Kekule form
(`CN1C=NC2=C1C(=O)N(C(=O)N2C)C`) while the same compound is elsewhere written lowercase.
Measured on the tree as it stands - identical answers from both spellings:

| Compound | lowercase TPSA | Kekule TPSA | difference |
|---|---|---|---|
| caffeine | 60.260 | 60.260 | 0.0e+00 |
| aspirin | 63.600 | 63.600 | 7.1e-15 |
| furan | 13.140 | 13.140 | 0.0e+00 |
| pyridine | 12.890 | 12.890 | 0.0e+00 |

`crippen().logP` was already spelling-independent for the same reason and stays at -1.029 for
both caffeine spellings. **Historical note, because the earlier text of this page said otherwise:**
before that change the Kekule spelling gave TPSA 56.22 against the perceived 60.26, a 4.04 A^2
swing caused purely by how the SMILES was typed. The fix went into the descriptor functions rather
than into their callers, which is the right place: a rule every caller must remember is a rule
some caller will forget.

## 4. Canonical SMILES and `graphHash`

`src/chem/Canonical.h` / `.cpp`. Without a canonical form, the same molecule entered two ways is
two unrelated strings: a library cannot be deduplicated, a cache cannot be keyed, and "are these
the same structure?" has no answer. The only other writer in the tree, `sketchToSmiles` in
`src/ui/Panels.cpp`, walks the graph in input order and is explicitly not canonical.

`canonicalRanks()` is Morgan-style refinement: initial invariants (element, degree, charge, total
H, aromaticity, ring membership) are refined by iterated neighbour-rank hashing until the
partition stops changing, then remaining ties are broken one atom at a time (lowest class first,
then lowest atom index) with a refinement after each break, so the result is always a total
order. `canonicalSmiles()` does a DFS from the lowest-ranked atom with neighbours visited in
canonical-rank order, allocates ring-closure digits in order of first need, and joins
disconnected components with `.` ordered by their own canonical strings.

**What is claimed.** `canonicalSmiles()` is order-independent: it roots its walk at the canonical
rank rather than at atom 0, so atoms the refinement cannot separate are interchangeable and
produce the same string. The contract is round-tripping -
`parseSmiles(canonicalSmiles(m))` yields a graph with the same canonical string.
`tests/test_chem_canonical.cpp` verifies this experimentally over re-rooted inputs; the probe in
section 7 reproduces it for three hand-written spellings of aspirin.

**What is not claimed.** Tie-breaking by atom index makes the ranking canonical *for a given
atom ordering*. **Full graph-automorphism canonicalisation** - a ranking provably invariant
under every relabelling of the input, which a proper refinement/backtracking canonicaliser
(nauty, RDKit's canonical ranker) provides - is a strictly stronger property that this engine
does **not** claim. The order-independence that callers depend on is demonstrated
experimentally, not proved.

Two further caveats that follow from the code:

- The writer **does not perceive** aromaticity; it writes what the graph is annotated with. A
  Kekule graph is written Kekule. Two spellings of the same molecule therefore converge only if
  both were perceived (or normalised) first.
- `graphHash()` is a **64-bit hash, not an identity proof**. `Canonical.h` states the
  consequence directly: a collision is a wrong cache hit - one molecule silently answering for
  another. The hash folds the atom and bond counts in alongside the canonical string so the
  cheapest structural mismatches cannot collide at all, but any caller that cannot tolerate a
  wrong answer must compare `canonicalSmiles()` as well.

## 5. SMARTS - parser plus VF2 matcher

`src/chem/Smarts.h` / `.cpp` (801 lines). This exists because every published rule set BioCAD
wants to consume as *data* - GLORYx and SyGMa biotransformations, PAINS / Brenk / Glaxo
structural alerts, SMARTCyp site-of-metabolism rules - is written in SMARTS. The alternative,
hand-coding each rule in C++, is exactly the uncitable logic this project is removing.

### Supported primitives

Atom primitives inside `[ ]`, all of them enforced (`Smarts.h`, `evalAtom()` in
`Smarts.cpp:533-563`):

| Primitive | Meaning |
|---|---|
| `*` | any atom |
| `a` / `A` | aromatic / aliphatic |
| `C` / `c` / `#n` | element, aliphatic / aromatic / by atomic number |
| `Hn` | total H count (`H` == `H1`) |
| `Xn` | total connectivity (heavy degree + H) |
| `Dn` | explicit heavy degree |
| `+`, `-`, `+n`, `-n`, `++`, `--` | formal charge |
| `R` / `Rn` / `R0` | in >= 1 SSSR ring / in exactly n rings / in none |
| `r` / `rn` | in a ring / in a ring of size n (`rn` rejected below 3, `Smarts.cpp:471`) |
| `$(...)` | recursive SMARTS anchored on this atom |
| `:n` | atom map, retained for a later SMIRKS phase, never matched on |

Bond primitives: `-` `=` `#` `:` `~` `@`, plus the same logical operators. `/` and `\` are read
as plain single bonds because the graph stores no stereochemistry (`Smarts.cpp:214-217`). An
omitted bond means "single or aromatic", as SMARTS requires (`BondPrimKind::Default`). `.`
separates components, and the match order builder treats a disconnected component by falling
back to the least-constrained unpicked atom (`buildOrder`, `Smarts.cpp:670`).

Logical precedence is the SMARTS ladder **`!` > `&` > `,` > `;`**, implemented as one
recursive-descent layer per level for atoms and mirrored for bonds
(`parseBondLow` / `parseBondOr` / `parseBondAnd` / `parseBondUnary`, `Smarts.cpp:160-220`).
Adjacency with no operator is the implicit high-precedence `&`.

### Unsupported primitives are hard errors

A primitive that parsed but was not enforced would silently match the wrong atoms, which is
worse than refusing the pattern, so anything unsupported fails the parse with a message naming
itself and its offset (`fail()`, `Smarts.cpp:46`). Measured:

```
[13C] -> isotope specifications are not supported (the graph stores no isotopes) at offset 1 in "[13C]"
[C^3] -> hybridisation primitive ^n is not supported at offset 2 in "[C^3]"
```

`parseSmarts()` returns `std::nullopt` and fills the caller's `error` string, so a bad rule pack
can name its bad rule; **a malformed pattern is never a pattern that matches nothing.**

The one deliberate exception is tetrahedral chirality: `@` / `@@` is **accepted and ignored**
(`Smarts.cpp:436`), because `chem::Molecule` stores no stereochemistry at all - there is nothing
to match against. The documented consequence is that such a pattern's match set is a **superset**
of the correct one. Measured: `[C@H](N)(C)C(=O)O` parses with no error.

### Match semantics

- **Subgraph monomorphism, not isomorphism.** The molecule may carry extra bonds between mapped
  atoms; only the bonds the *pattern* names are checked (`bondsOk`, and the comment at
  `Smarts.cpp:718-720`). That is what SMARTS means, and isomorphism would make ring and
  fused-ring patterns fail.
- **Uniqueness convention, as actually implemented:** two matches are distinct iff their *sets*
  of matched molecule atoms differ. `extend()` sorts the mapping into a key and inserts it into a
  `std::set<std::vector<int>>` before emitting (`Smarts.cpp:690-706`), so pattern automorphisms
  do not multiply the count. Measured: `c1ccccc1` matches benzene **once** and naphthalene
  **twice** (its two six-cycles). `Match::atoms` still holds the *first* mapping found for that
  set, in pattern-atom order, so a caller can tell which pattern atom hit which molecule atom -
  but only one of the symmetric mappings is reported.
- **`matchAll` is bounded** by a `limit` parameter defaulting to 1000; a pathological pattern
  truncates rather than running away.

### Search strategy and cost

Per molecule the matcher caches ring-membership counts and a combined (neighbour, bond) adjacency
list (`Matcher` constructor, `Smarts.cpp:511-523`). Per search it prefilters: each pattern atom's
whole expression is evaluated against every molecule atom once, giving one candidate list per
pattern atom, then candidates with fewer incident bonds than the pattern atom needs are dropped
(monomorphism allows extra bonds, never fewer). The VF2 extension order is most-constrained-first
and always grows along an existing pattern bond, so every extension is bond-checked immediately
(`buildOrder`, `Smarts.cpp:637-682`).

Recursive `$()` is **memoised**: each subpattern carries a process-unique id, and
`recursive()` keys a cache on `(subpattern id, atom index)` because the same `$()` is re-tested
against the same atom once per candidate-filter pass and again on every extension
(`Smarts.cpp:565-573`).

**Measured cost** on this machine (AMD Ryzen 5 7600X3D, `g++ -O2`, 200,000 iterations, full
`matchAll` including the prefilter and the result vector), pattern `[CX3](=O)[OX2H1]` against
ibuprofen: **0.73 us per (pattern, molecule) call**. The recursive spelling of the same query,
`[$([CX3](=O)[OX2H1])]`, costs 1.48 us. Both agree with the roughly-1-us order of magnitude the
engine was designed around; the explosion risk in rule-pack work is a metabolite tree, not the
matcher.

Ring-dependent primitives (`R`, `Rn`, `r`, `rn`, `@`) need ring perception and aromaticity must
already be perceived: **the matcher never mutates the molecule it is asked about.**
`prepareMolecule()` (`Smarts.h`) is the convenience that copies, perceives rings and perceives
aromaticity, returning a `PreparedMolecule { Molecule mol; RingInfo rings; }` so uppercase SMILES
input is matched correctly by aromatic patterns.

### What the matcher drives today

With a matcher in the tree, three things that used to be hand-coded C++ conditions are now data,
which is the whole point of building it:

| Consumer | Pack | Count | Failure behaviour |
|---|---|---|---|
| `chem::detectGroups()` (`src/chem/Analysis.cpp:159`) - the `FunctionalGroups` flags the pharmacology, ADMET and stability heuristics key on | `assets/packs/rules/functional-groups.json` | 24 group rules | `groupPackErrors()` names a missing pack, an unknown group key or an unparsable SMARTS. A flag left false by a broken rule is indistinguishable from an honestly false flag, so the breakage is surfaced rather than swallowed |
| `chem::screenAlerts()` (`src/chem/Alerts.h`) - bioactivation liability flags | `assets/packs/rules/alerts-bioactivation.json` | 15 alerts, each with a `mechanism` stated as a metabolic route and its own citation | `alertPackErrors()` lists any rule whose SMARTS failed to parse, with the parser's message; a rule that will not parse is dropped and reported, never treated as a rule that matches nothing. `alertRuleCount() == 0` with errors means the screen is inoperative and the UI must say so |
| `chem::crippen()` atom typing (section 6) | `assets/packs/descriptors/crippen.json` | 106 pattern entries over 67 heavy-atom classes | `CrippenResult::ok == false` with a `note` naming the file |

Measured: the packs load clean (`alert rules loaded=15 errors=0 groupPackErrors=0`), and
paracetamol raises exactly one alert - `para-aminophenol`, "para-Aminophenol / anilide
(quinone-imine former)", `warn=1`, 8 matched atoms, cited to Stepan et al. 2011 - while
`detectGroups` reports `phenol`, `anilide`, `amide` and `aromaticRing`. An empty alert result is
**not** a safety claim: the pack is a short in-house list, and absence of a listed motif says
nothing about the motifs it does not list. `Alerts.h` also records why PAINS, Brenk and the
ChEMBL-derived alert sets are deliberately *not* bundled: their BSD-3 and CC BY-SA 3.0 obligations
need their own pack and their own `NOTICE` entry.

The metabolite side is deliberately **not** a rule set: `assets/packs/rules/metabolism-facts.json`
carries 28 curated biotransformations with 19 distinct citations and enumerates nothing, because
rule-based metabolite enumeration has published precision in the single-to-low-double digits.

## 6. Descriptors - which are exact, which are close, and which are neither

`src/chem/Descriptors.h` / `.cpp`.

**Faithful to their definition** (they are counting rules, and the counting rule is the
definition):

| Descriptor | Implementation | Note |
|---|---|---|
| `molecularFormula` | Hill order, C then H then alphabetical (`:52-77`) | includes implicit H |
| `molecularWeight` | sum of standard atomic weights incl. H (`:44-50`) | average mass, from `src/chem/Periodic.cpp` |
| `hbdCount` / `hbaCount` | N/O bearing H; N + O count (`:79-92`) | the Lipinski counting convention, which is itself crude - `hbaCount` counts every N and O, amides and nitro groups included |
| `rotatableBondCount` | acyclic single bonds between two atoms of degree >= 2, amides excluded (`:94-105`) | the Veber definition |
| `ringCount` | `E - V + C` cyclomatic number (`:107-129`) | agrees with the SSSR size by construction |
| `formalCharge`, `heavyAtomCount`, `aromaticAtomCount`, `fractionCsp3` | direct counts (`:37-42`, `:131-152`) | `fractionCsp3` and `aromaticAtomCount` depend on aromaticity having been perceived |

**TPSA is close, not exact, and the convention matters more than the residual.**
`tpsa()` (`src/chem/Descriptors.cpp`) implements the Ertl fragment scheme over
nitrogen and oxygen. Cross-checked against RDKit 2026.03.5 across all 69 shipped
library compounds:

| Convention | Mismatches (>0.1 A^2) | Worst |
|---|---:|---:|
| N,O only (the default) | 10 of 69 | 1.56 A^2 (caffeine) |
| Including S and P | 11 of 69 | 13.5 A^2 (N-acetylcysteine) |

An earlier revision summed sulfur and phosphorus too, which is defensible - Ertl
tabulates them - but it put famotidine at 237.75 A^2 against a reference 175.83 A^2.
That 62 A^2 gap straddles Veber's 140 A^2 oral-bioavailability cutoff, and every
published TPSA threshold was derived on the N,O-only sum. So `tpsa()` is N,O-only
and `tpsaIncludingSulfurAndPhosphorus()` is a separate, separately-labelled number.
Two tools reporting "TPSA" can differ by this much purely on convention, which is
why any readout has to state which one it is.

The 10 residual mismatches are all nitrogen typing and all under 2 A^2: exact
multiples of 0.52 A^2 (celecoxib, indomethacin -0.52; phenazone, ondansetron,
theobromine, theophylline -1.04; caffeine -1.56) or +0.22 A^2 (NMN, thiamine,
berberine). They come from the aromaticity model over-assigning aromaticity to fused
pyrimidinedione and pyrazolone rings - the exact limitation `Aromaticity.h` documents
about itself - which then picks an aromatic instead of an aliphatic nitrogen
contribution. A known model boundary showing up downstream, not an arithmetic error.

**logP is now faithful to the method - which is not the same as being right.** `crippenLogP()`
(`src/chem/Descriptors.h:21-24`) delegates to `chem::crippen()` (`src/chem/Crippen.h`,
`src/chem/Crippen.cpp`), which is a full implementation of Wildman & Crippen's atomic-contribution
scheme (the header cites J Chem Inf Comput Sci 1999;39(5):868-873, doi:10.1021/ci990307l). Three
properties of that implementation matter to a reader deciding whether to trust it:

- **The parameters are data, not code.** The pack file
  `assets/packs/descriptors/crippen.json` holds **106 ordered pattern entries covering the
  method's 67 distinct heavy-atom classes** - 22 of those classes need more than one SMARTS to
  express (C1 takes three, C6 and C26 four each) - with each entry carrying its SMARTS, its logP
  contribution and its molar-refractivity contribution. Counted from the file, not asserted.
  Matching runs through the SMARTS engine of section 5.
  Classification is **first match wins in the array order the file ships**, which is how the
  paper's decision table works, so re-ordering the array changes the answer. The pack's own
  `order` field says so, and flags two deliberate order flips (O12 before O7, S2 before S1) that
  come from the reference file.
- **Hydrogens are derived, not matched**, because `chem::Molecule` carries no explicit H atoms:
  `Crippen.h` documents the exact H1-H4 mapping from the heavy atom's environment, itself first
  match wins (a phenol OH is H2, a carboxylic acid OH is H4).
- **Failure is visible.** `CrippenResult::ok` is false with a `note` naming the pack when the
  descriptor file is missing or malformed, `atomTypes` records the class chosen for every heavy
  atom for audit, and `unclassified` lists any atom no class matched. A missing pack yields an
  explicit failure, not a silent zero. `crippenCitation()` returns the exact string every
  `Quantity` built from this number must carry:
  `"Wildman-Crippen (J Chem Inf Comput Sci 1999;39:868-873), method RMS ~0.67 log units"`.

**What the numbers now look like.** Measured over the ten-compound reference set that
`tests/test_chem.cpp:107-132` uses (its `ref` column is experimental / consensus logP, not
Crippen-computed values):

| Compound | `crippen().logP` | Experimental reference | Deviation |
|---|---|---|---|
| amphetamine | 1.576 | 1.76 | -0.184 |
| methamphetamine | 1.837 | 2.07 | -0.233 |
| MDMA | 1.566 | 2.15 | -0.584 |
| MDA | 1.305 | 1.64 | -0.335 |
| phenethylamine | 1.188 | 1.41 | -0.222 |
| tyramine | 0.893 | 0.86 | +0.033 |
| ephedrine | 1.328 | 1.13 | +0.198 |
| cathinone | 1.217 | 0.59 | +0.627 |
| acetaminophen | 1.351 | 0.46 | +0.891 |
| 4-fluoroamphetamine | 1.715 | 1.90 | -0.185 |

**Measured MAE 0.349, RMS 0.430 log units** against experiment over that set - inside the
method's own published RMS of about 0.67, which is the most that faithfulness to Wildman-Crippen
can buy. The suite enforces both bounds in the method's own terms
(`REQUIRE_THAT(got, WithinAbs(c.ref, 1.35))` per compound - twice the published RMS - and
`REQUIRE(rms < 0.67)` in aggregate, `tests/test_chem.cpp:107-132`), with the reason stated in the
source: the previous +/-0.7 band only held because the ad-hoc estimator it replaced had been
tuned to look good on these ten rows.

**The stronger claim, checked against the reference implementation.** RDKit 2026.03.5 is installed
on this development machine purely as an *oracle* - it is not a dependency, is not linked, and is
not in `vcpkg.json` - and all 69 compounds of the shipped packs were compared value for value.
Reproduced independently for this page (harness five below): **0 mismatches on logP** (tolerance
0.005), **0 on molar refractivity** (0.02) and **0 on molecular formula** (exact string). Exact-value
oracle cases from that check are committed in `tests/test_chem_crippen.cpp`. That is a stronger
claim than agreeing with a handful of literature values: the implementation *is* the published
method, and whether the published method is right about a given compound is the separate question
the table above answers.

The same run produced the TPSA figures earlier in this section, also reproduced here: 10 of 69
compounds deviate by more than 0.1 A^2 on the N,O-only default (caffeine -1.56; phenazone,
ondansetron, theobromine, theophylline -1.04; celecoxib, indomethacin -0.52; NMN, thiamine,
berberine +0.22) and 11 of 69 with sulfur and phosphorus included (worst N-acetylcysteine -13.5,
then NMN -9.59). Every N,O-only residual is nitrogen typing that traces to the aromaticity model's
fused-ring boundary (section 3), and every one is under 2 A^2.

One consequence of getting the atom typing right was a real parser fix in
`chem::assignImplicitHydrogens` (`src/chem/Smiles.cpp`): counting an aromatic bond as order 1.5
gave thiophene's `s` a phantom implicit hydrogen, so `c1ccsc1` and `C1=CSC=C1` were different
molecules. An aromatic bond now contributes one sigma plus one pi for aromatic C/N, and aromatic
atoms no longer promote to a higher valence. Verified after the fix: `c1ccsc1` and `C1=CSC=C1`
both give `C4H4S` at 84.136, `c1ccccc1` and `C1=CC=CC=C1` both give `C6H6` at 78.114, pyrrole
`c1cc[nH]c1` gives `C4H5N` and pyridine `c1ccncc1` gives `C5H5N`.

**A finding about the older accuracy claim.** The figures this page was originally briefed with -
"caffeine -1.269 vs -0.07, ibuprofen 2.618 vs 3.97, dopamine +0.298 vs -0.98, divergence up to
1.35 log units" - compare a Crippen calculation against **experimental** logP, which is not what
the method predicts. Measured on the landed implementation: caffeine -1.029, ibuprofen 3.073,
dopamine 0.599, aspirin 1.310, benzene 1.687. Benzene and aspirin land on their reference values
to within 0.003, while caffeine (0.96 off experiment) and dopamine (1.58 off) are exactly the two
families the test file excludes as unfair oracles for an additive model, with its reason stated
in the source: the additive scheme overweights fused polar heteroaromatics and underweights the
intramolecular hydrogen bonding of a free catechol. Those residuals are **the published method's
error, not an implementation defect** - and they are why a Crippen logP is a `Predicted` quantity
that must be rendered with the method's RMS attached, never a value to quote to two decimals.

Consequences that remain, and are the project's rule:

1. A logP from this engine is only as good as the method: about 0.43 RMS on the tested set and
   about 0.67 by the method's own published benchmark. Efficiency metrics that consume logP or
   TPSA (BEI, SEI, LLE, LELP) must state which implementation produced their inputs.
2. Do not report the number as "AlogP". It is Wildman-Crippen as implemented here, with the
   parameter table transcribed from RDKit's BSD-3-Clause `Data/Crippen.txt` - an obligation
   recorded in `NOTICE`, section 4a.

## 7. What the engine still cannot do

| Gap | Evidence | Consequence |
|---|---|---|
| **No stereochemistry, anywhere** | `chem::Molecule` (`src/chem/Molecule.h`) has no field for it: `Atom` carries z, charge, H counts, aromatic and ring flags only, `Bond` carries order plus two flags | SMARTS `@`/`@@` is accepted and ignored (`Smarts.cpp:436`), `/` and `\` are plain single bonds, and a stereoisomer pair is **indistinguishable**. Measured: `C[C@H](N)C(=O)O` and `C[C@@H](N)C(=O)O` produce the identical canonical string `OC(=O)C(C)N` and therefore the identical `graphHash`. Any cache, dedup or library keyed on that hash will treat enantiomers and diastereomers as one compound |
| **No tautomer handling** | nothing in `src/chem` enumerates or standardises tautomers; the canonical writer preserves the graph it is given | Two tautomers of the same substance are two different molecules with two different hashes, two different descriptor sets and two different SMARTS match sets. Whichever tautomer a pack author typed is the one that is analysed |
| **No protonation / pH model in the graph** | charge is whatever the SMILES said | logP is a neutral-species estimate; there is no logD without the ionization work described in the project plan |
| **3D geometry is distance geometry, and is not a protein method** | `src/chem/Embed3D.h` describes it as "pragmatic distance geometry": 1-2 and 1-3 ideal distances from covalent radii and VSEPR angles, classical-MDS embedding via Eigen, then short steepest-descent relaxation on a bonded + angle + van-der-Waals term | It produces one plausible conformer for a viewer or a docking ligand-prep. It is **not** a conformational-ensemble generator, it does not minimise a real force field, and it must never be fed to a protein metric (TM-score, lDDT, GDT): a small-molecule distance-geometry embedding is not a protein structure. `embed3D` is deterministic given its seed and never emits NaN/Inf, which is a robustness guarantee, not an accuracy one |
| **No isotopes** | `Smarts.cpp:377` rejects isotope specifications outright because the graph stores none | No exact-mass work on labelled compounds; `molecularWeight()` is average mass only |
| **Fingerprints are BioCAD's own** | `morganFingerprint()` (`src/chem/Analysis.cpp:173`) is an ECFP-like circular hash to sorted unique 32-bit features, with `tanimoto()` over them | Similarity values are internally consistent and are **not** comparable to ECFP4/Morgan values from another toolkit, because the feature hashing differs |

## Measured evidence

Every number above came from throwaway programs compiled against the real sources and run on
this machine (AMD Ryzen 5 7600X3D, `g++ -O2`). Harness one - descriptors, aromaticity, canonical
form, SMARTS semantics and timing:

```
g++ -std=c++20 -O2 -Isrc /tmp/chemprobe.cpp \
    src/chem/Smiles.cpp src/chem/Rings.cpp src/chem/Aromaticity.cpp \
    src/chem/Canonical.cpp src/chem/Descriptors.cpp src/chem/Smarts.cpp \
    src/chem/Periodic.cpp src/chem/Crippen.cpp src/core/Assets.cpp \
    -I<nlohmann-json single_include> -o /tmp/chemprobe && /tmp/chemprobe
```

Output, verbatim:

```
caffeine   tpsa= 60.260 logP= -1.029 mw= 194.194 rings=2 arom=9
ibuprofen  tpsa= 37.300 logP=  3.073 mw= 206.285 rings=1 arom=6
dopamine   tpsa= 66.480 logP=  0.599 mw= 153.181 rings=1 arom=6
aspirin    tpsa= 63.600 logP=  1.310 mw= 180.159 rings=1 arom=6
benzene    tpsa=  0.000 logP=  1.687 mw=  78.114 rings=1 arom=6
benzene arom lower=6 kekule=6
naphthalene arom lower=10 kekule=10 rings=2
azulene aromatic atoms=0 rings=2
2-pyridone aromatic atoms=6
canon CC(=O)Oc1ccccc1C(=O)O          -> c1cc(c(OC(=O)C)cc1)C(=O)O  hash=eebf0e0fbc923ce3
canon O=C(O)c1ccccc1OC(C)=O          -> c1cc(c(OC(=O)C)cc1)C(=O)O  hash=eebf0e0fbc923ce3
canon c1ccc(OC(C)=O)c(C(O)=O)c1      -> c1cc(c(OC(=O)C)cc1)C(=O)O  hash=eebf0e0fbc923ce3
alanine L/R canonical equal=1  (OC(=O)C(C)N)
benzene ring matches=1 naphthalene=2
chiral smarts parsed=1 err=
isotope smarts parsed=0 err=isotope specifications are not supported (the graph stores no isotopes) at offset 1 in "[13C]"
hybridisation parsed=0 err=hybridisation primitive ^n is not supported at offset 2 in "[C^3]"
recursive parsed=1 err=
matchAll carboxylic-acid on ibuprofen: 0.745 us/call (hits=1)
recursive variant: 1.609 us/call (hits=1)
ALL PROBE ASSERTS PASSED
```

(An earlier run of the same harness measured 0.728 us and 1.477 us for the two timing lines; the
per-call cost is stable at roughly 0.7-0.8 us for the plain pattern and 1.5-1.6 us for the
recursive spelling.)

Harness two - the Crippen pack and the rule packs. Both spellings of caffeine now agree on logP
*and* on TPSA, because `crippen()` and `tpsa()` each prepare their own copy:

```
caffeine         ok=1 logP=  -1.029 MR=  51.20 ref= -0.07 dev=-0.959 tpsa= 60.260 note=assets/packs/descriptors/crippen.json
caffeine-kekule  ok=1 logP=  -1.029 MR=  51.20 ref= -0.07 dev=-0.959 tpsa= 60.260 note=assets/packs/descriptors/crippen.json
ibuprofen        ok=1 logP=   3.073 MR=  61.03 ref=  3.97 dev=-0.897 tpsa= 37.300 note=assets/packs/descriptors/crippen.json
dopamine         ok=1 logP=   0.599 MR=  42.53 ref= -0.98 dev= 1.579 tpsa= 66.480 note=assets/packs/descriptors/crippen.json
aspirin          ok=1 logP=   1.310 MR=  44.71 ref=  1.31 dev= 0.000 tpsa= 63.600 note=assets/packs/descriptors/crippen.json
benzene          ok=1 logP=   1.687 MR=  26.44 ref=  1.69 dev=-0.003 tpsa=  0.000 note=assets/packs/descriptors/crippen.json
citation: Wildman-Crippen (J Chem Inf Comput Sci 1999;39:868-873), method RMS ~0.67 log units
alert rules loaded=15 errors=0  groupPackErrors=0
paracetamol alert: para-aminophenol | para-Aminophenol / anilide (quinone-imine former) | Stepan et al. 2011, Chem Res Toxicol 24:1345-1410, doi:10.1021/tx200168d | warn=1 atoms=8
paracetamol groups: phenol=1 anilide=1 amide=1 aromaticRing=1
```

(The `ref` column there is experimental logP, which is why the deviations are large and why
section 6 treats them as the method's error rather than the implementation's. An earlier run of
this same harness, before `tpsa()` began perceiving, printed `tpsa= 56.220` on the
`caffeine-kekule` row.)

Harness three - the ten-compound accuracy set of `tests/test_chem.cpp`, aggregated:

```
amphetamine            ref= 1.76 got=  1.576 dev=-0.184
methamphetamine        ref= 2.07 got=  1.837 dev=-0.233
mdma                   ref= 2.15 got=  1.566 dev=-0.584
mda                    ref= 1.64 got=  1.305 dev=-0.335
phenethylamine         ref= 1.41 got=  1.188 dev=-0.222
tyramine               ref= 0.86 got=  0.893 dev=+0.033
ephedrine              ref= 1.13 got=  1.328 dev=+0.198
cathinone              ref= 0.59 got=  1.217 dev=+0.627
acetaminophen          ref= 0.46 got=  1.351 dev=+0.891
4-fluoroamphetamine    ref= 1.90 got=  1.715 dev=-0.185
n=10 MAE=0.349 RMS=0.430
```

Harness four - spelling independence of `tpsa()` and of the implicit-hydrogen rule, each compound
parsed twice (aromatic spelling, then Kekule):

```
caffeine   lower= 60.260 kekule= 60.260 diff=0.0e+00
aspirin    lower= 63.600 kekule= 63.600 diff=7.1e-15
furan      lower= 13.140 kekule= 13.140 diff=0.0e+00
pyridine   lower= 12.890 kekule= 12.890 diff=0.0e+00
thiophene  lower=  0.000 kekule=  0.000 diff=0.0e+00
CC(N)Cc1ccccc1                 tpsa= 26.020 logP=  1.576
NCCc1ccc(O)c(O)c1              tpsa= 66.480 logP=  0.599
CC(=O)Nc1ccc(O)cc1             tpsa= 49.330 logP=  1.351

c1ccsc1        formula=C4H4S      mw=  84.136
C1=CSC=C1      formula=C4H4S      mw=  84.136
c1cc[nH]c1     formula=C4H5N      mw=  67.091
c1ccccc1       formula=C6H6       mw=  78.114
C1=CC=CC=C1    formula=C6H6       mw=  78.114
c1ccncc1       formula=C5H5N      mw=  79.102
```

Thiophene's 0.000 is the N,O-only TPSA convention, not a bug (section 6); the aspirin 7.1e-15
difference is floating-point summation order, not a perception difference.

Harness five - the RDKit oracle comparison, reproduced for this page. The 69 pack SMILES were
extracted from `assets/packs/*.json`, run through the C++ engine
(`crippen()`, `tpsa()`, `tpsaIncludingSulfurAndPhosphorus()`, `molecularFormula()`) into a TSV,
then compared in Python against `rdkit.Chem.Crippen.MolLogP` / `MolMR`,
`Descriptors.TPSA(m)` / `Descriptors.TPSA(m, includeSandP=True)` and
`rdMolDescriptors.CalcMolFormula`. RDKit is an oracle only: it is not linked, not in
`vcpkg.json`, and no BioCAD build depends on it.

```
rdkit 2026.03.5

n=69  logP mismatches=0  MR=0  formula=0
TPSA (N,O only) mismatches>0.1: 10
    caffeine -1.56
    phenazone -1.04
    ondansetron -1.04
    theobromine -1.04
    theophylline -1.04
    celecoxib -0.52
    indomethacin -0.52
    nmn 0.22
    thiamine 0.22
    berberine 0.22
TPSA (incl S,P) mismatches>0.1: 11
    n-acetylcysteine -13.5
    nmn -9.59
    caffeine -1.56
    phenazone -1.04
    ondansetron -1.04
    theobromine -1.04
```

Formula comparison strips RDKit's charge suffix before comparing, because
`molecularFormula()` reports the Hill formula without a charge annotation.

Related reading: [provenance.md](provenance.md) for the tier a descriptor may carry,
[limitations.md](limitations.md) for the project-wide honesty rules, and
[packs.md](packs.md) for how rule and descriptor data reach the engine.
