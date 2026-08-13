# Variant analysis and point-mutation modelling

The least certain area BioCAD touches, and the one where a plausible-looking
number does the most damage. Everything here is shaped so the uncertainty cannot
be dropped on the way to the screen:

- a conservation score carries the homolog **count** and the **median identity**
  of the alignment it came from, and is refused outright below a minimum;
- a rebuilt side chain is `Provenance::Model` - a constructed artefact with **no
  energy claim**, and `RotamerRebuild` has no field to put one in;
- there is **no pathogenicity verdict, no clinical significance, no binder
  design**, and those absences are deliberate.

Contract: `IVariantModule` (`src/contracts/IModules.h`).
DTOs: `src/data/Variants.h`. Implementation: `src/bio/Conservation.*`,
`src/bio/Rotamer.*`, `src/modules/VariantModule.*`. Panel: `Variants`
("Variant Analysis", group `Discover`). Agent tools: `conservation_profile`,
`score_variant`, `rebuild_side_chain`.

## 1. Conservation

### The homolog minimum, and why there is one

`ConservationProfile::usable` is **false below 15 supplied homologs**
(`bio::kMinimumHomologs`), and then the profile has no columns at all and every
dependent `Quantity` on a `VariantScore` is `notComputed("at least 15 homologs
(N supplied)")`. The panel prints that sentence instead of drawing a track.

A conservation score is an estimate of a column's amino-acid distribution. With
five sequences there are at most five distinct observations, the pseudocount
dominates the posterior, and the resulting SIFT ratio is largely a function of
the pseudocount scheme rather than of the protein - yet it gets read against a
published 0.05 threshold as if it were the quantity SIFT computes over dozens to
hundreds of PSI-BLAST hits. PROVEAN's own threshold was calibrated over clusters
of up to 45 supporting sequences per query.

**15 is a BioCAD judgement, not a published constant**, and the profile's own
warning text says so. It is a floor, not a sufficiency test: 15 near-identical
sequences are still a bad set, which is why the effective sequence count after
identity clustering is displayed next to the raw count.

### What is computed

| Number | How | Tier |
|--------|-----|------|
| Alignment | `bio::alignGlobal` - the tree's one Gotoh affine-gap aligner, BLOSUM62, gap open 11 extend 1. Each homolog is aligned to the query and projected onto QUERY columns | - |
| Sequence weighting | homologs are single-linkage clustered at 80% identity and each member carries `1/clusterSize`; the sum is `effectiveSequenceCount` | - |
| `shannonEntropy` | entropy in bits of the **observed** weighted column distribution, no pseudocount. A perfectly conserved column is exactly 0; a column holding each of the 20 residues once is exactly log2(20) = 4.321928094887362 | Measured |
| `frequencies` | the pseudocounted posterior (Dirichlet 0.05 per residue), which is what a ratio and a log-odds must be computed from | - |
| `pssm` | log2(frequency / background) against Robinson & Robinson 1991 (PNAS 88:8880-8884), the composition NCBI BLAST pairs with BLOSUM62. The set used is named in `backgroundFrequencySource` | - |
| `blosum62Delta` | `score(wt, mut) - score(wt, wt)`, an exact table lookup difference | **Measured** |
| `siftScore` | `p(mutant) / max_y p(y)` over the column posterior. Exactly 1.0 when the mutant is the column consensus. Deleterious below **0.05** (Ng & Henikoff 2003) | Heuristic, no unit |
| `proveanScore` | the change in summed alignment score against the supporting set. Under a fixed alignment a single substitution telescopes exactly to `sum_y f_y (score(mut,y) - score(wt,y))` - a gap column contributes zero because the gap penalty does not depend on the query residue. Deleterious below **-2.282** (Choi et al. 2012, PLoS ONE 7:e46688; balanced accuracy 79.05% on their 58,684-variant human validation set) | Heuristic, no unit |

SIFT and PROVEAN are `Provenance::Heuristic` with an **empty unit**, which
`makeQuantity()` enforces: they rank, they are dimensionless, and no published
model of ours ran to produce them. They are *style-of* implementations - the
published pipelines differ in how they gather and weight sequences, so a number
here is not interchangeable with a number from the SIFT or PROVEAN servers, and
the `source` string on every Quantity says which set it came from.

### What this is not

`VariantScore::interpretation` states, in the DTO itself, that the numbers
describe how unusual the residue would be in that alignment, and that they are
**not a pathogenicity call, not a clinical interpretation, and not a verdict
about a person**. There is no ACMG classification, no ClinVar lookup, and no
"likely pathogenic" string anywhere in this subsystem.

## 2. Point-mutation rebuild

### The rotamer library: what ships, and what does not

The Dunbrack 2010 smoothed backbone-dependent rotamer library is CC BY 4.0 and
*would* be redistributable in a commercial product with attribution. It is not
here, for a reason that has nothing to do with the licence: the only distribution
channel is a download page issued after a **manually approved licence
application** (`https://dunbrack.fccc.edu/lab/bbdep2010`, "Example of Download
Page"), and the historical direct URLs under `dunbrack.fccc.edu/bbdep2010/` return
404. It could not be obtained offline, and inventing rotamer angles instead is not
an option.

**What ships instead** is `assets/packs/rotamers/rotamers-pdb-derived-2026.json`,
generated by `scripts/build-rotamer-pack.py`, in which every number is measured:

- **Dataset**: 219 RCSB PDB X-ray entries at <= 1.4 A resolution, 150-600 polymer
  residues, occupancy >= 0.5, altLoc blank or A. 47,831 side chains. The entry ids
  are listed in the pack's `dataset.entries`. Coordinates are CC0.
- **Included**: for each of the 18 residue types with a chi angle, the rotamer
  populations with their counts, probabilities, circular mean chi and circular
  standard deviation, both backbone-independent and binned at **60 degrees** in
  phi and psi; plus the side-chain internal coordinates (mean bond length, mean
  bond angle, and either the rigid torsion or the offset from the chi that turns
  the atom) needed to build the side chain.
- **Not included**: the Dunbrack library's adaptive kernel density estimates and
  its 10-degree phi/psi resolution. A 60-degree bin with fewer than 20
  observations is **not emitted at all**; the consumer falls back to the
  backbone-independent set for that residue and records it as an assumption on
  the result. Nothing is smoothed, interpolated or extrapolated into an empty bin.
- **Not included**: proline. Its side chain closes a ring onto the backbone
  nitrogen and constrains phi, so a rebuild to proline is **refused** with that
  reason rather than approximated. Glycine has no side chain and returns nothing
  to build, which is stated rather than treated as an error.
- The pack states all of this in its own `method` and `notDunbrack` fields, and
  `NOTICE` records both the CC0 coordinate provenance and the fact that the
  Dunbrack library is cited for its binning convention and **not redistributed**.

Measured accuracy of the construction: rebuilding all 285 complete side chains of
the 1VFB fixture from **their own** chi angles reproduces the deposited heavy-atom
positions with a median deviation of 0.158 A, a 95th percentile of 0.587 A, and a
worst case of 2.09 A (an arginine NH, five bonds out, where mean bond angles
accumulate). That residual is exactly why a rebuilt side chain is `Model` and not
`Measured`. The assertion is in `tests/test_bio_rotamer.cpp`.

### Selection: Goldstein dead-end elimination

The mutated position plus every neighbouring side chain whose CB is within
**6.0 A** (CA for glycine), nearest first, **capped at 8** neighbours, are made
flexible; the backbone is held fixed. Rotamer r at position i is eliminated when
some competitor t is better regardless of what the other positions do:

```
self[i][r] - self[i][t] + sum_{j != i} min_s ( pair[i][j][r][s] - pair[i][j][t][s] ) > 0
```

(Goldstein 1994, Biophys J 66:1335-1340.) Elimination repeats until a pass removes
nothing, then the survivors are repacked exhaustively when the remaining space is
at most 200,000 combinations and by iterative single-position descent otherwise -
which the result reports, because a solver that hides which one it ran is not
trustworthy.

**The selection objective is unit-free and is never reported.** It is the
heavy-atom clash count plus `-ln(rotamer probability)`; a clash is a heavy-atom
pair closer than 0.75 x the sum of their van der Waals radii, with hydrogens,
waters, ligands and the bonded backbone atoms of the packed residue and its
sequence neighbours excluded. There is no solvation term, no electrostatics, no
torsional term and no force field. `RotamerRebuild` has **no energy field**, and
`tests/test_variant_module.cpp` asserts at compile time that it never gains one.

## 3. Stability prediction (ddG) - not shipped, and why

`StabilityPrediction` exists in `src/data/Variants.h` and is **not produced by
any code path**. Nothing in this build predicts a ddG.

The tier is real and would be legal: **ThermoMPNN is MIT-licensed, 2.7M
parameters and CPU-feasible**, and ONNX Runtime is in vcpkg. What is missing is
concrete: **the model weights are not in this tree**, and there is no honest way
to produce a ddG without a model actually running. A regression invented here to
fill the gap would be exactly the failure this whole document is written against.

If that tier is ever added, the DTO already forces the two numbers that make it
readable, and neither is optional:

- the benchmark RMSE goes **in the value**, so a ddG renders as
  `-1.4 +/- 1.55 kcal/mol (ThermoMPNN, Fireprot-HF RMSE)` rather than as `-1.4`;
- `positivePredictiveValuePct` is a required field, and for the current best
  CPU-feasible models the positive predictive value for **stabilising** mutations
  is **46-56%** - roughly a coin flip - which the panel must print beside the
  number rather than in a footnote.

Until then the panel says, in place of a value: *"Stability change (ddG): not
computed. No model weights ship in this build."*

## 4. Not shipped, with reasons

Reproduced verbatim from [limitations.md](limitations.md) section 4, which is the
single list this section must not drift from:

| Not shipped | Reason |
|-------------|--------|
| Rosetta | Cannot be redistributed, and serving it to third parties is prohibited |
| FoldX | Commercial fee *and* worst-in-class: CAPRI Kendall 0.12, MDM2-p53 rho -0.14 |
| AlphaFold3 weights | Non-commercial only, and must not be shared |
| SignalP 6.0 | DTU commercial licence |
| Chai-1 | Apache-2.0, but Linux-only, bfloat16, and 48-80 GB VRAM |
| ColabFold / AF2 locally | 940 GB MSA database, and the public server's terms forbid fan-out from a distributed app |
| ESMFold public API | Unreliable, and capped at 400 residues |
| mCSM-PPI2 | Web-server-only; BioCAD links out instead of pretending to host it |
| RFdiffusion binder design | The authors' own pipeline is ~10,000 backbones to ~20,000 ProteinMPNN sequences to an AF2 `pae_interaction < 10` filter to wet lab. Rendering one backbone as "a designed binder" misrepresents a 10,000-to-1 funnel |
| De novo enzyme / catalytic design | No honest success rate to report from a desktop tool |
| Epitope / immunogenicity prediction | Linear B-cell prediction is around AUC 0.6. Propensity plots only, never a call |

## 5. The panel

`Variants` / "Variant Analysis", group `Discover`. Three sections, in this order,
because the order is the argument:

1. **Query and homologs**, with the supplied count, the required minimum, the
   effective count after clustering, and the median / min / max identity as stat
   cards - beside, not behind, everything below them.
2. **Conservation track**: per-column entropy on a fixed 0 to log2(20) axis so two
   proteins are comparable and a flat track cannot come from autoscaling. **Below
   the minimum the track is not drawn at all** - a greyed-out plot invites the
   reader to squint at it; a sentence reading "4 homolog(s) supplied, 15 required"
   does not.
3. **Substitution**, with each score's published threshold printed next to it, and
   **Rebuild**, badged `model`, listing the library, the chi angles, the library
   probability, the clash count, the repacked neighbours and every assumption -
   ending with the explicit `not computed` line for ddG.

## 6. Tests

| File | Asserts |
|------|---------|
| `tests/test_bio_conservation.cpp` | the shortfall refusal names its numbers; a conserved column is exactly 0 bits and a uniform one exactly log2(20) to 1e-12; SIFT of the consensus is exactly 1.0; ten BLOSUM62 deltas against direct lookups; the PROVEAN delta against hand arithmetic; the background sums to 1; near-duplicates collapse to effective count 1 |
| `tests/test_bio_rotamer.cpp` | the pack loads and declares it is not Dunbrack; LEU's dominant rotamer is `mt` with chi1 < 0 (a torsion sign error reads `pt`); a helical bin is used and a sparse one falls back; `placeAtom`/`dihedralDegrees` are exact inverses; 285 real side chains rebuild to median 0.158 A; the hand-built DEE case eliminates exactly the dominated rotamer; 200 random problems never eliminate a rotamer in the brute-force optimum |
| `tests/test_variant_module.cpp` | every dependent Quantity is `NotComputed` below the minimum while the BLOSUM62 lookup survives; units and tiers; proline, unknown residues and missing chains are refused by name; a rebuild is `Model`, names its library, and `RotamerRebuild` has no `energy`, `score` or `deltaDeltaG` member (compile-time); `Services::valid()` requires the module |
