# DNA / RNA workbench

The `NucleicAcid` panel ("DNA / RNA Workbench", group `Workspace`) and the
`INucleicAcidModule` contract behind it cover sequence handling that a bench
scientist actually does by hand: read a FASTA or GenBank file, look at its
features, cut it with restriction enzymes, read the six frames, find the ORFs,
compute an oligo's melting temperature, design a primer pair, measure codon usage,
and enumerate CRISPR protospacers with an honest off-target scope.

Implementation lives in `src/bio/NucSeq.*`, `src/bio/NucIo.*`,
`src/bio/Restriction.*`, `src/bio/OligoThermo.*` and `src/bio/Codon.*` (sequence
arithmetic and thermodynamics), composed by `src/modules/NucleicModule.*` (the
adapter, which owns primer design and the guide search outright).

## Biosecurity boundary

This boundary is permanent and is stated identically in `src/data/Nucleic.h`,
`src/modules/NucleicModule.h`, the panel, the agent tools and here:

- **No synthesis-vendor integration.** There is no vendor client, no vendor
  sequence format, and no code path that transmits a sequence anywhere.
- **No order-sheet export.** Export is **FASTA and GenBank only**. There is no
  plate map, no oligo order table, no CSV of constructs to buy.
- **No pathogen-driven batch gene design.** The workbench operates on one record
  the user pasted; there is no batch designer and no pathogen-derived workflow.
- **No therapeutic or germline CRISPR framing.** Guide search is sequence
  analysis. It emits no protocol, no delivery method, no clinical framing.
- **Codon optimization is constraint satisfaction, never an expression
  prediction.** `CodonOptimizationResult` deliberately has no yield, titre or
  expression field to put one in. The only guarantees it makes are that the output
  translates to exactly the input protein and contains none of the forbidden
  patterns the user named.

## Tm parameter set and salt convention

- Nearest-neighbour parameters are the SantaLucia unified set, loaded from
  `assets/packs/nucleic/nn-thermodynamics.json` with its citation. The parameter
  set name travels on every `Quantity::source`, so a Tm can never be quoted
  without it.
- `OligoThermo` echoes back the **exact** conditions used: `naMolar`, `mgMolar`,
  `oligoMolar`, `dntpMolar`. A Tm without its salt and strand concentration is not
  reproducible at a bench, so the panel prints those four numbers on the same
  screen as the Tm and the agent tool returns them in the same object.
- **The monovalent salt correction applies to Tm only.** `deltaH`, `deltaS` and
  `deltaG37` are the 1 M Na+ standard-state values.
- The correction that is applied is SantaLucia's entropy-per-phosphate form,
  `dS([Na+]) = dS(1 M) + 0.368 * (N - 1) * ln([Na+] / 1 M)` with `N` the number of
  nucleotides, and it appears verbatim in the returned `Quantity::source`.
- Defaults, used when a concentration arrives as zero: **0.05 M Na+** and
  **5e-8 M total strand concentration Ct** (total, not per-strand).
- **`mgMolar` is echoed but not applied.** BioCAD assumes no divalent-to-monovalent
  equivalence, because the published equivalences disagree with each other by more
  than the effect they model. `OligoThermo::assumptions` says so explicitly, and
  the panel renders that list next to the Tm.
- **Hairpin, self-dimer and hetero-dimer `deltaG37` carry no salt correction
  at all** and are the 1 M Na+ standard-state values. The published monovalent
  correction is defined for Tm, not for the free energy of a looped or mismatched
  structure; applying it there would be an invention. `selfStructures()` still
  accepts the salt concentration its caller has, and documents that it does not use
  it, rather than hiding the argument and inviting a second entry point later.
- Zero for any concentration means "the caller did not say", and the documented
  default stands in and is reported back on the result - never a Tm computed at
  0 M salt.

## Primer design

`designPrimers(record, begin, end, targetTmC)` enumerates candidates whose 5' ends
walk outward from the requested interval, so **the product always contains the
interval**: the interval is the contract, the flank is the freedom.

Every limit below is a **hard rejection, not a weight**. A pair that violates one
is absent from the result at any rank, because ranking a primer that will dimerise
with itself is worse than returning nothing. The defaults are in
`PrimerDesignLimits`:

| Limit | Default | Why |
|---|---|---|
| length | 18-30 nt | conventional working range |
| GC | 40-60% | outside it, Tm and specificity both drift |
| Tm window | target +/- 4 degC | each primer individually |
| \|Tm(F) - Tm(R)\| | < 2.0 degC | one primer annealing while the other does not is not fixable by ranking |
| hairpin / self-dimer dG37 | > -6.0 kcal/mol | roughly a 3-4 bp duplex; more stable competes with the template |
| cross-dimer dG37 | > -8.0 kcal/mol | a cross-dimer costs both primers at once |
| homopolymer run | <= 4 | five identical bases slips |
| 3'-terminal G or C | required | a GC clamp anchors extension |
| G/C in the last 5 | <= 3 | a GC-saturated 3' end primes mismatched templates |

Ranking among survivors uses Tm match, |dTm|, 3'-end stability (the dG37 of the
terminal pentamer, scored by distance from -7.5 kcal/mol - too weak and the
polymerase has nothing to hold, too strong and it extends off a mismatch) and the
worst remaining structure. Cross-dimer evaluation is the expensive step, so only
the best 16 candidates per side are paired: 256 hetero-dimer evaluations rather
than tens of thousands.

Every returned `PrimerPair` carries its `tmDifference`, its liabilities sorted
most-stable-first, and a warning listing the thresholds that were applied - so a
reader knows an empty list means "no acceptable pair exists", not "the search
failed".

## Guide search, and the honesty rule

`findGuides(target, reference, pam)` enumerates PAM-adjacent protospacers in the
target (`pam` is an argument, matched under IUPAC, default `NGG`) and counts
off-targets in **the reference the caller actually supplied** - which may be a
plasmid, a contig, or a genome. Circular records are linearised across the origin
so a site straddling position 0 is not silently dropped.

**Algorithm.** Pigeonhole seed index, not a naive O(n*m) scan. Any alignment of a
20-mer with at most 2 mismatches leaves at least one of three disjoint blocks
(7/7/6 nt) completely intact, because two mismatches cannot hit three blocks. The
reference's PAM-adjacent sites are indexed once by (block index, block sequence)
into a 64-bit key at 2 bits per base; a query looks up its own three blocks, unions
and dedupes the candidate positions, and only those are compared base-by-base with
an early exit past 2 mismatches. Cost is O(reference) once plus O(candidates) per
guide. The guide's own on-target site is located in reference coordinates and
excluded, so a guide never counts itself as its own off-target; when the target
does not occur in the reference at all, that is reported as a warning rather than
papered over.

**The scope rule.** `basesSearched` and `referenceName` are mandatory fields, and
`scopeStatement` is a sentence a user can read. No surface may render an off-target
count without them: the panel prints the scope statement *above* the guide table
rather than in a tooltip, and the `search_guides` agent tool refuses genome-wide
specificity questions outright and returns the scope statement alongside every
count.

Guide enumeration is capped (200 by default, `GuideSearchLimits::maxGuides`), and
when the cap is reached the scope statement says so: a truncated list is a scope
limit like any other, and leaving the reader to infer it from a round number would
be the same failure as reporting a count without its reference.

**`genomeWideClaimPossible`** is false by default and true only when both:

1. the reference is not the target itself (different id **and** strictly longer),
   and
2. `basesSearched >= 1,000,000`.

That criterion is deliberately conservative, and it is **necessary, not
sufficient**. The smallest free-living bacterial genomes are around 0.6 Mb, so a
1 Mb floor already excludes every plasmid, amplicon and contig anyone will paste
in - but a 1 Mb contig passes the floor and is still not a genome, and BioCAD
cannot verify that any input is complete. So even when the flag is true the scope
statement says the reference *could* be a genome but completeness is unverifiable,
and the prose statement, not the flag, remains authoritative.

### Measured numbers from a real run

From `/tmp/nm.cpp` against this implementation (the same assertions as
`tests/test_nucleic_module.cpp`):

- The real pUC19 (`tests/fixtures/pUC19.gb`, L09137.2, 2686 bp, circular, 4 features)
  searched against itself with `NGG`: `basesSearched == 2686`,
  `genomeWideClaimPossible == false`, the enumeration cap of 200 protospacers
  reached, and the scope statement reading in full:
  *"Searched 2686 bases of 'pUC19' (circular), both strands, for NGG-adjacent
  20-mer sites with up to 2 mismatches. Off-target counts below are counts WITHIN
  those 2686 bases only; the guide's own on-target site is excluded. No sequence
  outside this reference was examined, and this reference is far too small to be a
  genome, so these counts are NOT a genome-wide specificity claim. Guide
  enumeration stopped at the first 200 protospacers in the target, so the list
  below is a prefix and not every candidate site; the off-target counts for the
  guides shown are unaffected."*
  A truncated list is itself a scope limit, so it is stated in the scope statement
  rather than left for the reader to infer from a round number.
- The same PAM machinery with `TTTV` (Cas12a) over pUC19: 76 guides, every PAM
  `TTT[ACG]`, every protospacer 20 nt.
- The first 400 bp of pUC19 as the target, against a reference carrying one planted
  exact copy of guide `CCTCTGACACATGCAGCTCC` + `CGG` and one planted 2-mismatch
  copy (476 bases searched): `exactOffTargets == 1`,
  `oneMismatchOffTargets == 0`, `twoMismatchOffTargets == 1` - the on-target itself
  excluded, the planted 2-mismatch site counted in its own bucket and not as an
  exact hit.
- A 500 bp requested product on a 900 bp template at a 60 degC target Tm yields
  5 pairs. Best pair: `TTGCCAAGGTCATGCGATCAGCTTG` (Tm 60.28 degC, GC 52.0%) with
  `GCAATCGATGCCTAGCTACGATGGC` (Tm 60.12 degC, GC 56.0%), |dTm| 0.166 degC, product
  100-600 (500 bp), worst structure dG37 -4.48 kcal/mol. Across all five pairs
  |dTm| <= 0.721 degC, every GC in 50.0-56.0%, and no structure at or below
  -6.0 kcal/mol.
- Oligo `GTAAAACGACGGCCAGT` at 0.05 M Na+ and 2.5e-7 M strand: Tm 51.67 degC
  (`Predicted`, unit `degC`), dG37 -21.73 kcal/mol, GC 52.9%. With every
  concentration left at zero the documented defaults stand in and are reported
  back: 0.05 M Na+, 5e-8 M total strand (Ct), Tm 49.17 degC.

## Agent tools

`parse_sequence`, `translate_sequence`, `restriction_map`,
`oligo_thermodynamics`, `design_primers`, `optimize_codons`, `search_guides`.

Each tool description carries the discipline the number needs: `optimize_codons`
forbids any expression claim, `oligo_thermodynamics` requires the conditions to be
quoted with the Tm, `design_primers` states that violations are rejected rather
than ranked, and `search_guides` refuses the genome-wide question.
