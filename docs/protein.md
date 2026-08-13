# The protein core

Phase 5 added a `bio::` namespace beside `chem::`: structure and sequence I/O, alignment,
superposition and structure comparison. **No external library was taken.** `gemmi` and
`chemfiles` are not in this vcpkg registry (verified with three independent 404s each against
a passing `eigen3` control), `OpenStructure` is LGPL-3.0 and `ESBTL` is GPL. Eigen was already
linked and is the only maths dependency.

| Path | Contents |
| --- | --- |
| `src/bio/Structure.{h,cpp}` | Model -> chain -> residue -> atom, plus `sequenceOf` |
| `src/bio/Annotations.h` | SEQRES / HELIX / SHEET / SSBOND, held out of band |
| `src/bio/PdbReader.{h,cpp}` | Fixed-column PDB v3.3 |
| `src/bio/CifReader.{h,cpp}` | Minimal STAR/mmCIF parser |
| `src/bio/Fasta.{h,cpp}` | FASTA with the UniProt and RCSB header grammars |
| `src/bio/Align.{h,cpp}` | Gotoh affine-gap DP, BLOSUM matrices, Karlin-Altschul E-values |
| `src/bio/Superpose.{h,cpp}` | Kabsch with the reflection correction |
| `src/bio/Score.{h,cpp}` | lDDT and Shrake-Rupley SASA |
| `src/modules/BioModules.{h,cpp}` | `RealSequence` and `RealStructure`, thin adapters over `bio::` |

## Coordinates are not conformers

`chem::Conformer` is the output of a small-molecule distance-geometry embedder. It is not a
protein, and a TM-score or lDDT computed from one is meaningless.

**This is enforced at compile time, not by a comment:** `bio::Structure` has no constructor,
conversion or overload taking a `chem::Conformer`, and `Superpose.{h,cpp}` depends on neither
`bio::Structure` nor `chem::Conformer` - it is pure `Point3` geometry. There is no expression
that feeds embedded small-molecule coordinates into a structure-comparison number.

## Why the mmCIF reader is mandatory

The PDB file format has been **frozen at v3.30 since 21 November 2012**. No PDB file is issued
for an entry with more than 62 chains, more than 99999 atoms, multi-character chain ids, or a
B-factor above 999.99. A PDB-only reader therefore cannot open modern structures, and fails by
silently not finding a file rather than by reporting anything - which is the worst failure mode
available. Both readers are first-class, and the test fixtures describe the *same* entry in
both formats so they can be diffed against each other.

The mmCIF reader resolves columns **by tag name only**. The committed fixture deliberately
uses a non-canonical column order, so a positional parser would fail the test.

## PDB column traps

Each of these is a real bug that a naive reader ships with:

| Trap | Handling |
| --- | --- |
| Element symbol | Columns 77-78 when present; otherwise inferred positionally, where a one-letter element's name starts at column 14 and a two-letter element's at column 13. `" CA "` is carbon-alpha; `"CA  "` is calcium. Getting this wrong turns a protein into a metal. |
| Negative `resSeq` | Legal, and round-trips. |
| Insertion codes | Part of residue **identity**: residues 100 and 100A are different residues, and pairing by number alone silently merges them. |
| Columns 67-76 | Undefined in v3.3 and never parsed. |
| `altLoc` | Retained, not collapsed. |
| Malformed lines | Appended to `Structure::warnings` and skipped. The reader never throws on a truncated file, and the UI shows the warnings. |

`auth_seq_id` and `label_seq_id` are both carried, because they disagree and papers cite the
author numbering. **Any UI that shows a residue number states which numbering it used** - the
Protein Structure panel prints it explicitly.

## Alignment

Gotoh affine-gap dynamic programming with the three matrices, global (Needleman-Wunsch) and
local (Smith-Waterman), with full traceback. BLOSUM62 loads from
`assets/packs/matrices/blosum62.json` (canonical NCBI, Henikoff & Henikoff 1992, half-bit
units, 24 symbols plus `*`), so BLOSUM45/80 need no rebuild. Defaults are gap open 11 and
extend 1, the BLASTP defaults for BLOSUM62, with the NCBI convention that a gap of length k
costs `open + k*extend`.

**Identity and similarity are two separately defined numbers**, both reported, both defined in
the header:

- identity = identical columns / aligned columns
- similarity = columns with a positive substitution score / aligned columns

Conflating them is the standard way an alignment is oversold. A K->R substitution scores +2 in
BLOSUM62: 91.7% identical, 100% similar, and the panel shows both.

**E-values are local-alignment only.** An E-value on a global alignment is meaningless, so
`GlobalAlignment` has no `evalue` member and `evalueOf` only overloads on `LocalAlignment` -
a compile-time distinction, not a runtime warning. The test contains a `static_assert` proving
`evalueOf` is not invocable with a global alignment. The Karlin-Altschul lambda/K/H/alpha/beta
parameters live in the matrix JSON and the `BLAST_ComputeLengthAdjustment` edge correction is
applied; `statisticsFor` returns null rather than a guessed lambda for an unparameterised gap
cost.

## Superposition

Kabsch via Eigen `JacobiSVD`, **including the determinant reflection correction**: if
`det(V * U^T) < 0`, the third column of V is negated. Without it, Kabsch cheerfully returns a
mirror image with a deceptively good RMSD - the single most common silent error in structural
superposition code. The test explicitly superposes a reflected point cloud and asserts the
RMSD is large; `Superposition::reflectionCorrected` surfaces that the correction fired.

Measured: a 37-degree rotation plus translation over 40 points is recovered with a fitted RMSD
of 2.2e-15, and a mirrored cloud gives RMSD 6.68 with `reflectionCorrected = true`.

## Structure comparison

**lDDT** is superposition-free: for every reference atom pair inside the 15 A inclusion radius,
a distance is preserved if it changes by less than each of 0.5, 1, 2 and 4 A, and the score is
the mean preservation fraction over the four tolerances. Because it never superposes, it
measures local quality and is insensitive to domain motions. `lDDT(X, X)` is **exactly 1.0** by
construction, and lDDT is invariant under a rigid rotation and translation of the whole model -
both are asserted, and both hold exactly in the measured run.

A consequence worth knowing: displacing one atom by 0.3 A leaves lDDT at exactly 1.0, because
0.3 A is below the smallest tolerance. That is the metric behaving correctly, not a bug.

**SASA** is Shrake-Rupley with 92 golden-spiral points and a 1.4 A probe, over an explicit
Bondi radii table, ignoring hydrogens (united-atom, because most crystal structures have none).
Relative accessibility uses the Tien et al. 2013 theoretical maxima, and a non-standard residue
gets `NotComputed` rather than a relative value against a guessed maximum.

**Every SASA readout carries its algorithm, probe radius, point count, radii set and hydrogen
policy** in `Quantity::source`. A bare SASA number is not reproducible: two tools will
routinely disagree by ten percent on the same structure, and the difference is entirely in
these parameters.

Residues are paired between structures by `(chain id, authSeqId, insertionCode)`, and the
unmatched count is returned so a caller can refuse a comparison that matched almost nothing.

**`StructureComparison::tmScore` is always `NotComputed("vendored TM-align")`.** TM-align is
planned as a vendored single permissively-licensed source file; until it is actually there,
the field reports that it was not computed rather than carrying a plausible substitute.

## Contracts, panels and tools

`ISequenceModule` (`alignGlobal`, `alignLocal`) and `IStructureModule` (`load`, `compare`,
`sasa`) join `Services`, implemented in `src/modules/BioModules.cpp` as thin adapters over the
pure `bio::` functions. There is one implementation, and the tests run it.

Two panels:

| Id | Title | Group |
| --- | --- | --- |
| `Sequence` | Sequence Compare | Discover |
| `Structure3D` | Protein Structure | Workspace |

The existing `Structure` id is the small-molecule Structure Workbench and was neither reused
nor renamed; panel ids are persisted in `imgui.ini` and are part of the app's API.

Agent tools: `align_sequences`, `fetch_structure` (honest that it reads a local file or a
cached download, and states the numbering scheme), and `compare_structures`.

## Verified numbers

Measured by an executed harness against the real sources, not asserted from a test file:

```text
pdb atoms=24  cif atoms=24  warnings=0/0     both readers agree on the fixture
global id=91.6667% sim=100.0000% score=65.0  evalue=not computed  (global has none)
kabsch rmsd(rotated)  = 2.174e-15  reflectionCorrected=0
kabsch rmsd(mirrored) = 6.676135   reflectionCorrected=1
lDDT(X,X)             = 1.000000000000000
lDDT(rigidly moved)   = 1.000000000000000
SASA(single C)        = 120.762822 A^2 vs analytic 4*pi*(1.70+1.40)^2 = 120.762822
compare: tmScore = not computed
```
