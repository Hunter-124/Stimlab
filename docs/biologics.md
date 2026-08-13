# Biologics and antibody engineering

The `Antibody Workbench` panel, `IBiologicsModule`, `src/bio/{Imgt,Liabilities,Interface}.*`
and the packs under `assets/packs/biologics/`.

Four things happen here, in order of how much they can be trusted: IMGT numbering
(exact, or refused), sequence liabilities and descriptors (cited motifs and exact
arithmetic), biologics masses (composition arithmetic over NIST SRD 144), and
interface geometry with a geometric alanine scan (measured geometry, plus one
deliberately unit-free rank ordering).

## The do-not-ship list, verbatim

**No humanness score and no humanization suggestion. No immunogenicity prediction.
No affinity-maturation suggestion. No aggregation free energy. No viscosity. No
expression titre. No shelf-life prediction.**

Their absence is the feature, not a gap. There is no field for any of them in
`src/data/Biologics.h`, no code path that produces one, and the agent tool
descriptions forbid improvising one from the numbers that are present. If you want
one of these, the honest answer is an experiment, not a number from a sequence.

Two more that follow from the same rule:

- **`closestGermlineSet` is never a species identification.** It is the best-scoring
  entry in a reference pack. A camelid VHH scores best against human `IGHV3-23`; that
  is a similarity result. The panel header says so and there is no species field to
  render.
- **The alanine scan is never a ddG.** See below.

## 1. Numbering

`bio::numberDomain()` numbers one V-DOMAIN in the IMGT unique numbering (Lefranc et
al., Dev Comp Immunol 27:55-77, 2003).

There is no Python runtime and no HMMER. What ANARCI consumes is IMGT's own reference
data, and IMGT publishes V-REGION amino-acid sequences **already gapped into the unique
numbering**: in `assets/packs/biologics/imgt-reference.json`, index *i* of a `gapped`
string *is* IMGT position *i*, and every functional entry carries Cys at index 104.
Position transfer is therefore a table lookup through one alignment, and the aligner is
Phase 5's Gotoh (`bio/Align.h`) - there is exactly one aligner in this tree.

How a query is numbered:

1. Local-align it against all 550 functional germline V alleles with BLOSUM62 and
   gap costs 11/1. The best and runner-up **gene sets** and their bit scores are
   reported.
2. Transfer the winning germline's IMGT columns onto the query to locate the anchors.
3. Locate IMGT 118 through the J reference: the J-PHE/J-TRP that opens FR4. The matched
   J also fixes FR4's length - 11 residues (118-128) for a heavy chain, 10 for a light
   one - which is what stops the first constant-domain residue being numbered as 128.
4. Number the frameworks from the germline profile and the CDRs from the IMGT layout in
   the pack: `deletionOrder` empties slots centre-outward as a loop shortens, and
   `insertionAnchors` carry insertion codes when it is longer.

### A failed anchor returns nothing

IMGT 23 (1st-CYS), 41 (CONSERVED-TRP), 89 (hydrophobic), 104 (2nd-CYS) and 118
(J-PHE/J-TRP) must **all** be satisfied. If one is not, the alignment found something
that is not a V-DOMAIN in the register it thinks it is, and every downstream number -
CDR lengths, liability positions, scheme conversion - would be confidently wrong.
So `AbDomain::numbered` stays false, **`residues` stays empty**, and every failure is
named in `anchorFailures`. Wrong numbering is worse than none.

A T-cell receptor is refused by identity, before the anchor test, because a TCR *is* a
V-set domain and would pass all five anchors: if the closest germline set is `TRBV` or
`TRAV`, numbering stops and says so. A sequence whose best germline scores below the
40-bit floor is refused the same way - that floor is also what catches a TCR alpha
chain, since `TRAV`/`TRAJ` could not be retrieved when the pack was built.

### Scheme views

IMGT is canonical internally. Kabat and Chothia are **table-driven views**, produced by
`convertScheme()` from `assets/packs/biologics/scheme-maps.json`, transcribed cell by
cell from the IMGT Scientific chart's numbering-correspondence table. The landmarks
that make that transcription checkable are asserted in `tests/test_biologics_numbering.cpp`:
IMGT 23/41/104 are Kabat H22/H36/H92 and L23/L35/L88.

What is **not** shipped, and is therefore refused rather than guessed:

| View | Status |
|------|--------|
| Kabat | IGHV, IGKV, IGLV - shipped |
| Chothia | IGKV, IGLV - shipped as overrides on Kabat for IMGT 1-104 |
| Chothia, heavy chain | **refused.** The IMGT chart publishes no Chothia (Chothia & Lesk 1987) row for IGHV; only Kabat and Al-Lazikani/AbM are tabulated, and the per-CDR VH page was unreachable |
| Martin (Abhinandan & Martin 2008) | **refused.** Not in the IMGT chart, and no other transcribable source was obtained |
| AHo (Honegger & Pluckthun 2001) | **refused.** Same reason |

A refused conversion returns `numbered == false` with a warning naming the missing
table. A guessed correspondence would renumber every residue plausibly and wrongly,
which is the same failure mode as a bad anchor.

The chart stops at IMGT 115 and tabulates `-` across the heavy CDR3, because Kabat
numbers a CDR3 by its **observed length**, not by a fixed position. CDR3 and FR4 are
therefore produced by the length rule in `cdr3Rule` (heavy H3 95-102 with insertions at
100, light L3 89-97 with insertions at 95, FR4 continuing at H103 / L98), which
reproduces the light-chain cells the chart *does* tabulate (IMGT 112 -> L95A ...
IMGT 115 -> L95D) exactly. A Chothia view says in its warnings that its CDR3 and FR4
came from the Kabat rule.

## 2. Liabilities and developability

`assets/packs/biologics/liabilities.json` holds one cited rule per motif: N-X-S/T
glycosylation, NG/NS/NT deamidation **in that published risk order**, DG/DS/DT
isomerization, DP fragmentation, exposed Met and Trp oxidation, unpaired cysteine,
N-terminal pyroGlu, C-terminal Lys clipping, glycation and RGD. A rule with no citation
is dropped at load with a named error.

A flag is a motif match, nothing more. The **rate** of any of these reactions depends on
local structure, formulation pH, excipients, light and temperature, none of which a
sequence contains. `risk` is a rank index inside a motif family, never a rate.

Exposure weighting is applied **only** when a structure is supplied, through Phase 5's
Shrake-Rupley SASA and the Tien et al. 2013 maxima. Without coordinates,
`exposureKnown` is false and `relativeSasa` is zero - it does **not** default to
"exposed", which would turn every buried methionine into a finding.

Descriptors, all exact arithmetic on the sequence and therefore `Measured`:

| Descriptor | Definition | Source |
|-----------|------------|--------|
| pI, net-charge curve | Henderson-Hasselbalch over the Bjellqvist pKa set, termini included | Bjellqvist et al. 1993 |
| eps280 | `5500*nTrp + 1490*nTyr + 125*nCystine`, cystines counted as pairs | Pace et al. 1995 |
| GRAVY | mean Kyte-Doolittle hydropathy | Kyte & Doolittle 1982 |
| Aliphatic index | `100*(A% + 2.9*V% + 3.9*(I%+L%))` | Ikai 1980 |
| Instability index | `10/L * sum DIWV` over consecutive dipeptides, per chain | Guruprasad et al. 1990 |
| Fv charge symmetry | `q(VH) * q(VL)` at pH 5.5 | Sharma et al. 2014 |

The instability index is a sequence statistic correlated with in-vivo half life in its
1990 training set. It is not a shelf life and nothing here converts it into one.

### TAP-style metrics

`tapPsh`, `tapPpc`, `tapPnc` and `tapSfvcsp` are `NotComputed` unless a structure **and**
its stated origin are supplied. That is not fussiness: the published TAP thresholds were
derived on homology models, and applying them as a verdict to a crystal structure - or
the reverse - is the exact error the refusal prevents. When an origin is given, the four
are computed under an in-house protocol that is printed in full in every `source` string
(CDR vicinity = any residue with a heavy atom within 4 A of a CDR heavy atom;
hydrophobicity = Kyte-Doolittle; weights = relative side-chain SASA) and the report says
outright that the published TAP thresholds are **not** applicable to them. No verdict is
emitted either way.

## 3. Mass and peptide mapping

Every mass comes from `chem::parseFormula` / `monoisotopicMass` / `averageMass` over the
NIST SRD 144 pack. **No mass is a literal.** The data are residue and monosaccharide
*compositions*; the numbers are computed from them, so a corrected isotope table
corrects every mass in the panel.

The same rule covers the deltas that are usually written as constants:

- a disulfide is **two hydrogens removed from the formula**, which is
  2 x 1.007825 Da per bond from the table (the familiar -2.015650);
- pyroglutamate is a cyclisation: Gln loses NH3, Glu loses H2O, applied as composition
  changes;
- deamidation is `O - NH` = 0.984016 Da and the 13C spacing is `[13C] - C` =
  1.003355 Da, both evaluated against the table.

The ladder reports intact, reduced, deglycosylated, per-chain, Fab-like, glycoform
(G0F/G1F/G2F/G0/G2F+NeuAc), pyroGlu and C-terminal-Lys-clipped species.

**Above ~10 kDa the average mass is the reported one** and the monoisotopic column is
secondary, because the all-light-isotope peak of a 150 kDa protein carries a negligible
fraction of the envelope and is not resolved on a routine instrument. Every entry above
10 kDa carries that note.

`MassLadder::requiredResolvingPower` is `R = m/dm` with
`dm = |(13C-12C) - (O-NH)| = 0.019339 Da`: the resolving power needed to tell a
deamidation from an isotope peak at that mass. At 150 kDa that is
`150000 / 0.019339 = 7.76e6`. Below it, **a +1 Da shift is not evidence of
deamidation**, and printing the number is how a mass difference stops being mistaken for
an identification.

Peptide mapping supports trypsin (not before Pro), Lys-C, Glu-C, Asp-N and chymotrypsin
(not before Pro), with missed cleavages, singly-charged b/y ions and variable
modifications expressed as formula arithmetic. Glu-C is applied at Glu only and says so:
in phosphate buffer it also cleaves after Asp, and that variant is not assumed.

## 4. Interfaces and the geometric alanine scan

`bio::interfaceOf()` measures, over the coordinates it is given:

- `BSA = SASA(A) + SASA(B) - SASA(AB)`, the **total** area buried on both sides, so the
  per-side figure that papers usually quote is half of it;
- residue contacts within 4.5 A with their atom-pair count and minimum distance;
- hydrogen bonds (N/O to N/O within 3.5 A), salt bridges (4.0 A between charged groups),
  hydrophobic contacts, pi-pi (5.5 A centroid to centroid), cation-pi (6.0 A) and
  disulfides (2.5 A SG-SG);
- Levy's support / core / rim partition, which needs both the isolated and the complexed
  accessibility, so the report carries both rather than one burial number;
- for a chain that numbers as a V-DOMAIN: the CDR contacts, the paratope and the
  epitope. Framework contacts are **reported**, not filtered out, and their count is
  named in a warning.

Every area carries its SASA parameter string (algorithm, probe radius, point count,
radii set, hydrogen policy). Two tools disagree on an interface area by ~10% from those
choices alone, so an area quoted without them is not reproducible.

### The scan is a rank ordering, not an energy

`bio::alanineScan()` truncates each interface side chain beyond C-beta - keeping
N, CA, C, O, CB - and re-measures the interface, reporting the buried area, atom
contacts, hydrogen bonds and salt bridges that disappear. Gly and Ala are skipped
because there is nothing beyond CB to remove. The remaining structure is **not**
re-minimised, so the lost area is an upper bound on what a real mutant would lose.

`AlanineScanPosition::impact` is therefore a `Provenance::Heuristic` Quantity with an
**empty unit**: `0.5 * lost BSA + 0.3 * lost contacts + 0.2 * lost (H-bonds + salt
bridges)`, each normalised to the largest loss in the scan. There is no solvation term,
no entropy, no relaxation and no electrostatic model in it. `makeQuantity()` throws if
anyone attaches kcal/mol to a Heuristic, and that throw is the enforcement.

`benchmarkSpearman` is `notComputed("a measured benchmark subset")`. A geometric
hotspot proxy is worth exactly its measured rank correlation, and this build ships no
measured antibody/antigen mutation set. Hard-coding a remembered correlation would be
the failure this whole file exists to prevent; when a real set is shipped, the
correlation gets computed and this field stops being NotComputed.

## Data and licences

| Pack | Content | Licence |
|------|---------|---------|
| `imgt-reference.json` | 550 functional IMGT-gapped V alleles (IGHV/IGKV/IGLV/TRBV) and 44 J genes, the region bounds, the five anchors, the CDR layout and the camelid FR2 hallmark | **CC BY 4.0** (IMGT). Attribution, retrieval date and the licence sentence are in the pack |
| `scheme-maps.json` | IMGT -> Kabat / Chothia position tables and the CDR3 length rules, plus the explicit `unavailable` list | **CC BY 4.0** (IMGT) |
| `liabilities.json` | 15 cited motif rules and the exposure threshold | authored here; every rule cites its primary paper |
| `descriptors.json` | Kyte-Doolittle, Guruprasad DIWV and Bjellqvist pKa tables | **BSD-3-Clause** (Biopython), scales cited to their primary papers |

See `NOTICE` section 4a for the obligations these carry.

## Verification

`tests/test_biologics_numbering.cpp`, `tests/test_biologics_mass.cpp` and
`tests/test_biologics_interface.cpp`, tagged `[biologics]`:

- trastuzumab VH (RCSB 1N8Z chain B) numbers with all five anchors and its IMGT CDRs
  `GFNIKDTY` / `IYPTNGYT` / `SRWGGDGFYAMDY` (8/8/13);
- trastuzumab VK, D1.3 VK (1VFB) and Kol VL-lambda (2FB4) number and are typed apart;
- the 1MEL camelid VHH's 26-residue CDR3 numbers with 13 insertion codes on 111 and 112
  only, and 105-110 and 113-117 keep their positions;
- the 1TCR beta chain is rejected with `residues` empty and its reason listed;
- IMGT -> Kabat -> IMGT round-trips identically on all five, and Martin / AHo are
  refused;
- a chain mass agrees computed two independent ways (whole-chain formula versus a
  per-residue sum), reduction adds exactly `2 x nSS x m(H)`, G1F - G0F is exactly one
  Hex per site, and the resolving power matches `150000/0.019339`;
- `1BRS` buries 783 A^2 per side (published ~780), `1VFB` 686 and `1DQJ` 888, with the
  barnase Arg59 / barstar Asp39 hotspots at the top of the geometric ranking;
- **every** alanine-scan `impact` has an empty unit and `Provenance::Heuristic`.
