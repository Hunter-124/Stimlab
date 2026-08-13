# Assay and biophysics workbench

Two panels, `Assay` ("Assay Workbench") and `AssayDesign` ("Assay Design"), and one contract,
`IAssayModule` (`src/contracts/IModules.h`, implemented by `src/modules/AssayModule.*` over
`src/assay/{Dataset,Qc,Fits,Biophysics,Design}.*`).

This is the only part of BioCAD whose input is data somebody measured. Everything else turns a
structure into a prediction; this turns a plate reader export into numbers with error bars. That
makes the provenance rule sharper here than anywhere else: **a raw well value is `Measured`, a
fitted parameter is `Model`, and the two are never rendered in the same colour.**

Scope: this is analysis and experimental design. A concentration ladder for a plate is not a
dose for a person, and a dilution volume is not a synthesis procedure.

## Import format

The canonical internal representation is long, one row per well. Columns are matched
case-insensitively and in any order; unknown columns survive as dataset metadata rather than
being dropped, because the instrument knows things about the run that BioCAD does not.

| Column | Meaning |
| --- | --- |
| `plate_id` | Groups rows into plates. |
| `well` | `A1` .. `AF48`; the letter block gives the row, the number the column. |
| `role` | `sample`, `positive`, `negative`, `blank`, `reference`, `empty`, `unknown`. |
| `sample_id` | What is in the well. |
| `series_id` | Groups the wells of one concentration series. For an `[S] x [I]` matrix this is the numeric inhibitor concentration. |
| `concentration`, `conc_unit` | The nominal concentration. For a sensorgram this is the analyte concentration. |
| `replicate` | Replicate index. |
| `readout`, `readout_unit` | The measurement. |
| `time_s` | Sensorgram / ITC time axis. |
| `temperature_c` | Melt-curve axis. |
| `excluded`, `exclusion_rule` | An exclusion and the rule id that made it. |

A grid sniffer recognises 96/384/1536 reader exports by a `1..12/24/48` header row with `A..AF`
row labels and melts each block into the same representation. `AssayDataset::detectedLayout`
reports which path ran.

**An exclusion never deletes a point.** `Well::excluded` hollows the marker and
`Well::exclusion_rule` records which rule hollowed it, so the decision is auditable after the
fact. Outlier rules (Grubbs, Dixon Q for n = 3..10, Tukey fences) are opt-in.

## Plate QC and its published bands

`IAssayModule::qc` returns a `QcReport` of `Quantity` fields, each of which is either right or
`NotComputed` naming the prerequisite it lacked.

| Statistic | Formula | Bands |
| --- | --- | --- |
| Z-prime | `1 - 3(sp + sn) / abs(mup - mun)` | `>= 0.5` excellent, `0` to `0.5` marginal, `<= 0` unusable |
| Robust Z-prime | the same with median and `1.4826 * MAD` | same bands |
| SSMD | `(mup - mun) / sqrt(sp^2 + sn^2)` | reported, not banded |
| Signal / background | `mup / mun` | `NotComputed` when the background mean is zero |
| Signal / noise | `(mup - mun) / sn` | `NotComputed` without negative-control variance |
| `%CV` | `100 * sd / mean` | positive ratio-scale data only |

Two behaviours are deliberate and not negotiable:

- **Z-prime is `NotComputed` without BOTH a positive and a negative control.** Computing one
  from whatever the extreme wells happened to be produces a number that looks like a QC verdict
  and is not one.
- **Edge, row and column effects are reported, never corrected.** The outer ring is compared to
  the interior with a Mann-Whitney test and rows/columns with Kruskal-Wallis, and the p-values
  are printed. Median-polishing a gradient away hides the pipetting problem that caused it.
  B-score normalisation exists as an explicit user choice, not as an automatic repair.

A plate with `Z' <= 0` reads, verbatim from `QcReport::interpretation`:

```
Z-prime <= 0: unusable, the control distributions overlap (bands: >= 0.5 excellent,
0 to 0.5 marginal, <= 0 unusable). Edge, row and column effects are reported as
p-values and never auto-corrected.
```

## Which models exist, and which engine owns them

`IAssayModule::fit` dispatches by model, because a concentration series and a trace are not the
same data:

| Model | Fitted from | Engine |
| --- | --- | --- |
| 4PL, 5PL | a concentration series | `assay::fitSeries` (`Fits.cpp`) |
| Michaelis-Menten, Hill, substrate inhibition, Morrison tight binding | a concentration series | `Fits.cpp` |
| competitive / uncompetitive / noncompetitive / mixed inhibition | the full `[S] x [I]` matrix | `Fits.cpp` |
| Langmuir 1:1, two-compartment mass transport | sensorgrams (response vs time, per analyte concentration) | `Biophysics.cpp` |
| Boltzmann melt, two-state thermodynamic melt | signal vs temperature | `Biophysics.cpp` |
| Wiseman one-set-of-sites isotherm | an ITC experiment | `Biophysics.cpp` |

Notes that change what you may quote:

- **The `[S] x [I]` global fit is BioCAD's ONE producer of `InhibitionModality`.** It ranks the
  modalities by AICc and answers `Unknown` when the top two differ by under 2 units. Phase 4's
  Cheng-Prusoff consumes that answer; a Ki computed with the wrong modality is off by a factor
  of `[S]/Km`, which is 10x at `[S] = 10*Km`. A Lineweaver-Burk plot is available as a
  diagnostic only - no fitter in BioCAD attaches to a transformed axis, because the transform
  distorts the error structure it would be fitting.
- **An ITC fit requested from a bare well list is refused by name.** `n`, `K` and `dH` are only
  meaningful against a known cell volume, macromolecule and titrant concentration and a blank
  titration, none of which is a property of a well, so the module returns `converged = false`
  with the missing metadata named rather than inventing a cell volume.
- **A sensorgram's injection stop is inferred** from the response maximum when the export does
  not carry it, and the inference is returned as a fit warning, because it biases `kd`.
- **Steady-state KD** appears as its own parameter and is `NotComputed` when the association
  phase never reached equilibrium.

### The 5PL rule: EC50 is not C

For the five-parameter logistic

```
y = D + (A - D) / [1 + (x/C)^B]^G
```

`C` is the inflection parameter, **not** the half-maximal concentration. The midpoint is

```
EC50 = C * (2^(1/G) - 1)^(1/B)
```

which equals `C` only when `G = 1`, i.e. when the curve is a 4PL. `FitResult::derivedEc50`
always carries the corrected value and the panel labels it "EC50 / midpoint"; the raw `C`
appears in the parameter table under its own name `log10C` so the two can never be confused.
Reporting `C` as an EC50 is the most common 5PL error in the literature and BioCAD has no code
path that can do it. A 5PL fit with fewer than eight concentrations is refused: the asymmetry
parameter is not identifiable from a short ladder.

An EC50 outside the tested concentration range sets `FitResult::extrapolated`, and the panel
renders it **grey with the reason** instead of as a result. It is a bound, not a potency.

`nH` is labelled "empirical slope", never "cooperativity": in a functional assay, amplification
and receptor reserve bend the slope.

## Design simulation (10.4)

`IAssayModule::simulate` takes an `AssayDesignSpec` - a truth model, its parameters, the ladder,
the replicate count, and the error structure you believe you have - and answers what the design
would actually recover.

The simulated plate carries, in this order:

1. the truth value at the ACTUAL concentration in the well,
2. a linear plate gradient (`plateGradientPct`, multiplicative in row and column),
3. DMSO tolerance signal loss scaled to the top concentration,
4. proportional readout noise (`proportionalNoiseCv`),
5. additive detector noise (`additiveNoiseSd`).

**Pipetting error compounds.** `pipettingCv` is a per-transfer lognormal factor and the
accumulated factor multiplies down the ladder, so the bottom of a 10-step ladder carries
`sqrt(10)` times the CV of one transfer. Drawing an independent error per well would understate
the uncertainty at exactly the concentrations that set the EC50. The REPORTED concentration is
the nominal one, because that is what the experimenter writes on the plate map - the pipetting
error is invisible to the analysis, which is the whole reason it hurts.

Every simulated plate is serialised to the long CSV above and pushed through
`assay::importText` -> `assay::plateQc` -> `assay::fitSeries`. A design simulator with a private
analysis path measures the private path, not the product.

### Reproducibility

The RNG is **PCG64-DXSM** (128-bit LCG state, the 2019 DXSM output permutation) with a
Box-Muller normal transform, both implemented in `src/assay/Design.cpp`. `std::mt19937_64` is
specified bit-exactly but `std::normal_distribution` and `std::lognormal_distribution` are not -
the standard fixes their distributions, not their algorithms - so a report built on `<random>`
distributions would change with the standard library. Each run's seed is derived arithmetically
from `spec.seed`, so run 700 is reproducible without replaying runs 0-699, and the same spec
yields a byte-identical `recoveredEc50` vector.

### What it reports, and the numbers measured here

Over `replicateRuns` seeded repetitions the report carries median Z-prime, median recovered
EC50, the median log10 CI width, the convergence rate, and - the headline - **empirical CI
coverage**: the fraction of runs whose OWN reported 95% interval contained the truth. A fitter
whose nominal 95% interval covers the truth 70% of the time is not reporting a 95% interval, and
counting is the only way to find that out.

Measured on this tree (`src/assay/Design.cpp` + the landed `Dataset/Qc/Fits`), truth
`A = 100, B = 1, C = 1e-7 mol/L, D = 0`, a ten-point half-log ladder from 1e-5 mol/L, three
replicates, additive noise SD 2.0 readout units, seed 20260813, 1000 runs:

```
coverage: 94.4000% (model) convergence 100.00% medianZ'=0.8840 medianEC50=1.002499e-07
CIwidth(log10)=0.07384
```

94.4% against a nominal 95% over 1000 runs is inside the Monte Carlo error of the nominal level
(the standard error of a coverage estimate at p = 0.95, n = 1000, is 0.69 points), so the Wald
interval this design reports means what it says. The report warns whenever coverage falls
outside 90-99%, and the fix is a wider ladder or more replicates - not reporting the interval
anyway.

### D-optimal concentration placement

`assay::dOptimalLadder` runs Fedorov coordinate exchange on the log-determinant of the
information matrix `F'F`, with sensitivities taken by central differences so the criterion works
for every truth model the simulator supports. The candidate set is **not a continuous interval**:
a bench ladder can only reach `top / sqrt(fold)^j`, i.e. the entered ladder plus the one
intermediate point a second half-step dilution can make. Optimising over a continuum would
return concentrations nobody can pipette.

Measured for the spec above (300 runs per arm):

```
D-optimal gain over the entered ladder: 1.17110
CI width log10: entered uniform ladder 0.07423   D-optimal ladder 0.06653
```

The D-optimal design piles points near the asymptotes and the midpoint, which is what a
D-criterion does; it is optimal for the parameters, not for seeing the curve, so the panel prints
both ladders rather than replacing yours.

### Dilution and mass calculator

`assay::massForStock(molarMass, molarity, volumeL)` and `assay::serialDilution(...)` are exact
arithmetic, not models. The ladder warns when a transfer volume falls below the stated floor
(2 uL by default), because that is where the pipetting CV you assumed stops applying, and it
reports the `sqrt(n)` compounding factor at the last step.

## Agent tools

| Tool | What it does |
| --- | --- |
| `import_assay` | Parses an export and reports the layout, plates, roles and every warning. |
| `assay_qc` | QC for one imported plate, with the Z-prime bands stated in the result. |
| `fit_dose_response` | 4PL. Takes either a `points` array (Phase 4 path) or a `plate_csv` plus `series_id`, which routes through the assay import/fit path and therefore also reports the extrapolation flag. |
| `fit_binding_kinetics` | Global 1:1 Langmuir or mass-transport fit over a sensorgram set. |
| `fit_enzyme_inhibition` | Global `[S] x [I]` fit, ranking modality by AICc; `Unknown` under 2 AICc units. |
| `design_assay` | The forward simulation, returning the coverage-first summary. |

`design_assay` returns the summary quantities and drops the per-run vector: thousands of numbers
say nothing a model can use, and the coverage figure is the result.
