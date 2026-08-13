# Population PK, uncertainty and mechanistic drug interactions (Phase 13)

This is the layer that turns Phase 4's single exposure curve into a band, reads a real
concentration-time series without assuming a compartment, and computes what a perpetrator
drug does to a victim's exposure. It adds no dependency: the maths runs on Eigen, the one
ODE integrator in `src/numeric/Ode.h`, and the one PK engine in `src/pkpd/PkEngine.h`.

**The permanent boundary.** Everything here is an *exposure scenario* or an *exposure
ratio*. There is no dose, no dose adjustment, no dosing interval and no per-patient
prediction in any DTO, any panel or any agent tool, and the three agent tools refuse the
conversion twice over - in the tool description the model reads, and in the handler, which
gates on `assumptions_acknowledged` and welds a disclaimer into the returned JSON.

---

## 1. The reproducible sampler (`src/sim/Random.{h,cpp}`)

A percentile band, a coverage figure or a CI width is only checkable if it can be
reproduced. `std::normal_distribution` and `std::uniform_real_distribution` are **not
algorithmically specified** by the C++ standard, so the same engine and the same seed give
different numbers on libstdc++ and MSVC. Every number in this phase would therefore change
with the toolchain. The algorithms are pinned instead:

| Piece | Algorithm | Verification |
|---|---|---|
| Bit generator | PCG64-DXSM (128-bit LCG with the cheap multiplier `0xda942042e4dd58b5` and the 2019 double-xorshift-multiply output permutation) | Known-answer against NumPy 2.4.6's `PCG64DXSM.random_raw` from a placed raw state |
| Uniform | top 53 bits x 2^-53, so every representable double in `[0,1)` is equally likely and 1.0 is unreachable | range and 10-bin equidistribution over 2e5 draws |
| Normal | Wichura AS241 (PPND16) inverse CDF, one uniform per variate | eight quantiles against `scipy.stats.norm.ppf` to 1e-14 relative |
| Correlated normals | Eigen `LLT` of Omega, `eta = L z` | recovered covariance over 1e5 draws |
| Stratified sampling | Latin hypercube, one point per stratum per dimension | every stratum occupied exactly once |
| Correlation induction | Iman-Conover on van der Waerden scores | target Spearman 0.7 recovered to 0.02, marginals bit-identical |

Two details are load-bearing. Inverse transform is used for normals rather than Box-Muller
or a ziggurat because it consumes **exactly one uniform per variate**, so the stream
position depends only on how many variates were drawn. And `imanConover` inverts the
Spearman target through `rho = 2 sin(pi * r_s / 6)` before driving the Pearson correlation
of the scores - skipping that step is the classic error where you ask for 0.70 and measure
0.68.

`Pcg64Dxsm::fromRawState` exists only so the known-answer test can start where the
reference implementation was placed: seeding schemes differ between libraries, the state
transition and output function do not, and it is those the contract is about.

## 2. The population layer (`src/sim/Population.{h,cpp}`)

`simulatePopulation(spec, regimen, variability)` integrates every subject through
`pkpd::simulate`. There is no second PK model, so a band can never disagree with the
typical-value curve drawn beside it.

Three layers, separately toggled, because they answer different questions:

1. **Between-subject variability.** `eta ~ MVN(0, Omega)` on the log scale, so
   `Pi = theta_i * exp(eta_i)` and no parameter can go negative. Omega is the covariance of
   the log-scale random effects; an SD of 0.30 is roughly a 30% CV on the parameter. A
   non-positive-definite Omega is **reported, not repaired** - a nudged Omega is no longer
   the variability that was entered.
2. **Parameter uncertainty.** Latin hypercube over the unit hypercube, rank-correlated by
   Iman-Conover to the correlation implied by the fit covariance, then mapped to log-normal
   multipliers.
3. **Residual error.** Proportional then additive, clamped at zero, because a negative
   observed concentration is not a measurement.

Concentrations are stored in **one allocation** indexed `[timeIndex * subjects + subject]`.
Time is the slow axis, so every subject at a given time is contiguous and the percentile
pass sorts a contiguous span instead of gathering across 1000 separate vectors. At
1000 subjects x 4000 points that is a single 32 MB block.

Percentiles are the linear-interpolated "type 7" definition (R's `quantile()`, NumPy's
`percentile()`), stated explicitly because nearest-rank and interpolated definitions differ
by a whole order statistic at small N. At most 50 individual trajectories are kept for the
faint-line overlay; the percentiles use every subject, and the panel says which is which.

`provenanceStatement` is a required sentence, not decoration:

> This band is the variability that was entered - Omega, the fit covariance and the
> residual error above - propagated through the stated PK model with seed N; it is not a
> prediction about any individual, and it is not a dose.

## 3. Noncompartmental analysis (`src/sim/Nca.{h,cpp}`)

NCA assumes no structural model, so it has to be strict about the two places it can still
go wrong.

- **AUC** is linear-up/log-down. A linear trapezoid over a decaying exponential
  systematically *overestimates* the area; the log-down form
  `dt*(c1-c2)/ln(c1/c2)` is the exact integral of the segment. The test asserts equality
  with `(C0/ke)(1 - exp(-ke*T))` to 1e-11 at three sampling densities.
- **lambda_z** is fitted only to points **strictly after Tmax**, over at least three of
  them, and the window is chosen by **adjusted** R-squared so an added point has to earn
  its degree of freedom. Fewer than three such points yields `NotComputed`, and with it no
  half-life, no AUCinf, no Vz.
- **Clast is the regression-predicted value**, not the last observed point, so the whole
  extrapolated tail does not inherit one measurement's noise.
- **Above 20% extrapolated**, `extrapolationUnreliable` is set and each derived quantity -
  AUCinf, CL(/F), Vz(/F), Vss, MRT, AUMC - gets its own warning naming it.
- **Vss is IV-only.** After an extravascular dose MRT contains the mean absorption time, so
  `CL * MRT` is not a volume at all; the field reads `NotComputed` with that reason.
- At steady state (`tauH > 0`) the exposure per dose is AUCtau, so clearance is
  `dose/AUCtau`, not `dose/AUCinf`, and Cavg and swing come from the same interval.

## 4. Drug interactions (`src/sim/Ddi.{h,cpp}`)

**One source of physiology.** `assets/packs/physiology.json` carries `Qh = 97 L/h`,
`Qen = 18 L/h`, GFR and per-enzyme `kdeg`, and `core::physiology()` is the only reader.
Phase 1.7's well-stirred model in `src/chem/AdmetModel.h` used to hard-code `Q_H = 90.0`
while this pack said 97 - a constant that appears twice has already diverged. It now reads
`core::physiology().hepaticBloodFlowLPerH`, and a missing pack makes the prediction report
the missing pack rather than fall back to a built-in number.

**Three models, in increasing order of what they claim.**

*FDA basic-model R-values* are screening numbers and need no victim parameters:

```
R1     = 1 + [I]h,u / Ki           R1,gut = 1 + [I]g / Ki
R2     = (kdeg + kobs) / kdeg      kobs   = kinact*[I] / (KI + [I])
R(ind) = 1 / (1 + d*Emax*[I]/(EC50 + [I]))
```

*The mechanistic static AUCR* combines all three mechanisms over both organs:

```
A = kdeg/(kdeg + kobs)   B = 1 + d*Emax*[I]/(EC50+[I])   C = 1/(1 + [I]/Ki)
AUCR = 1/((A*B*C)_h * fm + (1 - fm))  *  1/((A*B*C)_g * (1 - Fg) + Fg)
```

- **No fm, no AUCR.** `aucRatio` is `notComputed("fm")`. Defaulting fm to 1 is exactly what
  turns a 1.3-fold interaction into a 5-fold one.
- **No Fg** leaves `gutIncluded = false` and reports `aucRatioHepaticOnly`.
- `theoreticalCeiling = 1/(1-fm)` is displayed beside the ratio, because it is the number
  that says whether a predicted ratio is even attainable.
- `dominantMechanism` names which factor moved furthest from 1, so the reader knows whether
  to expect an immediate onset (reversible) or one over several enzyme half-lives.

*The dynamic enzyme model* integrates, through `numeric::rk4Integrate`:

```
dE/dt    = kdeg*(1 + d*Emax*I/(EC50+I)) - kdeg*E - kinact*I/(KI+I)*E
CLint(t) = CLint0 * E / (1 + I/Ki)
```

At constant `I` its steady state must **equal** the static model, and `EnzymeTimeCourse`
carries both numbers plus their difference so the agreement is visible instead of asserted
in a comment. The test requires agreement below 1e-9 and, separately, that the A*B*C factor
implied by the reported AUCR is the same number again. The steady state is integrated out
to 40 time constants rather than to the display horizon: a 24 h plot of a 36 h-turnover
enzyme has not converged, and comparing an unconverged value with the static model would
report a disagreement that is the plot's fault.

**Organ impairment** is an editable exposure ratio, never a formula from a clinical score:

```
CL_impaired/CL_normal = fe*renalFunctionRatio + (1 - fe)*hepaticClintRatio
AUC ratio             = 1 / that
```

with `fe` a victim input (absent -> `NotComputed`) and the low-extraction limit of the
well-stirred model stated as an assumption, together with the reminder that a
high-extraction drug's clearance is bounded by Qh and would move far less. There is
deliberately **no Child-Pugh-to-clearance mapping and no creatinine-clearance-to-dose
formula**: those are dosing decisions.

**Whole-body PBPK is absent** for the reason stated at the top of `src/data/Population.h`:
its required fu, blood-to-plasma ratio, tissue partition, Papp and transporter inputs are
not derivable from anything BioCAD has, so a PBPK profile here would be precise fiction.

## 5. Panels and tools

| Surface | id | What it shows |
|---|---|---|
| Panel | `PopPk` | percentile bands with <= 50 faint trajectories under them, occupancy of the median profile, an Omega SD / correlation editor, the seed, the NCA table for the median profile, the parameter table coloured by provenance, and the assumption list |
| Panel | `InteractionScenarios` | R-values, the static AUCR with its ceiling and dominant mechanism, the dynamic enzyme-activity plot with the static comparison printed underneath, and the impairment sliders. No risk score and no red/green verdict anywhere |
| Tool | `simulate_population` | requires `assumptions_acknowledged`; returns a thinned band, the seed, Omega, and a disclaimer refusing dose conversion |
| Tool | `noncompartmental_analysis` | analyses supplied observations; the disclaimer refuses turning CL or Vz into a dose |
| Tool | `predict_interaction` | requires `assumptions_acknowledged`; returns the R-values, the AUCR (or `NotComputed` for a missing fm), the dynamic/static agreement and the kdeg used, with a disclaimer refusing conversion into a dose, a contraindication or a risk category |

The population panel simulates **on demand**, not every frame: hundreds of subjects over a
48 h horizon is tens of millions of RK4 steps, and re-running it at 60 Hz would both make
the UI unusable and make the displayed seed meaningless, because the band would be replaced
before anyone could read it.

## 6. Verification

`tests/test_sim_random.cpp`, `tests/test_sim_nca.cpp` and `tests/test_sim_ddi.cpp` - 21
cases. The load-bearing assertions:

- PCG64-DXSM reproduces NumPy's `PCG64DXSM.random_raw` exactly from a placed state.
- AS241 matches `scipy.stats.norm.ppf` to 1e-14 relative at eight quantiles.
- Iman-Conover reaches a target Spearman of 0.7 within 0.02 while leaving every marginal
  bit-identical.
- IV bolus D = 100 mg, CL = 5 L/h, V = 50 L gives AUCinf = 20 mg*h/L (1e-9),
  Vss = 50 L (1e-7), MRT = 10 h (1e-7).
- The log-down trapezoid is exact on monoexponential data to 1e-11.
- Truncating that profile at one half-life gives exactly 50% extrapolated and flags every
  derived quantity by name.
- The same seed gives bit-identical percentile bands; a different seed does not.
- The dynamic and static DDI steady states agree to better than 1e-9.
- A missing fm returns `NotComputed("fm")`; a missing Fg reports hepatic-only.
