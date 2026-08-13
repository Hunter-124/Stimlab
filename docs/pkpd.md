# Pharmacodynamics and PK/PD

The PK/PD subsystem answers one question honestly: **given parameters you supply, what
exposure and what target occupancy follow?** It does not answer "what dose should I take",
and there is no code path that could be made to.

> The panel carries a fixed banner saying so, the agent tool `simulate_exposure` refuses to
> run without an explicit `assumptions_acknowledged` flag, and the assistant's system prompt
> forbids recommending a dose, a dose change, or a personal regimen for anyone.

## Where the code lives

| Path | Contents |
| --- | --- |
| `src/numeric/Optimize.{h,cpp}` | Levenberg-Marquardt least squares with an analytic Jacobian |
| `src/numeric/Ode.{h,cpp}` | Fixed-step RK4 and trapezoidal integration |
| `src/pkpd/Fits.{h,cpp}` | 4PL dose-response, Cheng-Prusoff, Schild |
| `src/pkpd/PkEngine.{h,cpp}` | The four PK structural models, closed forms, occupancy |
| `src/pkpd/Pharmacodynamics.{h,cpp}` | `RealPharmacodynamics`, implementing `IPharmacodynamicsModule` |
| `src/data/Domain.h` | Every DTO named below |
| `tests/test_numeric.cpp`, `tests/test_pkpd_fits.cpp`, `tests/test_pkpd_engine.cpp` | The oracles |

There is exactly one fitter and one integrator in the tree. A second curve fit with slightly
different convergence behaviour is a bug factory, so anything that needs least squares goes
through `numeric::levenbergMarquardt`.

## Dose-response

Fitted in log-concentration form, which is the parameterisation whose Jacobian is well
conditioned across the asymptotes:

```text
E(x) = Bottom + (Top - Bottom) / (1 + 10^(nH * (log10 EC50 - x))),   x = log10([A])
```

With `u = ([A]/EC50)^nH` and `delta = Top - Bottom`, the analytic Jacobian is

```text
dE/dTop       = 1/(1+u)
dE/dBottom    = u/(1+u)
dE/dlogEC50   = delta * nH * ln10 * u / (1+u)^2
dE/dnH        = -delta * u * ln(u) / (1+u)^2
```

Initial guesses are `Top = max(y)`, `Bottom = min(y)`, `log10 EC50` at the x nearest the
response midpoint, and `nH = 1`. Optional `1/y^2` weighting is available. Parameter standard
errors come from the covariance matrix and land in each `Quantity::error`; when the system is
rank-deficient the errors are **empty**, not fabricated.

`nH` is labelled the **empirical slope**, never "cooperativity". In a functional assay,
signal amplification and receptor reserve bend the slope, so a value above 1 is not evidence
of cooperative binding, and calling it that would launder an assay artefact into a mechanism.

## Cheng-Prusoff

`ChengPrusoffInput` carries a modality and the fields that modality requires:

| Modality | Ki | Required |
| --- | --- | --- |
| `Competitive` | `IC50 / (1 + [S]/Km)` | `substrate`, `km` |
| `Uncompetitive` | `IC50 / (1 + Km/[S])` | `substrate`, `km` |
| `Noncompetitive` | `IC50` | - |
| `RadioligandBinding` | `IC50 / (1 + [L]/Kd)` | `radioligand`, `kdRadioligand` |

**There is no default modality and no silent fallback.** A missing required field returns a
`NotComputed` `Quantity` that names the field. This matters more than it looks: the competitive
and uncompetitive corrections diverge with `[S]/Km` - 10x apart at `[S] = 10*Km` and 100x apart
at `[S] = 100*Km` - and ChEMBL has no `UNCOMPETITIVE` action type to disambiguate an imported
record with. Guessing competitive would silently invent orders of magnitude.

When an enzyme concentration is supplied and Ki approaches it, the classic and the
depletion-corrected values are reported side by side rather than one replacing the other:
the classic Cheng-Prusoff relation assumes `[I] >> [E]t`, and that assumption is stated.

## Schild

`log10(DR - 1)` regressed on `log10([B])`; the report carries `pA2`, the slope, and the
slope's 95% confidence interval.

**If the CI excludes 1, `KB` is `NotComputed` and `kbUsable` is false.** A KB read off a
non-unit-slope Schild plot is meaningless - the antagonism is not simple competitive - and
printing one anyway is the single most common way this analysis is abused.

## PK engine

Four structural models (`PkModel`): IV bolus, IV infusion, one-compartment with first-order
absorption, and two-compartment with first-order absorption. Optional Michaelis-Menten
elimination. Arbitrary dose-event lists for multiple dosing.

The models are **integrated** with RK4 rather than evaluated only in closed form, for two
reasons that are not stylistic:

1. The closed forms are numerically singular at `ka == ke`. The Bateman function's `ka - ke`
   denominator vanishes; the true limit `C = F*D*ke*t*exp(-ke*t)/V` is finite, and the
   implementation returns it, but a naive closed-form-only engine produces infinity.
2. Nonlinear elimination has no closed form at all.

The closed forms still exist, as the **unit-test oracle** and the summary readout:
`batemanConcentration`, `accumulationRatio` (`Rac = 1/(1 - exp(-ke*tau))`) and
`steadyStateAverage` (`Cav,ss = F*D/(CL*tau)`). The integrated result is asserted against them
to 1e-6 in the linear case; a divergence means the integrator is wrong, and the test says so.

Flip-flop kinetics (`ka < ke`) are detected and flagged: the terminal phase then reflects
absorption, not elimination, and reading a "half-life" off it is a mistake.

Michaelis-Menten runs return `halfLife` as `NotComputed`, because a saturable system does not
have one.

## Occupancy is the headline output

```text
theta(t) = [A]_free(t) / (Kd + [A]_free(t))
```

This is the honest synthesis of PK and PD: it needs a free concentration and a Kd, and
nothing else - no Emax, no transducer ratio, no tissue assumption. A sigmoid-Emax plus
effect-compartment link is available but is explicitly a second-order, assumption-heavy view.

If the supplied Kd is `NotComputed` or non-positive, the curve is empty and names the missing
Kd. There is no default Kd.

## Every parameter is tagged

`PkModelSpec` holds a `Quantity` per parameter, so each one carries its own provenance:
`Measured` (entered from a paper), `Predicted` (a model ran), or `Heuristic`/assumed (a stated
default). The panel prints one assumption line per non-measured parameter, verbatim, under
the plot.

**F, ka and fu have no credible structure-only predictor.** They default to assumed and the
panel says exactly that. Any predicted clearance additionally displays the PKSmart
generalisation figure: only 10.4% of high-clearance compounds are predicted within 2-fold
(Seal et al. 2025, *J Cheminform* 17:147). A number that is right one time in ten should look
like one.

## Agent tools

| Tool | Behaviour |
| --- | --- |
| `fit_dose_response` | Fits a 4PL to supplied points and returns the parameters with their errors |
| `convert_ic50_to_ki` | Cheng-Prusoff, refusing with a named missing field rather than assuming a modality |
| `simulate_exposure` | Requires `assumptions_acknowledged: true` in its arguments, and errors without it |

The acknowledgement flag exists so the model cannot slide from "here is a curve" into a
dosing narrative without the caller having explicitly accepted that the whole thing rests on
assumptions.

## Acceptance oracles

These are the numbers the test suite pins, and they are all checkable by hand:

- 4PL recovers `Top=100, Bottom=0, EC50=1e-7, nH=1.2` from noiseless synthetic data to 1e-9.
- Oral one-compartment (`D=100 mg, F=0.8, ka=1.2/h, ke=0.15/h, V=50 L`) matches the Bateman
  closed form to 1e-6 at t = 0.5, 1, 2, 6 and 24 h.
- The `ka == ke` limit is finite and matches the analytic limit.
- IV bolus `D=100 mg, CL=5 L/h, V=50 L` gives `AUC(inf) = 20 mg*h/L` and `t1/2 = 6.93 h`.
- `Rac` over 10 q12h doses matches `1/(1 - exp(-ke*tau))`.
- Cheng-Prusoff returns `NotComputed` naming `Km` when `Km` is absent.
- A Schild fit whose slope CI excludes 1 marks `KB` unusable.
