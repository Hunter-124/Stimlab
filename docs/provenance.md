# Provenance: every number carries its origin

A number without its origin is a number that lies. `-9.3 kcal/mol` from AutoDock Vina and
`-9.3` from a descriptor regression render identically, are stored identically, and are
copied into a slide deck identically - but one is a scoring-function evaluation of a
constructed pose and the other is a rank-ordering heuristic with no physical meaning at
all. BioCAD refuses to let those two values look the same anywhere: not in the type
system, not in the JSON, and not on screen.

The mechanism is a five-tier `Provenance` enum plus a `Quantity` type that carries it,
both in `src/data/Domain.h`, with the one dangerous combination made unrepresentable by a
throwing factory in `src/data/Domain.cpp`.

## The five tiers

Defined at `src/data/Domain.h:39-45`; labels at `src/data/Domain.cpp:18-27`; colours at
`src/ui/Theme.cpp:137-146`.

|Tier|Exact meaning|Units|Colour (RGB)|
|---|---|---|---|
|`Measured`|Exact geometry or statistics computed from the input, or an experimental value retrieved with a citation|Physical|Green `52, 211, 153`|
|`Predicted`|A published model actually ran. Physical units, and a benchmark error is mandatory|Physical|Blue `96, 165, 250`|
|`Model`|A constructed artefact - a built structure, a docked pose - with no energy claim attached|Physical, but describing an artefact|Purple `192, 132, 252`|
|`Heuristic`|Rank ordering only. Arbitrary scale. Physical units are forbidden|**Must be empty**|Amber `251, 191, 36`|
|`NotComputed`|A prerequisite was missing; `source` names what|None|Grey `138, 149, 167`|

The provenance palette is deliberately distinct from the `Verdict` palette
(`src/ui/Theme.cpp:124-132`): provenance says how a number was obtained, never whether the
news is good. A green `Measured` value can be a terrible result.

The enum serializes to stable lowercase strings - `"measured"`, `"predicted"`, `"model"`,
`"heuristic"`, `"not computed"` - via `NLOHMANN_JSON_SERIALIZE_ENUM` at
`src/data/Domain.h:47-53`, so a persisted `Quantity` keeps its tier across runs and across
the run store. `tests/test_provenance.cpp:50-58` pins that round trip.

## `Quantity`

`src/data/Domain.h:67-74`.

```cpp
struct Quantity {
    double      value = 0.0;
    std::string unit;          // empty string is REQUIRED when provenance == Heuristic
    double      error = 0.0;   // 0 means "no error bar available"
    Provenance  provenance = Provenance::NotComputed;
    std::string source;        // citation, engine name, or model+benchmark
};
```

|Field|Contract|
|---|---|
|`value`|The number itself. Nothing else in the struct is optional decoration around it.|
|`unit`|Free-form physical unit string (`"kcal/mol"`, `"nM"`, `"L/h"`). Empty is required for `Heuristic` and is also how a genuinely dimensionless ratio is expressed.|
|`error`|The error bar, in the same unit as `value`. Zero means "no error bar available", not "exact". For a `Predicted` value the originating model's benchmark error belongs here.|
|`provenance`|Defaults to `NotComputed`, so a default-constructed `Quantity` claims nothing.|
|`source`|Who says so: a citation, an engine version string (`"AutoDock Vina 1.2.5"`), a model plus its benchmark (`"ThermoMPNN (Fireprot-HF RMSE 1.55)"`), or - for `NotComputed` - the name of the missing prerequisite.|

Two helpers accompany it, both declared at `src/data/Domain.h:79-83`:

```cpp
Quantity makeQuantity(double value, std::string unit, double error, Provenance p,
                      std::string source);
Quantity notComputed(std::string missingPrerequisite);
```

`notComputed("Km")` produces a tier-`NotComputed` quantity whose `source` is the missing
field name, which is what the UI prints (`src/data/Domain.cpp:45-50`,
`tests/test_provenance.cpp:36-41`).

## The invariant: a heuristic cannot carry a unit

This is enforced, not documented. `makeQuantity` throws `std::invalid_argument` when a
`Heuristic` arrives with a non-empty unit (`src/data/Domain.cpp:29-43`):

```cpp
if (p == Provenance::Heuristic && !unit.empty()) {
    throw std::invalid_argument(
        "heuristic quantity must not carry a physical unit (got \"" + unit +
        "\" from " + source + ")");
}
```

The aggregate is still constructible directly - this is plain data with no private state -
but every construction site in the tree goes through `makeQuantity`, and
`tests/test_provenance.cpp:21-34` locks the behaviour in. The point is that
"kcal/mol on a rank-ordering score" is not a style violation to be caught in review; it is
an exception at the point of construction.

## `weakest()` - and the ordering trap

`src/data/Domain.h:60-62`:

```cpp
constexpr Provenance weakest(Provenance a, Provenance b) {
    return static_cast<int>(a) >= static_cast<int>(b) ? a : b;
}
```

A derived quantity inherits the least trustworthy of its inputs. A ligand efficiency
computed from a `Measured` dG and a `Heuristic` heavy-atom estimate is `Heuristic`; there
is no averaging and no partial credit.

The trap: the enumerators are declared **strongest-first**
(`Measured, Predicted, Model, Heuristic, NotComputed`), so the weakest tier is the
**maximum** integer value, not the minimum. `std::min` over provenance values silently
produces the *most* trustworthy input and is exactly backwards. The header says so at
`src/data/Domain.h:35-37`; always call `weakest()` rather than open-coding the comparison.
`tests/test_provenance.cpp:43-48` covers the four interesting pairs, including
`weakest(Heuristic, NotComputed) == NotComputed`.

## How it renders

`drawQuantity` (`src/ui/Panels.cpp:43-63`, declared `src/ui/Panels.h:16`) is the single
widget every panel uses for a derived number:

```
Best affinity  -9.30 +/- 0.42 kcal/mol   (model - AutoDock Vina 1.2.5)
```

The whole string is drawn in `theme::provenanceColor(q.provenance)`. Three deliberate
choices:

- The error bar is **part of the value**, appended as `+/- error` when `error > 0`, not
  hidden in a tooltip a screenshot will never show.
- The tier word and the `source` are appended in parentheses, inline. A number can be read
  out of a screenshot together with its justification.
- `NotComputed` short-circuits to `not computed - needs <source>`, so a missing
  prerequisite is a named gap rather than a zero.

The Presets panel reuses the same palette for a non-numeric fact: a pack target is drawn
`dockable` in `Model` purple or `no box` in `NotComputed` grey (`src/ui/Panels.cpp:1549-1551`).

## Worked example: why `kdFromDeltaG` exists and is never wired to a dock

AutoDock Vina's reported standard error is 2.85 kcal/mol (Trott & Olson 2010, PMC3041641).
At 298.15 K, `RT = 1.987204259e-3 * 298.15 = 0.5925 kcal/mol`, so a 2.85 kcal/mol error is
a factor of `exp(2.85 / 0.5925) = ~123` in Kd. A docked score presented as a nanomolar
affinity is therefore a number that could be anywhere from 8 nM to 120 uM. Converting it
manufactures precision that does not exist.

So the conversions exist, exactly and reversibly, and take a **measured** dG or pActivity:

```cpp
inline constexpr double kGasConstantKcal = 1.987204259e-3;   // kcal/(mol.K)
inline double kdFromDeltaG(double dG, double T = 298.15);    // Kd = exp(dG / RT)
inline double deltaGFromKd(double kd, double T = 298.15);    // dG = RT.ln(Kd)
```

`src/chem/AdmetModel.h:167-187`, with the reasoning in the header comment at `:170-173`.
`tests/test_provenance.cpp:60-72` checks `dG = -9.0 kcal/mol -> Kd ~ 2.5e-7 M`, the
standard state `Kd = 1 M -> dG = 0`, and that temperature enters through RT.

The docking side is tiered accordingly. `DockJobResult::provenance` replaced what used to
be a `bool real` (`src/contracts/IDockingBackend.h:51-57`), and `fromEngine()` is now
literally `provenance == Provenance::Model` (`:72`). Backends must degrade to
`Provenance::Heuristic` rather than throw when their engine is unavailable (`:80-81`).

The docking panel then does the only honest thing (`src/ui/Panels.cpp:1052-1062`):

```cpp
drawQuantity("Best affinity",
             d.fromEngine()
                 ? makeQuantity(d.bestAffinity(), "kcal/mol", d.affinitySpread,
                                Provenance::Model, d.engine)
                 : makeQuantity(d.bestAffinity(), "", 0.0, Provenance::Heuristic,
                                "descriptor estimate - rank ordering only"));
```

A real dock is `Model` - a constructed pose, not a measurement - and keeps kcal/mol with
the pose spread as its error bar. The descriptor fallback is `Heuristic` and therefore
**cannot** keep the unit: writing `"kcal/mol"` on that branch throws.
`tests/test_provenance.cpp:95-118` asserts precisely that, which is the mechanism that
stops a fallback score from ever being presented as an energy.

## Worked example: the rehabilitated bioavailability model

The old code computed `firstPassSurvival = 1 / (1 + burden)` where `burden` accumulated
invented coefficients (`catechol += 15.0`, `ester += 2.20`, ...). That is algebraically the
well-stirred hepatic model with the units filed off. It is now written as the model it
always was (`src/chem/AdmetModel.h:125-154`):

```
F_H  = Q_H / (Q_H + fu.CLint)
CL_H = Q_H * fu.CLint / (Q_H + fu.CLint)
```

- `Q_H = 90.0` L/h, adult human hepatic blood flow - a population average, not a patient
  (`HepaticAssumptions::hepaticBloodFlowLPerH`, `src/chem/AdmetModel.h:61-63`), and stated
  as an assumption in the UI.
- `fu.CLint` is **not** predictable from structure, so it is an explicit parameter.
  `unboundIntrinsicClearanceLPerH` defaults to `-1.0`, meaning "derive an assumed value
  from the perceived structural liabilities" (`:64-66`), and `clIntMeasured` is true only
  when the value came from an experiment (`:67-68`).
- `assumedUnboundIntrinsicClearance` (`:97-123`) keeps the same ordinal coefficients, but
  expresses them as a multiple of hepatic blood flow so the model stays dimensionally
  coherent. The header states outright that any result built on them is `Heuristic`.

The arithmetic is now hand-checkable: `Q_H = 90`, `fu.CLint = 30` gives `F_H = 0.75` and
`CL_H = 22.5 L/h` (`tests/test_provenance.cpp:120-138`), and `CL_H < Q_H` for any input,
because an extraction ratio cannot exceed blood flow. `tests/test_provenance.cpp:151-154`
confirms the rehabilitation changed the units and the honesty, not the ranking: with the
assumed `fu.CLint` the result is still exactly `1 / (1 + fu.CLint / Q_H)`.

The rationale text spells the assumption out and shouts when it is assumed
(`src/chem/AdmetModel.h:156-165`):

> Hepatic availability under the stated assumptions: absorbed fraction (78%) x F_H (6%),
> well-stirred model with Q_H = 90 L/h and fu.CLint = 1350.0 L/h (ASSUMED from structural
> liabilities - rank ordering only, not a percentage to quote). Limited by catechol
> COMT/MAO first-pass metabolism.

And the absorption panel labels the row for what it is - "Hepatic availability (rank
order)", `Heuristic`, no unit (`src/ui/Panels.cpp:792-797`). It no longer says "predicted
bioavailability", because it never was one.

## Adding a new number: checklist

1. **Pick the tier from how the number was actually obtained**, not from how confident you
   feel. Retrieved with a citation or computed exactly from the input: `Measured`. A
   published model ran: `Predicted`, and you must have its benchmark error. A constructed
   artefact with no energy claim: `Model`. Anything else that only orders candidates:
   `Heuristic`.
2. **Build it with `makeQuantity`**, never by aggregate-initialising `Quantity`. If the
   call throws, the tier and the unit disagree and one of them is wrong.
3. **A `Heuristic` gets an empty unit.** If losing the unit makes the number look useless,
   that is the correct signal - it was useless as a physical value.
4. **Fill `source`.** An engine version, a citation, or a model plus benchmark. A quantity
   whose source is empty renders bare parentheses and is a review finding.
5. **Fill `error` when one exists**, in the same unit as `value`. For `Predicted`, the
   model's published benchmark error is mandatory.
6. **Derived from other quantities? Fold the tiers with `weakest()`** - never `std::min`,
   never a hand-rolled comparison.
7. **Missing a prerequisite? Return `notComputed("<the missing thing>")`** rather than a
   zero, a default, or a silently degraded estimate.
8. **Render it with `drawQuantity`.** Do not print a bare `%.2f` next to a label; the tier
   colour and the inline source are the deliverable.
9. **Test the boundary**, not the arithmetic alone: that the tier is what you claim, and
   that the illegal combination throws. `tests/test_provenance.cpp` is the pattern.
