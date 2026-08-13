# Reaction networks, chemical kinetics, flux and pathway enrichment

Phase 14. Everything here is in `src/sim` (target `biocad_sim`), reached through
`ISimulationModule`, `IEnrichmentModule` and - only in a `BIOCAD_ENABLE_FBA` build -
`IFluxModule`. Panels: **Reaction Network** (`Networks`), **Metabolic Flux** (`Flux`),
**Pathway Enrichment** (`Enrichment`). Agent tools: `simulate_reaction_network`,
`fit_stability_kinetics`, `pathway_enrichment`, and `metabolic_flux` when FBA is on.

The single rule this phase exists to enforce: **a simulation is a property of the
mechanism and the numbers you entered.** A time course is not a measurement of a cell,
a computed objective is not a growth rate measured in an organism, and no docking score
may ever enter a flux, a pathway or a network number.

## One representation

`sim::Network` compiles a `NetworkSpec` into a dense species x reaction stoichiometric
matrix `S` and evaluates `dc/dt = S.v(c)`. Five rate laws: mass action, reversible mass
action, Michaelis-Menten, reversible Michaelis-Menten, Hill. Boundary species are
clamped pools with an exactly-zero derivative row.

**The analytic Jacobian never divides by a concentration.** For a mass-action rate
`v = k * prod(c_i^a_i)` the tempting form is

```
dv/dc_j = v * a_j / c_j
```

which is algebraically right and numerically indefensible: the moment a species is
depleted it divides by zero, and the derivative it should return there is generally
finite and nonzero. `Network::rateJacobian` recomputes the product over the *other*
species and multiplies by `a_j * c_j^(a_j - 1)`, which is exact at `c_j = 0`. The test
suite compares the analytic Jacobian with a central difference at a state where a
species is exactly zero, so the shortcut cannot creep back in.

### Structure from the matrix

* **Conservation laws** are the LEFT null space of `S` restricted to non-boundary rows
  (Eigen `FullPivLU::kernel`). Each basis vector is a conserved moiety, reported with a
  readable label (`y1 + y2 + y3`), and its value over time is stored in every
  `TimeCourse` as an audit. A boundary species is excluded because a moiety flowing into
  a clamped pool is not conserved.
* **Wegscheider cycle conditions** are the RIGHT null space. For a cycle supported
  entirely on reversible mass-action steps, `prod (kf/kr)^k_i` must be 1. A violation is
  reported in `NetworkSpec::wegscheiderViolations` as an **error in the parameters**, and
  `integrate()` refuses to run: rate constants that do net work around a closed cycle
  produce a smooth, plausible, thermodynamically impossible steady state.

## Solvers

| method | what it is | when |
|---|---|---|
| `rk4` | `numeric::rk4Integrate`, fixed step, explicit | smooth, well-scaled networks |
| `rosenbrock` | ROS3P (Lang & Verwer), 3 stages, 3rd order, L-stable, adaptive | anything stiff |
| `gillespie` | exact stochastic simulation algorithm | low copy number, where the variance *is* the result |
| `tau-leap` | explicit tau-leaping, Cao-Gillespie-Petzold step selection | high copy number |

**One LU factorization per Rosenbrock step.** `(1/(gamma*h) I - J)` is factorized once
and all three stages solve against it; the test suite asserts
`jacobianEvaluations == acceptedSteps + rejectedSteps`.

**Why the error estimate is against the first-order solution.** For this table, applied
to a linear problem, the stage vectors satisfy `(a21 - 1) k1 + k2 == 0` identically, and
every second-order set of weights is forced to the same `m3`. Any second-order embedded
solution therefore differs from `m` only along a direction that vanishes on a *linear*
problem - so its error estimate is exactly zero for `A -> B`, and the step would grow
without bound on precisely the problems whose analytic solutions we check. The estimate
is taken against the provably first-order solution `y + (1/gamma) k1` (the
linearly-implicit Euler step embedded in stage one). It is `O(h^2)` on linear and
nonlinear problems alike, which makes the step control **conservative**: asking for
`relTol = 1e-8` on `A -> B` delivers about `1e-13`.

**A Rosenbrock step conserves exactly, not to the tolerance.** With `M = aI - J` and a
conservation law `y` (so `y'S = 0` and hence `y'J = 0`), `y'M = a y'`, so
`y'M^-1 = y'/a` and every stage vector satisfies `y'k_i = 0`. The reported
`worstConservationDrift` is therefore a genuine audit of the arithmetic - Robertson
comes back at ~1e-14 - and not a restatement of the tolerance.

**Every result carries a `SolverReport`**: method, both tolerances, accepted and
rejected steps, Jacobian evaluations, non-negativity clips and CPU seconds. A nonzero
clip count is a WARNING in `TimeCourse::warnings`, never a silent repair: a solver that
quietly clamps a negative concentration is hiding a step-size failure.

**Stochastic state is a copy number, not a concentration.** The SSA counts molecules in
unit volume and applies the combinatorial factor `prod (n_i choose a_i) * a_i!`. A
saturable rate law is not an elementary reaction channel, so a network containing one is
**refused by name** rather than approximated by its mean-field rate - which would put a
deterministic rate inside an exact stochastic algorithm and destroy the variance it was
run for. Tau-leaping is approximate: its bias is bounded by its own `epsilon`, not by the
replicate count, and the suite measures that the bias tracks `epsilon` (2.4% at 0.03,
0.32% at 0.005 on a 200000-molecule decay) rather than pretending it vanishes.

Each replicate draws from its own `sim::Pcg64Dxsm` stream derived from the run seed and
the replicate index, so a replicate is reproducible alone and adding replicates never
perturbs the ones already computed.

## Chemical kinetics, and what was deleted

`predictPhWindow`, `predictThermalWindow` and `shelfLifeEstimate` are **gone** from
`src/chem/AdmetModel.h` and `StabilityReport`. They worked by adding invented per-group
interval bounds ("ester: pH 3.0-6.0", "catechol: store below 8 C"), averaging five
0-100 flag counts, and mapping the average onto one of four strings such as
`"~24 months @ 25C/60%RH"`. Every number in that chain was authored, and the output
looked exactly like a stability study.

`StabilityReport::shelfLife` is now a `Quantity`, and with only a structure in hand it
reads `notComputed("measured degradation rate constants at three or more temperatures
...")`. The Stability panel's thermal and pH rows print the **mechanism** and say
plainly that a window requires measurement.

What replaces them, in the Reaction Network panel's kinetics section and in
`fit_stability_kinetics`:

* `sim::arrhenius` - `k = A exp(-Ea/RT)`, fitted through `numeric::levenbergMarquardt`
  in `(ln A, Ea)` with the analytic Jacobian. Recovers `A` and `Ea` from noiseless data
  to better than 1e-9 relative.
* `sim::eyring` - `k = kappa (kB T/h) exp(dS/R) exp(-dH/RT)`, fitted in
  `(dH kJ/mol, dS J/(mol K))`, with the **joint 95% confidence ellipse** for `(dH, dS)`
  from the 2x2 covariance block. The two are strongly correlated (the isokinetic
  artefact), so two separate error bars overstate what the experiment determined.
* `sim::shelfLife` - `t = -ln(1 - f) / k(T)` from the Arrhenius fit, carrying the fit's
  relative prediction width as its error.
* `sim::phRate` - `k_obs = kH[H+] + k0 + kOH[OH-]`, with the minimum reported from its
  **closed form** `pH = 0.5 (pKw + log10(kH/kOH))` and `k_min = k0 + 2 sqrt(kH kOH Kw)`,
  not by scanning a grid (a scanned minimum inherits the grid spacing as a fake
  precision). `pKw = 14.0` is an assumption recorded in the result: at 37 C it is 13.62
  and the reported minimum moves by 0.19 pH units.

**Fewer than three distinct temperatures refuses to extrapolate.** Two points fit a
line exactly and say nothing about how wrong it is, so `extrapolationSupported` is
false, `predictedRateAt25C` is `notComputed`, and no shelf life is produced. Fewer than
three pH values likewise refuses a three-parameter pH-rate fit.

## SBML

`sim::importSbml` / `sim::exportSbml` implement the SBML **Core** subset over
`pugixml` (MIT). **libSBML is deliberately not linked, statically or otherwise**: it is
LGPL, and the relinking obligation is incompatible with distributing the single fully
static Windows binary this project ships.

Supported: `model`, `listOfCompartments`, `listOfSpecies`, `listOfParameters`,
`listOfReactions` with reactants / products / modifiers, and a `kineticLaw` whose
Content MathML is one of the five rate laws. Level 3 Version 2 is the export target;
**Level 2 is also read**, because the curated CC0 BioModels corpus is overwhelmingly
Level 2 and a Level-3-only reader could not open a single real fixture.

Refused **by name**, through the `error` out-parameter: `functionDefinition`,
`assignmentRule`, `rateRule`, `algebraicRule`, `event`, `constraint`,
`initialAssignment`, `stoichiometryMath`, an unsupported SBML level, any MathML operator
outside `times plus minus divide power`, and any well-formed law that is not one of the
five (typically an enzyme-modified `kcat * E * S / (Km + S)` where `E` is a species and
not a reactant). Silently dropping an assignment rule would change the model's meaning
and still produce a curve, which is worse than an error, so `importSbml` returns
`std::nullopt`.

### The committed fixtures, and a finding

`tests/fixtures/sbml/` holds three real curated BioModels entries, released **CC0**, so
no test touches the network:

| fixture | model | size |
|---|---|---|
| `BIOMD0000000052.xml` | Brands 2002, monosaccharide-casein Maillard kinetics | 11 species, 11 reactions |
| `BIOMD0000000050.xml` | Martins 2003, Amadori degradation | 14 species, 16 reactions |
| `BIOMD0000000035.xml` | Vilar 2002, circadian oscillator | 10 species, 16 reactions |

All three are pure mass action, all three import, all three integrate, and
`BIOMD0000000052` round-trips through the writer with every species, stoichiometry and
rate constant identical *and* the same trajectory to the last bit.

**Finding.** Most curated BioModels entries are **not** importable, and this is a
property of the corpus rather than a gap in the reader. A survey of the first 79 curated
entries found that the majority use `assignmentRule`, `event`, `functionDefinition`, or
enzyme-modified rate expressions such as
`uVol * V1 * MKKK / ((1 + (MAPK_PP/Ki)^n) * (K1 + MKKK))` (BIOMD0000000010) that are not
any of BioCAD's five closed-form laws. Six of the 79 were pure mass action. Supporting
the rest means either an expression evaluator with per-reaction ASTs or new `RateLaw`
enumerators; both change the `Systems.h` DTO contract and were left out of scope rather
than approximated. The reader names the exact construct in every refusal, which is the
useful behaviour in the meantime.

## Metabolic flux - `BIOCAD_ENABLE_FBA`, default OFF

**Which LP solver shipped, and why.** The plan offered HiGHS from vcpkg. HiGHS is not
present in this development environment, and this project builds only on Windows while
it is developed on Linux, so adding it would have meant shipping an LP path never once
executed by the person who wrote it. `src/sim/Flux.cpp` therefore contains a **dense
two-phase primal simplex written here**, pivoted by Bland's rule (FBA problems are
massively degenerate and Dantzig's rule can cycle on them forever; Bland's rule provably
terminates). It is checked against a hand-solved 2-variable LP, an infeasible LP and an
unbounded LP, and against a hand-solved toy metabolic model. It is **not** a
genome-scale solver: a dense tableau over a 2000-reaction reconstruction is out of
reach, and `FluxSolution::warnings` says so above 300 reactions.

**`balance()` is a gate, not advice.** Every reaction's two sides are compared element by
element through `chem::parseFormula`, plus net charge. A species with no formula is
**UNCHECKED, and unchecked counts as unbalanced** - never as "probably fine". `fba()`
refuses to run until the whole network balances, because a reaction that does not
conserve elements can carry flux out of nothing and an objective optimised over it
measures the arithmetic, not the metabolism. (`SpeciesSpec` gained `formula` and `charge`
fields for this; the DTO as landed had nowhere to state either.)

Provided: `fba`, `fva` (min and max per reaction with the objective pinned at a stated
fraction of its optimum), `parsimonious` (minimum total `|flux|` among the optima, with
`|v|` linearised exactly), and `deletions` of order 1 or 2. **Every `FluxSolution`
carries its objective reaction and every bound in force**, not only the bounds the caller
passed, because a flux distribution without the medium it was allowed is not
interpretable.

The toy model in the suite is hand-solvable: uptake capped at 10, two pyruvate per
glucose, a maintenance demand of 3 on the lactate branch, so the acetaldehyde export
optimum is exactly `2*10 - 3 = 17`. Deleting the maintenance reaction *raises* the
objective to 20, which is the kind of result a reader should be able to check by hand.

## Sensitivity and metabolic control analysis

`sim::steadyState` integrates with the Rosenbrock solver and then polishes with a
least-squares Newton step (the Jacobian is singular by exactly the number of
conservation laws, so a plain solve would fail). `sim::sweep` scales one reaction's rate
over a list of factors. `sim::controlAnalysis` computes flux and concentration control
coefficients by central log-log difference on a rate-scaling perturbation that moves
**both** directions of a reversible step, so the equilibrium constant is untouched and
the perturbation is an enzyme-amount change rather than a thermodynamic one.

Both theorem residuals are reported and shown first in the panel:

* **summation** - flux control coefficients sum to 1, concentration control coefficients
  to 0;
* **connectivity** - `sum_k C^J_{r,k} * eps_{k,j} = 0` for every species.

These are theorems, not fit statistics. A nonzero residual means the steady state was not
reached, the Jacobian is wrong, or the perturbation was too large - and reporting it is
how that gets caught instead of shipping a confident table of wrong coefficients. The
suite checks a two-reaction network whose coefficients are hand-computable (2/3, 1/3 for
flux; 2/3, -2/3 for concentration; elasticities -0.5 and 1).

## Pathway enrichment

Hypergeometric upper tail with log-gamma binomials, then **Benjamini-Hochberg** step-up
q-values with monotonicity enforced. The implementation is checked against the worked
example in Benjamini & Hochberg (1995), *J R Statist Soc B* 57:289-300: their 15
p-values, four rejections at alpha = 0.05, and the exact step-up values.

**The background is a required argument with no default.** A hypergeometric p-value is a
function of the universe the query was drawn from; using "every gene in the database"
when the experiment could only detect a few thousand transcripts inflates every result.
`enrich()` intersects both the query and each pathway with the background, reports query
identifiers that are not in the background (they cannot be tested), and records the
universe size, the number of query identifiers and the number of pathways tested - the
last of which is needed to recompute the q-values at all.

### The shipped pack: what is and is not in it

`assets/packs/pathways/reactome-human.gmt`, 1.0 MB, **2868 sets**.

* **Source**: `https://reactome.org/download/current/ReactomePathways.gmt.zip`,
  Reactome release 97, retrieved 2026-08-13.
* **Licence**: CC0 1.0 Universal. Reactome places its data in the public domain, so the
  file is redistributable verbatim with no share-alike obligation. It is attributed
  anyway - and the release is stored on every hit - because a pathway result without its
  database release cannot be reproduced.
* **Included**: every human (`R-HSA-`) pathway in release 97 with at least one gene
  product, as pathway name, stable id and HGNC symbols. Reactome's GMT is
  name-then-id, the reverse of the Broad convention; `sim::parseGmt` detects the
  `R-HSA-` stable id and swaps them.
* **NOT included**: other species; the pathway **hierarchy** (parent/child relations),
  so a hit on "Metabolism" and a hit on "Glycolysis" are reported as two independent
  tests and the redundancy is visible rather than corrected; reaction-level detail;
  small-molecule participants; and Gene Ontology terms. KEGG is absent deliberately -
  its REST API is academic-use-only and its data cannot be redistributed here. STRING
  and DrugBank are absent for the same class of reason.

There is **no pathway impact score**. No database supports propagating a rate, a flux or
a docking score through a pathway graph, and such a number would be fabrication with a
scientific veneer.

## Graph metrics

Degree, connected components (BFS), **Brandes betweenness** (halved, because the
accumulation counts each unordered pair twice on an undirected graph) and **Louvain**
communities with the partition's modularity reported beside them. Betweenness is checked
against a hand-counted five-node path graph (0, 3, 4, 3, 0) and Louvain against two
disjoint triangles (two communities, modularity exactly 0.5).

**Every edge must carry evidence.** An edge with an empty `evidence` string or
`Provenance::NotComputed` is **dropped** and named in `GraphMetrics::warnings`; so is a
self-loop, which changes no shortest path and inflates a degree. Degree, betweenness and
community membership are only as trustworthy as the edge set, and an interaction network
assembled from unattributed edges produces confident numbers about nothing. The
Enrichment panel builds its edges from the loaded reaction network, each carrying the
reaction id as its evidence and `Provenance::Model` - the network is a written-down
mechanism, not a measurement.
