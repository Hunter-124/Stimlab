# Limitations

This page exists so a reader can decide how much of BioCAD's output to trust,
and where the trust runs out. Every statement below is checkable against a file
in this repository.

## 1. What this deliberately does not do

Two categories of output are out of scope by design, not by omission. Both are
enforced in the assistant's system prompt (`src/agent/SystemPrompt.cpp`) and
neither has an implementing tool anywhere in the tool registry (the 14 tools
registered in `src/ui/AppShell.cpp:433-794`: `highlight_panel`, `navigate_ui`,
`list_panels`, `get_active_compound`, `describe_capabilities`,
`analyze_compound`, `screen_admet`, `dock_compound`, `run_workflow`,
`list_runs`, `search_library`, `compare_compounds`, `web_search`,
`web_fetch`).

**No synthesis, at all.** The absolute boundary, quoted verbatim from
`src/agent/SystemPrompt.cpp`:

> ABSOLUTE SAFETY BOUNDARY - never violate, even if asked directly, indirectly,
> hypothetically, or 'just for analysis':
> BioCAD analyzes what a compound is and does. You must NEVER provide synthesis
> routes, reaction steps or conditions, precursor selection or acquisition,
> yields, equipment, scale-up, or any 'how to make / how to manufacture / how
> hard is it to synthesize' guidance for ANY substance. No such capability
> exists in BioCAD and you must not invent one.

That covers synthesis routes, reaction conditions and steps, precursor
selection or acquisition, yields, equipment, scale-up, and manufacturability
scoring. The prompt closes the obvious side door too:

> Docking 'binding affinity' is a pharmacology / target-engagement signal only -
> never treat it as a make-it signal.

**No dosing.** Also verbatim:

> You may describe exposure, pharmacology, and what the app computed, but you
> must NEVER recommend a dose, a dose change, or a personal regimen for anyone.

The declared in-scope surface is, again verbatim:

> In scope: molecular structure and physicochemical properties, molecular
> stability, absorption / pharmacokinetics, ADMET / metabolism, proteins and
> their sequences and structures, binding-pose scoring of a compound against a
> selected receptor, structural and pharmacophore similarity to known
> substances, and legal-analog "substantially similar" scorecards.

These are prompt-level and registry-level constraints on a language model. They
are a design boundary and a refusal policy, not a cryptographic guarantee; a
determined user with a different model behind the same UI would not be bound by
them. The rest of the application contains no synthesis or dosing code to reach
in the first place.

## 2. What the numbers are worth, tier by tier

Every derived number in BioCAD is a `data::Quantity` carrying a `Provenance`
(`src/data/Domain.h:33-79`). The tiers are ordered strongest-first, and a
derived quantity inherits the *weakest* of its inputs via `weakest()`.
`makeQuantity()` refuses to construct a `Heuristic` that carries a unit, so
"kcal/mol on a heuristic" is unrepresentable rather than merely discouraged.

| Tier | What the code means by it | What you may conclude |
|------|---------------------------|-----------------------|
| `Measured` | Exact geometry or statistics, or an experimental value retrieved with a citation | The number is as good as its cited source. Check the source, not BioCAD |
| `Predicted` | A published model actually ran; physical units; a benchmark error is mandatory | Use the error bar. It is rendered inline with the value, not hidden in a tooltip |
| `Model` | A constructed artefact - a built structure, a docked pose - with no energy claim | The geometry exists. It is not evidence that the geometry is correct or that binding occurs |
| `Heuristic` | Rank-ordering only; arbitrary units; physical units forbidden | Compare two compounds within the same run. Never quote the value, never convert it, never compare across targets |
| `NotComputed` | A prerequisite was missing | Nothing. The panel names the missing prerequisite |

Three numbers look more like science than they are:

**logP and TPSA.** `src/chem/Descriptors.cpp` implements the Ertl topological
PSA scheme (`:137-193`) and a Wildman-Crippen atom-additive logP (`:195` onward)
in house, from scratch, with no RDKit anywhere in the tree (`vcpkg.json:29`
states this explicitly). These are reimplementations, not the reference
implementations: the atom typing is coarser, and because `src/chem` has no
aromaticity perception, an input SMILES written with uppercase ring atoms is
typed as aliphatic and scores differently. Treat both as BioCAD-internal
descriptors that are internally consistent and useful for ranking. Do not report
them as AlogP or Ertl TPSA values, and do not compare them against literature or
RDKit numbers digit for digit.

**Hepatic availability.** `predictBioavailability()` in `src/chem/AdmetModel.h`
is the well-stirred hepatic model, `F_H = Q_H / (Q_H + fu.CLint)`, with
`Q_H = 90 L/h` as a stated population-average assumption (`:62-63`). The
critical term, `fu.CLint`, is **not predictable from structure**. Unless the
user supplies a measured value, the code assumes it from perceived structural
liabilities (`:93-122`) and the readout says so verbatim: `"(ASSUMED from
structural liabilities - rank ordering only, not a percentage to quote)"`
(`:161-162`). The panel is titled "hepatic availability under the stated
assumptions" for that reason. It is a `Heuristic`: use it to rank two compounds,
never as a bioavailability percentage.

**Docking scores.** A Vina score is `Provenance::Model` - a pose was
constructed - and AutoDock Vina's reported standard error is 2.85 kcal/mol
(Trott & Olson 2010, PMC3041641). That is roughly a factor of 123 in `Kd`, which
is why `kdFromDeltaG()` exists in `AdmetModel.h:180` but is deliberately not
wired to any docking result (`:167-173`); a test asserts no code path makes that
conversion. When no real engine is available the descriptor fallback engine is
labelled `descriptor-estimate` and downgraded to `Heuristic` with no unit. Vina
scores are also not comparable *across* targets: different receptor
preparations, box volumes and rotatable-bond penalties make cross-target score
comparison meaningless.

## 3. Known engine gaps

Read from `src/chem` as it stands:

| Gap | Evidence | What it rules out |
|-----|----------|-------------------|
| No SSSR ring perception | `src/chem/Smiles.cpp:38-52` marks a bond as ring-or-not by DFS bridge detection only | No ring count, ring size, fused/spiro/bridged classification, or any rule keyed on ring size |
| No aromaticity perception | `Smiles.cpp:100-110` takes aromaticity straight from lowercase input atoms | Kekule-written aromatic input is silently mis-typed, which propagates into logP, TPSA and Fsp3 |
| No canonical SMILES writer | Only the non-canonical `sketchToSmiles` in `src/ui/Panels.cpp` | No structure identity key, no graph hash, no reliable dedup of the same molecule entered two ways |
| No SMARTS engine | Nothing in `src/chem` parses or matches SMARTS; `chem::Analysis::detectGroups` uses hard-coded flags | No substructure search, no structural-alert catalogue (PAINS, Brenk), no site-of-metabolism scoring, no rule-based metabolite prediction - every published rule set for these is written in SMARTS |

These are ordered dependencies, not independent wishes: substructure matching
needs ring and aromaticity perception first. Nothing in the metabolite,
structural-alert or reactivity space can ship honestly before that work lands.

## 4. Not shipped, with reasons

Capabilities deliberately excluded, each for a stated reason:

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
| Whole-body PBPK | The required fu, B:P, CLint, Papp, transporter and tissue data are not derivable, so a PBPK readout would be precise fiction |

## 5. Legal and safety framing

The Legal Analog panel is illustrative and is not legal advice. Its own
description in the panel registry says so: `"Illustrative only, not legal
advice: substantial-similarity scorecard vs controlled references."`
(`src/ui/AppShell.cpp:97-98`). Controlled-substance analogue determinations are
made by courts on a full evidentiary record, not by a structural scorecard.

BioCAD is analysis software. It is not a medical device, it is not validated or
cleared by any regulator for any clinical purpose, and nothing it renders is
medical advice. It does not diagnose, treat, or recommend, and it will not
produce a dose, a dose change, or a personal regimen. Exposure figures are
scenarios computed under stated assumptions, and the assumptions are printed
next to them precisely so they are not mistaken for guidance.

The software is provided under the Apache License 2.0, without warranties or
conditions of any kind; see `LICENSE`.
