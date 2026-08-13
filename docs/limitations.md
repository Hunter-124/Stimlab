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

**logP and TPSA.** `chem::tpsa()` (`src/chem/Descriptors.cpp:146-215`) is an
in-house implementation of the Ertl topological PSA scheme; `chem::crippenLogP()`
delegates to `chem::crippen()` (`src/chem/Crippen.{h,cpp}`), a full
Wildman-Crippen atomic-contribution implementation whose parameter table - 106
ordered pattern entries covering the method's 67 heavy-atom classes - is data in
`assets/packs/descriptors/crippen.json`, with the four hydrogen classes derived in
`Crippen.cpp` because the graph carries no explicit H atoms. There is no RDKit
anywhere in the tree (`vcpkg.json:29` states this explicitly), so both are
reimplementations of published methods. logP and molar refractivity are now
verified EXACT against a reference implementation; TPSA is within 1.6 A^2 of it on
the N,O-only convention; logP is faithful to *its method*, and that method carries a
published RMS of about 0.67 log units against experiment, which every `Quantity`
built from it must display. Measured deviations, and the rules that follow, are in
section 3 below and in [cheminformatics.md](cheminformatics.md). Do not report
either as "AlogP" or as an "Ertl TPSA" value from another toolkit, and do not
diff them digit for digit against one.

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

## 3. Descriptor fidelity

`src/chem` now has SSSR ring perception (`src/chem/Rings.cpp`), graph
aromaticity perception (`src/chem/Aromaticity.cpp`), a canonical SMILES writer
plus graph hash (`src/chem/Canonical.cpp`) and a SMARTS parser + VF2 matcher
(`src/chem/Smarts.cpp`), so the earlier "no rings, no aromaticity, no canonical
form, no matcher" gap list is closed. What remains is a fidelity question, and it
is measurable. Values below were produced by compiling the real sources on a
Linux host and running them (harness and verbatim output in
[cheminformatics.md](cheminformatics.md#measured-evidence)):

| Quantity | Compound | BioCAD | Reference | Deviation |
|---|---|---|---|---|
| TPSA (A^2) | ibuprofen | 37.30 | 37.3 | 0.00 |
| TPSA (A^2) | aspirin | 63.60 | 63.6 | 0.00 |
| TPSA (A^2) | dopamine | 66.48 | 66.5 | 0.02 |
| TPSA (A^2) | caffeine | 60.26 | 58.44 | 1.82 |
| logP | benzene | 1.687 | 1.69 | 0.003 |
| logP | aspirin | 1.310 | 1.31 | 0.000 |
| logP | ibuprofen | 3.073 | 3.97 | 0.897 |
| logP | caffeine | -1.029 | -0.07 | 0.959 |
| logP | dopamine | 0.599 | -0.98 | **1.579** |

TPSA reference values are the ones PubChem publishes for these compounds
(PubChem's computed TPSA is the Ertl scheme). The logP reference values are
**experimental** octanol-water logP, not Crippen-computed values. No page number
or DOI is quoted for either set, because neither was read out of a primary
publication's tables.

(Those TPSA deviations are against PubChem's published values; the RDKit
atom-for-atom comparison immediately below reports caffeine at 1.56 A^2 rather
than 1.82 A^2 because it is a different reference, not a different calculation.)

**logP and molar refractivity are exact against the reference implementation.**
RDKit 2026.03.5 was installed on the development machine purely as an ORACLE (it is
not a dependency and is not in `vcpkg.json`) and every one of the 69 shipped library
compounds was compared atom-for-atom:

| Quantity | Compounds compared | Mismatches |
|---|---:|---:|
| Wildman-Crippen logP | 69 | **0** (tolerance 0.005) |
| Molar refractivity | 69 | **0** (tolerance 0.02) |
| Molecular formula | 69 | **0** (exact string) |

That is a much stronger claim than agreeing with a handful of literature values, and
it is the claim worth making: the implementation *is* the published method. Whether
the published method is right about a given compound is a separate question,
answered immediately below.

The formula result is not incidental. It was 69/69 only after a real bug was found
in `assignImplicitHydrogens` (`src/chem/Smiles.cpp`): aromatic bonds were counted as
order 1.5, which gave thiophene's sulfur a phantom implicit hydrogen and made
`c1ccsc1` and `C1=CSC=C1` different molecules.

**TPSA is close to the reference but not exact, and it has a convention choice.**
Cross-checked atom-for-atom against RDKit 2026.03.5 over all 69 shipped library
compounds:

| Convention | Mismatches vs RDKit (>0.1 A^2) | Worst deviation |
|---|---:|---:|
| N,O only (BioCAD default) | 10 of 69 | 1.56 A^2 (caffeine) |
| Including S and P | 11 of 69 | 13.5 A^2 (N-acetylcysteine) |

Two things came out of that check and both are now explicit in the code:

1. **The default sums N and O only.** Ertl's paper also tabulates sulfur and
   phosphorus, and an earlier revision included them - which put famotidine at
   237.75 A^2 against the reference 175.83 A^2, a 62 A^2 gap that straddles Veber's
   140 A^2 oral cutoff. Every published TPSA threshold was derived on the N,O-only
   sum, so that is what `tpsa()` returns;
   `tpsaIncludingSulfurAndPhosphorus()` is a separate function for the other
   convention, and a readout must say which one it used.
2. **The 10 residual mismatches are all N-typing, and all small**: deviations are
   exact multiples of 0.52 A^2 (celecoxib and indomethacin -0.52; phenazone,
   ondansetron, theobromine and theophylline -1.04; caffeine -1.56) or +0.22 A^2
   (NMN, thiamine, berberine). They trace to the aromaticity model over-assigning
   aromaticity to fused pyrimidinedione and pyrazolone rings - a limitation
   `src/chem/Aromaticity.h` already documents in its own terms - which then selects
   an aromatic rather than an aliphatic nitrogen contribution. It is a known model
   boundary surfacing in a descriptor, not an arithmetic error, and it is under
   2 A^2 in every case.

**logP is faithful to its method, and its method has real error.**
`chem::crippenLogP()` delegates to `chem::crippen()`
(`src/chem/Crippen.{h,cpp}`), which implements Wildman & Crippen's
atomic-contribution scheme, whose parameter table - 106 ordered pattern entries
covering the method's 67 heavy-atom classes, with H1-H4 derived in `Crippen.cpp` -
is data in `assets/packs/descriptors/crippen.json`. Against the ten-compound
accuracy set in `tests/test_chem.cpp:107-132` the measured error is
**MAE 0.349, RMS 0.430** log
units, inside the method's own published RMS of about 0.67 log units - which
`chem::crippenCitation()` states in every `Quantity` built from the number, and
which the suite enforces in the method's own terms (+/-1.35 per compound, aggregate
RMS < 0.67). Exact-value oracles against the reference implementation live in
`tests/test_chem_crippen.cpp`. The larger residuals above (caffeine, dopamine) are
the additive method's known weaknesses on fused polar heteroaromatics and free
catechols, documented as such in the test file, not implementation bugs. A logP
from BioCAD must therefore be rendered with that RMS attached; it is not a value
to quote to two decimals, and it is not "AlogP".

Two further honesty rules for this section:

**An in-house reimplementation is labelled as such.** Any descriptor, score or
metric BioCAD computes with its own code is named as a BioCAD-internal quantity
wherever it is surfaced - never as the toolkit-standard quantity of the same
name. "AlogP", "Ertl TPSA", "ECFP4" and "Morgan similarity" are other people's
validated implementations; BioCAD's `tpsa`, `crippen` and `morganFingerprint`
are its own implementations of the same published ideas and are comparable only
to themselves. A reimplementation is also a place a bug can hide that a
widely-used library would not have, and that is the honest reason for the label.
The one exception is recorded in `NOTICE` section 4a: the Crippen parameter table
is transcribed from RDKit's BSD-3-Clause `Data/Crippen.txt`, so those *values* are
the reference values, under that licence.

**Descriptors no longer depend on how the SMILES was typed.** `chem::crippen()`,
`chem::detectGroups()`, `chem::screenAlerts()` and `chem::tpsa()` each perceive
rings and aromaticity on a private copy, so no caller can hand them an
unperceived graph and none has to remember the rule - which matters because
`src/modules/RealBackend.cpp:41-53` passes the raw `chem::parseSmiles` output and
`assets/packs/cns-monoamine.json:10` stores caffeine in Kekule form. Measured:
caffeine TPSA is 60.260 from both spellings (difference 0.0e+00), logP -1.029 from
both. Before that fix the Kekule spelling gave TPSA 56.22, a 4.04 A^2 swing from
spelling alone; the fix went into the descriptor functions rather than into their
~15 call sites.

(Which atoms TPSA counts, and why thiophene therefore reads 0.00, is item 1
above.)

What the engine still cannot do at all - no stereochemistry anywhere in
`chem::Molecule` (so enantiomers share a canonical string and a graph hash), no
tautomer handling, no isotopes, and a 3D embedder that is distance geometry and
explicitly not a protein method - is enumerated with its evidence in
[cheminformatics.md](cheminformatics.md).

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
