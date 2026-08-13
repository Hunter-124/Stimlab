# Mechanism of action, off-target coverage, pathways, stacks and pharmacogenomics

Phase 7. One module, `biocad::RealMechanism` (`src/modules/MechanismModule.{h,cpp}`), behind
`IMechanismModule`, wired as `Services::mechanism`. Four panels: **Mechanism of Action**
(`Mechanism`), **Off-Target Panel** (`PanelScreen`), **Pathway Context** (`Pathways`) and
**Stack Checker** (`StackCheck`). Five agent tools: `retrieve_mechanisms`,
`screen_offtarget_panel`, `pathway_context`, `check_stack`, `pharmacogenomic_notes`.

Everything here **retrieves**. Nothing in this module derives a mechanism from a structure, a
fingerprint or a docking pose, and the DTOs in `src/data/Mechanism.h` offer no field in which to
put one.

## Licensing: the verdict for every source considered

Licensing is enforced in the implementation, not footnoted. Only the permitted rows below are
queried or bundled; the forbidden rows have no code path and no pack entry, and
`tests/test_mechanism.cpp` asserts that no mechanism pack names one as a data source.

| Source | Endpoint or artefact | Licence | Verdict and reason |
|---|---|---|---|
| ChEMBL | `https://www.ebi.ac.uk/chembl/api/data/mechanism.json?molecule_chembl_id=...`, `/data/target/{id}.json` | CC BY-SA 3.0 | **Used.** Queried live, cached under `%APPDATA%/BioCAD/cache/api`. Attribution is required and share-alike attaches to derived data, so every committed ChEMBL-derived pack carries the notice (see `assets/packs/mechanism/action-types.json`). |
| Reactome | `https://reactome.org/ContentService/data/mapping/UniProt/{acc}/pathways`, `/data/event/{stId}/ancestors` | CC0 | **Used.** The cleanest licence available; no obligation beyond courtesy attribution. |
| UniProt | `https://rest.uniprot.org/uniprotkb/{acc}` | CC BY 4.0 | **Permitted.** Accessions are the join key for pathway context. |
| RCSB PDB | `https://files.rcsb.org/download/{id}.pdb` | Open | **Permitted.** Already used by receptor preparation, which is what a panel screen docks into. |
| AlphaFold DB | `https://alphafold.ebi.ac.uk/api/prediction/{acc}` | CC BY 4.0, commercial use explicitly allowed | **Permitted.** |
| openFDA | `https://api.fda.gov/drug/{event,label}.json` | US Government work, public domain | **Permitted.** |
| CPIC | `https://api.cpicpgx.org/v1/...` and the published guidelines | CC0 | **Used, offline.** The phenotype vocabulary and CYP2D6 activity-score bands ship as a bundled pack; `pharmacogenomics()` therefore reports `networkAvailable = false` in every build and says so, rather than implying a live lookup. |
| FDA drug-interaction labeling tables | The FDA "Drug Development and Drug Interactions" substrate / inhibitor / inducer tables | US Government work, public domain | **Used.** The enzyme and transporter role assignments in `assets/packs/mechanism/interactions.json`. |
| KEGG | — | Non-academic use requires a paid licence; the REST API is academic-only | **NOT called, nothing bundled.** A deep link to a public `https://www.kegg.jp/pathway/hsa#####` page is offered instead: linking is not redistribution. `keggPathwayUrl()` returns a page URL and nothing in BioCAD fetches it. |
| STRING | — | Free for academics; commercial requires a licence | **Not used, not bundled.** |
| DrugBank | — | CC BY-NC 4.0, commercial requires payment | **Not used.** A non-commercial licence cannot back a portfolio/commercial-capable product. |
| PDBbind | — | Redistribution blocked | **Not used.** |
| BioLiP | — | Redistribution blocked | **Not used.** |
| Guide to PHARMACOLOGY | — | ODbL + CC BY-SA; API keys required from late August 2026 | **Blocked pending licence review.** ChEMBL plus Reactome cover the mechanism and pathway need without it. |

## The science gate

Every network-touching path is compiled behind the existing `BIOCAD_HAVE_SCIENCE` feature and
goes through the **one** HTTP client the assistant's web tools already own
(`biocad::agent::apiGet`, `src/agent/WebTools.cpp`) - one rate limiter, one cache, one
User-Agent. There is no second HTTP client.

With the science feature **off** (the default `windows` preset):

* `mechanisms()` returns `networkAvailable = false`, `retrievalAttempted = false`, no entries,
  and a warning stating that no transport is compiled in.
* `pathways()` likewise, plus its permanent "there is no pathway impact score" warning.
* `pharmacogenomics()` returns `networkAvailable = false` and the bundled CC0 notes.
* `screenPanel()` and `checkStack()` never needed a network: the docking hot path is
  cache-only and the interaction pack is local.

This is a build configuration, not an error, and no method fabricates content to fill the gap.

## What each method does, and what it refuses

### `mechanisms(compoundId)`

Resolves the compound's ChEMBL cross-reference from the loaded compound packs (or accepts a
`CHEMBL...` id directly), then retrieves ChEMBL's mechanism records: `action_type`,
`mechanism_of_action`, target name / accession / organism, and every reference. The
`action_type` is checked against the controlled vocabulary shipped in
`assets/packs/mechanism/action-types.json`; a value that is not in that list is reported
**verbatim** with a warning, never mapped onto a neighbouring value. The free-text
`mechanism_of_action` is displayed verbatim and is never parsed into a vocabulary, because
parsing it would invent structure the source does not have.

`coverageNote` is always populated and always says the same thing: an empty entry list is a
statement about this query and this database's curation, **not** about the compound. An
uncurated compound and a mechanism-free compound are indistinguishable here.

`InhibitionModality` is BioCAD's own orthogonal axis and stays `Unknown`: ChEMBL cannot express
it, so it is never guessed from the action type. Phase 10's enzyme-kinetics fits are its only
producer, and Phase 4's Cheng-Prusoff is its consumer.

### `screenPanel(molecule, panelId)`

Runs the **existing** docking module over every target in the named panel pack, in parallel on
the existing `workflow::JobSystem`. Three panels ship:

| Pack | `panelId` | Declared size | Rows enumerated |
|---|---|---:|---:|
| `safetyscreen44.json` | `safetyscreen44` | 44 | 44 |
| `safetyscreen87.json` | `safetyscreen87` | 87 | 44 |
| `cipa-currents.json` | `cipa-currents` | 5 | 5 |

`declaredSize` is what coverage is measured against, so an incomplete roster **inflates** the
unknown count instead of shrinking the denominator. SafetyScreen87 declares 87 and enumerates
only the 44 minimum-panel rows, so 43 targets are unknown to BioCAD before docking is even
attempted - and the coverage statement says exactly that.

The headline is `unscreened`, rendered first and largest. A row counts as screened only when a
**real engine** produced a pose (`Provenance::Model`); the labelled descriptor estimate is not a
screen, and recording it as one would move a target out of the unknown column on the strength of
a rank-ordering heuristic.

There is **no composite safety score and no cross-target comparison.** Receptor preparations,
box volumes and rotatable-bond penalties differ per row, so per-target scores are not on a
common scale. Each `PanelTargetResult` therefore carries its own `receptorPreparation` and
`boxDefinition`, and the DTO offers no rank, no aggregate and no comparison field.

The hERG margin is `measured IC50 / free Cmax`, flagged below 30-fold, and is computed **only**
from two user-supplied measurements (`RealMechanism::setHergInput`). With either missing it is
`notComputed` naming the missing one. A predicted hERG IC50 and any derived QT or TdP risk are
prohibited: there is no field for them and no code path that could fill one.

### `pathways(uniprotAccession)`

Reactome (CC0) UniProt-to-pathways plus per-event ancestors, rendered as a hierarchy with deep
links to `reactome.org`. **No pathway impact score.** No database supports propagating a docking
score, an affinity or an expression value through a pathway graph, and every report carries a
warning saying so rather than leaving the absence to be noticed.

### `checkStack(memberIds)`

`assets/packs/mechanism/interactions.json` carries 96 members - the FDA table's CYP and
transporter roles for drugs, plus more than 50 hand-curated supplement, food and lifestyle
entries (grapefruit furanocoumarins, St John's wort / hyperforin PXR induction, quercetin, EGCG,
piperine, curcumin, resveratrol, goldenseal, berberine, milk thistle, kava, licorice, harmala
alkaloids, tyramine-rich food, polyvalent cations, smoking, chronic ethanol, kratom and more).
Each member carries **its own citation**, enforced by the loader: a member without one is
dropped with an error.

Flags come from three derivations - inhibitor x substrate, inducer x substrate, and
mechanism-class pairs from `classRules` (MAO inhibition x monoamine release, MAO inhibition x
serotonin reuptake inhibition, P-gp inhibition x P-gp substrate, vitamin K supply x vitamin K
antagonism, polyvalent-cation chelation x chelation-susceptible absorption, and so on).

Every flag is a **mechanism with a citation and a direction**, never a severity score, and every
flag carries the pack's `boundaryNote` verbatim - a pack without that sentence is rejected at
load. A member that is not in the pack lands in `unknownMembers` and is reported as **not
screened**, because an unrecognised member silently dropped reads as a cleared stack.

### `pharmacogenomics(compoundId)`

CPIC (CC0) phenotype vocabulary - **UM, RM, NM, IM, PM**. "Extensive metabolizer" is deprecated
nomenclature: it appears nowhere, and a pack that ships it has that phenotype rejected with an
error. The CYP2D6 activity-score bands (PM = 0; IM 0 < AS <= 1.0; NM 1.25 <= AS <= 2.25;
UM > 2.25) follow Caudle and co-workers, *Clinical and Translational Science* (2020),
"Standardizing CYP2D6 genotype to phenotype translation", and ship as data so no panel restates
them from memory.

These are **conditional notes** - what a phenotype would imply. BioCAD interprets no genotype,
assigns nobody a phenotype, and emits no dose or dose adjustment.

## Agent tools

Each of the five states its boundary in its description **and** enforces it in the handler: the
returned JSON carries a `boundary` field, `screen_offtarget_panel` puts `unscreened` first in the
payload, and every handler refuses rather than guessing when the service or an argument is
missing. A description the handler does not back is decoration.

## Packs

`assets/packs/mechanism/`

| File | `kind` | Contents |
|---|---|---|
| `action-types.json` | `action-types` | ChEMBL `action_type` value names with BioCAD's own glosses; `declaredSize` 35, 33 transcribed. Carries the CC BY-SA notice. |
| `safetyscreen44.json` | `panel` | 44-row roster. |
| `safetyscreen87.json` | `panel` | Declares 87, enumerates 44, states the gap. |
| `cipa-currents.json` | `panel` | IKr, INa, ICaL, IKs, IK1. Kir2.1 has no receptor preset and is a declared gap. |
| `interactions.json` | `interactions` | 96 members, 11 class rules, one shared boundary note. |
| `pharmacogenomics.json` | `pharmacogenomics` | CPIC phenotypes, CYP2D6 bands, 12 conditional notes. |

Schema version 1 for every kind. An unknown version is an error surfaced in the panel, never a
silent skip.

### Known coverage gap in the target packs

`assets/packs/safety-offtarget.json` tags only 22 targets with `"panels": ["safetyscreen44"]`;
the other 22 rows of the roster live in `cns-monoamine.json` and `analgesics-otc.json` without
the tag. The mechanism panel pack is therefore the authority for the roster, and it resolves each
row to a receptor by target id. Adding a binding box to a pack target is all it takes for that
row to become screenable - no code change, and the unscreened count drops by one.
