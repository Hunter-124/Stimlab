# BioCAD

**A native workstation for molecular, protein, and pharmacological analysis.**

![BioCAD dashboard](docs/media/dashboard.png)

BioCAD is a single C++20 Windows desktop application that takes a compound from a SMILES
string to a scored binding pose without a browser, a Python runtime, a Docker daemon, or a
localhost server. Ligand and receptor preparation, an in-house cheminformatics engine, real
AutoDock Vina docking, ADMET and PK screening, a re-runnable workflow DAG, a SQLite run
history, and a tool-driven in-app assistant all live in one process and one executable.

Its distinguishing property is not the feature list. It is that **every derived number
carries the reason it may be trusted**, and the type system refuses to let a rank-ordering
score wear physical units. See [docs/provenance.md](docs/provenance.md) - that document is
the point of the project.

> **Research software notice.** BioCAD is for lawful research and educational use.
> Computational output is a hypothesis, not a clinical, safety, efficacy or regulatory
> conclusion. It does not recommend doses, and it is not a medical device.

## Honest numbers, by construction

Most tools of this kind print a number to two decimal places and let the reader supply the
confidence. BioCAD attaches a tier to every derived value:

| Tier | Meaning | Rendered |
| --- | --- | --- |
| `Measured` | Exact geometry or statistics, or an experimental value with a citation | green |
| `Predicted` | A published model actually ran; physical units; benchmark error mandatory | blue |
| `Model` | A constructed artefact (a built structure, a docked pose) with no energy claim | purple |
| `Heuristic` | Rank ordering only; arbitrary units; physical units are **forbidden** | amber |
| `NotComputed` | A prerequisite was missing, and the value names which one | grey |

The rule is enforced, not documented: `makeQuantity()` in `src/data/Domain.cpp` throws when a
`Heuristic` is handed a unit. So a real Vina result is `Model` and reads `-9.30 kcal/mol`,
while the descriptor fallback is `Heuristic` and *cannot be expressed* in kcal/mol at all.
That is why the app owns `kdFromDeltaG()` yet deliberately never applies it to a docking
score: Vina's reported standard error of 2.85 kcal/mol is a factor of about 123 in Kd, so a
nanomolar affinity derived from a docked pose would be a fabricated precision.

The error bar and the tier are part of the value, never a tooltip.

## Capabilities

| Area | What it provides | |
| --- | --- | --- |
| **Structure and properties** | SMILES parsing, molecular graph, descriptors, formula/MW/logP/TPSA, 3D conformer embedding, and a DX11 molecular viewport | ![Structure](docs/media/structure.png) |
| **Docking** | Ranked poses from real AutoDock Vina / smina, with Vina-GPU (OpenCL) and a first-party CUDA backend; auto-provisioned engines and on-demand receptor prep | ![Docking](docs/media/docking.png) |
| **Absorption and PK** | Absorbed fraction and hepatic availability under the well-stirred model, with the assumption set printed beside the number | ![Absorption](docs/media/absorption.png) |
| **ADMET and metabolism** | Structure-derived liability perception driving metabolism, stability and safety endpoints | ![Metabolism](docs/media/metabolism.png) |
| **Similarity and analogs** | Fingerprint and pharmacophore similarity against the loaded catalog, plus an analog sketch/compare workflow | ![Similarity](docs/media/similarity.png) |
| **Legal analog** | Substantial-similarity scorecard against controlled references - illustrative, explicitly not legal advice | ![Legal](docs/media/legal.png) |
| **Workflows** | A cancellable, content-cached prep-to-dock DAG with a live execution view | ![Workflows](docs/media/workflows-dag.gif) |
| **Data packs** | The compound and target catalog is versioned JSON, not C++; drop a pack in and it appears without a rebuild | ![Presets](docs/media/presets.png) |
| **Assistant** | A tool registry the model drives: read properties, dock, run a workflow, navigate and highlight the UI | ![Settings](docs/media/settings.png) |

## The catalog is data

There is no hard-coded compound table. Packs are versioned JSON documents under
`assets/packs/`, overridden by pack id from `%APPDATA%\BioCAD\packs`:

```jsonc
{
  "schemaVersion": 1,
  "id": "my-pack",
  "title": "My compounds",
  "compounds": [
    { "id": "ibuprofen", "name": "Ibuprofen", "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
      "xrefs": { "chembl": "CHEMBL521" } }
  ],
  "targets": [
    { "id": "PTGS2", "name": "COX-2 (prostaglandin G/H synthase 2)", "pdb": "5KIR",
      "box": { "cx": 0.0, "cy": 0.0, "cz": 0.0, "sx": 22.0, "sy": 22.0, "sz": 22.0 } }
  ]
}
```

An unknown `schemaVersion` is a load error shown in the Presets panel, never a silent skip. A
binding-site box requires a real PDB entry; a target without one is listed as an honest
coverage gap rather than given an invented site. Four packs ship built in: 67 compounds and 59
targets, 29 of them dockable. Full schema in [docs/packs.md](docs/packs.md).

## Architecture

Modules are compile-time dependency injection through a struct of interface pointers
(`src/contracts/Services.h`), populated by `RealBackend::services()` and mirrored by a
`RealBackend::services()`. There is no plugin loader, and there is deliberately no fake
backend: the test suite runs the shipping implementation, because a suite that validates a
double while the product ships the original proves nothing.

```text
┌───────────────────────────────────────────────────────────────────────────┐
│ UI thread: Win32 + Dear ImGui + ImPlot + DirectX 11                       │
│ navigation · panels · molecular viewport · assistant controls             │
└──────────────────────────────┬────────────────────────────────────────────┘
                               │ commands / snapshots / events
┌──────────────────────────────▼────────────────────────────────────────────┐
│ Worker services: workflow DAG · job system · docking · chemistry · agent  │
│ SQLite run store · artifact storage · provisioning · data packs           │
└──────────────────────────────┬────────────────────────────────────────────┘
                               │ native libraries and subprocess engines
┌──────────────────────────────▼────────────────────────────────────────────┐
│ Vina · Vina-GPU · CUDA backend · SQLite · JSON · WIC · optional curl      │
└───────────────────────────────────────────────────────────────────────────┘
```

Preparation, provisioning and docking run off the UI thread. Details, including the checklist
for adding a module, are in [docs/architecture.md](docs/architecture.md); the engine
provisioning, receptor-prep and scoring story is in [docs/docking.md](docs/docking.md).

### Documentation

| Document | What it covers |
|---|---|
| [docs/architecture.md](docs/architecture.md) | Layers, the `Services` compile-time DI seam, the job system and DAG cache, thread discipline |
| [docs/provenance.md](docs/provenance.md) | The five provenance tiers, the `Quantity` type, and why a `Heuristic` cannot carry a unit |
| [docs/cheminformatics.md](docs/cheminformatics.md) | The in-house chem engine: why there is no RDKit, SSSR rings, graph aromaticity, canonical SMILES and `graphHash`, the SMARTS parser and VF2 matcher, measured descriptor fidelity, and what the engine still cannot do |
| [docs/packs.md](docs/packs.md) | Pack JSON schema, resolution order, and how to author one |
| [docs/docking.md](docs/docking.md) | Engine locate/provision order, receptor prep, PDBQT, the box marker, what a Vina score is and is not |
| [docs/pkpd.md](docs/pkpd.md) | Dose-response fits, Cheng-Prusoff, Schild, the PK engine and occupancy |
| [docs/protein.md](docs/protein.md) | Structure and sequence I/O, alignment, superposition, structure scoring |
| [docs/data-sources.md](docs/data-sources.md) | Every external service, its endpoint, its licence, and what may be cached or committed |
| [docs/limitations.md](docs/limitations.md) | The do-not-ship list, descriptor fidelity, and the honesty rules |

## Quick start

BioCAD is built and supported on **Windows x64**. Install Visual Studio 2022 Build Tools
(MSVC v143 + Windows SDK), CMake 3.27+, Ninja, [vcpkg](https://github.com/microsoft/vcpkg)
with `VCPKG_ROOT` set, and PowerShell. The CUDA presets additionally need the NVIDIA CUDA
Toolkit (the checked-in preset targets `sm_86`).

```powershell
.\build.ps1                    # fast dynamic development build + tests
.\build.ps1 -Release           # static, self-contained release build + package
.\build.ps1 -Science           # development build with the live LLM and web tools
.\build.ps1 -Release -Science  # static release with the optional science features
```

The development executable lands at `build\windows\bin\BioCAD.exe`, with the data packs copied
beside it. Static presets produce an executable that needs no Visual C++ redistributable.

```powershell
cmake --preset windows
cmake --build --preset windows
ctest --preset windows --output-on-failure
```

| Preset | Purpose |
| --- | --- |
| `windows` | Fast dynamic development build. |
| `windows-static` | Static-CRT, portable, single-file release build. |
| `windows-science` | Development build with curl-backed assistant and web tools. |
| `windows-science-static` | Static release build with the science features. |
| `windows-cuda` | CUDA docking build for supported NVIDIA environments. |
| `windows-cuda-static` | CUDA build with a statically linked CUDA runtime. |

### Headless docking self-test

```powershell
.\build\windows\bin\BioCAD.exe --selftest-dock --smiles "CC(N)Cc1ccccc1" --target DAT --compute cpu
```

`--compute` takes `auto`, `cpu` or `gpu`. Exit `0` means a real engine produced the poses
(`Provenance::Model`); exit `2` means the run fell back to the labelled descriptor estimate
(`Provenance::Heuristic`). First-run provisioning needs network access and stores components
under `%APPDATA%\BioCAD`.

### Deterministic screenshots

The app renders its own back buffer to PNG through WIC, so documentation media can be produced
on a CI runner with no interactive desktop session:

```powershell
.\build\windows\bin\BioCAD.exe --shot docs\media\docking.png --shot-panel Docking
.\scripts\capture-docs.ps1                 # every panel, plus the GIF frame sequences
```

`--shot-frames N` writes `<stem>-0000.png` onwards for an animation; `--shot-warmup N` sets
how many frames render before the first capture; `--shot-size W H` fixes the client area
(1600x1000 by default). Exit `3` means a capture was requested and did not fully succeed.

## Local data and runtime layout

```text
%APPDATA%\BioCAD\
├── biocad.db        # run history (SQLite, WAL)
├── artifacts\       # content-addressed structures, poses, and reports
├── runtime\         # provisioned engines and prepared receptors
├── packs\           # user data packs; override a built-in pack by id
├── presets\         # user preset overrides
├── logs\            # rotating logs and crash minidumps
├── cache\           # download and retrieval cache
└── manifest.json    # provisioned component versions and checksums
```

An existing `%APPDATA%\StimLab` tree from the previous name is migrated once on first launch
so multi-gigabyte provisioned engines are not downloaded again; cached receptors are dropped
because they carry the old box marker. API keys are encrypted at rest with Windows DPAPI.

## Repository guide

```text
src/
├── agent/      providers, tool registry, prompts, and web tools
├── app/        Win32 entry point, CLI, and the capture path
├── bio/        protein structure/sequence I/O, alignment, superposition, lDDT, SASA
├── chem/       molecular graph, descriptors, ADMET/PK model, 3D embedding
├── contracts/  frozen service and backend interfaces (the Services seam)
├── core/       paths, config, logging, manifests, secrets, generated version
├── data/       domain DTOs, Provenance and Quantity
├── modules/    real backend services and docking backends
├── numeric/    the one Levenberg-Marquardt fitter and the one RK4 integrator
├── packs/      versioned JSON data-pack schema and loader
├── pkpd/       dose-response fits, Cheng-Prusoff, Schild, PK/occupancy engine
├── render/     DirectX 11 device, WIC capture, molecular viewport
├── storage/    SQLite run store
├── ui/         application shell, panels, theme
└── workflow/   DAG execution and job scheduling

assets/packs/   built-in compound and target catalogs
docs/           architecture, provenance, cheminformatics, packs, docking, pkpd, protein, data sources, limits
tests/          Catch2 suite
scripts/        build, CI, packaging, signing, capture helpers
```

## Testing and CI

Catch2 through CTest. The GitHub Actions workflow builds and tests the dynamic and static
presets, asserts that `vcpkg.json` and `project(BioCAD VERSION ...)` agree, packages the static
release, captures every panel with `--shot`, and encodes the animation sequences to GIF on a
Linux runner.

```powershell
ctest --preset windows --output-on-failure
.\scripts\ci.ps1
```

CUDA docking must be checked on a real NVIDIA machine with `--selftest-dock --compute gpu`;
hosted runners have no GPU.

## What this deliberately does not do

Out of scope by design, and enforced in the assistant's system prompt and tool registry:

- No synthesis routes, reaction steps or conditions, precursor selection or acquisition,
  yields, equipment, or scale-up guidance. Binding affinity is a target-engagement signal,
  never a make-it signal.
- No dose, dose change, or personal regimen recommendation - the app emits exposure scenarios
  with their assumptions listed, and nothing else.
- No nanomolar affinity derived from a docking score.
- No invented binding sites: a target without a real structure is a coverage gap.

The full list, including the third-party tools deliberately not integrated and why, is in
[docs/limitations.md](docs/limitations.md).

## Contributing

1. Build and test with the PowerShell wrappers or the matching CMake preset.
2. Keep the UI thread non-blocking; preparation, provisioning and docking belong on workers.
3. Never present a fallback estimate as a real engine result, and never emit a derived number
   without a `Provenance`.
4. Add deterministic tests for behaviour changes and run the relevant CTest preset.
5. Keep credentials, downloaded engines, local databases and generated artifacts out of commits.

## License

Apache-2.0. See [LICENSE](LICENSE) and [NOTICE](NOTICE); third-party data sources and their
terms are listed in [docs/data-sources.md](docs/data-sources.md).
