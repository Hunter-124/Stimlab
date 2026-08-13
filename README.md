# BioCAD

**A native Windows workstation for computational drug discovery.**

BioCAD brings ligand and receptor preparation, molecular analysis, docking, ADMET screening, reusable workflows, and an in-app assistant into one C++20 desktop application. It is designed around CNS-research use cases while keeping targets and presets configurable rather than hard-coded.

> **Research software notice**
> BioCAD is intended for lawful research and educational use. Computational predictions are hypotheses—not clinical, safety, efficacy, or regulatory conclusions—and must be validated by qualified professionals and appropriate experimental methods.

## Why BioCAD

- **Native and self-contained** — a Win32/DX11 desktop application with SQLite-backed local state; release builds can be packaged as a single portable executable.
- **Real docking paths** — provisions and invokes AutoDock Vina for flexible CPU docking, with Vina-GPU/OpenCL and a native CUDA rigid-body backend available for GPU workflows.
- **Practical molecular tooling** — SMILES parsing, 3D embedding, descriptors, similarity, structural alerts, property screening, and a 3D molecular viewport.
- **Reproducible workflows** — a cancellable, content-cached DAG engine models preparation and docking pipelines without blocking the UI.
- **Local-first data** — projects, artifacts, logs, engine runtimes, and configuration live under `%APPDATA%\BioCAD`.
- **Optional assistant and web tools** — multi-provider LLM support, encrypted API-key storage through Windows DPAPI, and optional web retrieval in science-enabled builds.

## Capabilities

| Area | What it provides |
| --- | --- |
| **Ligand preparation** | SMILES validation, molecular graph construction, 3D conformer embedding, descriptors, and format-oriented utilities. |
| **Receptor preparation** | PDB cleanup, protonation-oriented preparation, docking boxes, and configurable CNS-target presets. |
| **Docking** | Ranked poses through AutoDock Vina; optional Vina-GPU/OpenCL and first-party CUDA execution paths, with explicit fallback labeling. |
| **ADMET screening** | Descriptor- and rule-based property and metabolism screening, with optional science-stack extensions. |
| **Workflows** | Re-runnable, content-cached preparation-to-docking pipelines with progress, cancellation, and a live DAG view. |
| **Results and storage** | SQLite run history, artifacts, summaries, comparison views, and local project data. |
| **Assistant** | Tool-driven in-app assistance for navigating and operating supported application workflows; enabled only when configured. |

## Architecture

BioCAD keeps interactive work responsive by separating presentation, services, and native engines in a single process:

```text
┌───────────────────────────────────────────────────────────────────────────┐
│ UI thread: Win32 + Dear ImGui + ImPlot + DirectX 11                       │
│ navigation · panels · molecular viewport · assistant controls             │
└──────────────────────────────┬────────────────────────────────────────────┘
                               │ commands / snapshots / events
┌──────────────────────────────▼────────────────────────────────────────────┐
│ Worker services: workflow DAG · job system · docking · chemistry · agent  │
│ SQLite run store · artifact storage · provisioning                         │
└──────────────────────────────┬────────────────────────────────────────────┘
                               │ native libraries and subprocess engines
┌──────────────────────────────▼────────────────────────────────────────────┐
│ Vina · Vina-GPU · CUDA backend · SQLite · YAML · JSON · optional curl     │
└───────────────────────────────────────────────────────────────────────────┘
```

Long-running preparation, provisioning, and docking jobs run off the UI thread. The application surfaces an engine result as real only when a real engine completed the work; unavailable engines degrade to a clearly labeled estimate rather than silently claiming a docking result.

## Quick start

### Prerequisites

BioCAD is built and supported on **Windows x64**. Install:

- Visual Studio 2022 Build Tools with the MSVC v143 C++ toolchain and Windows SDK
- CMake 3.27 or newer and Ninja
- [vcpkg](https://github.com/microsoft/vcpkg), with `VCPKG_ROOT` set to its installation directory
- PowerShell 7 or Windows PowerShell
- NVIDIA CUDA Toolkit only for the CUDA presets (the checked-in CUDA preset targets `sm_86`; adjust it for another GPU)

The repository uses vcpkg manifest mode. Configure and build through the supplied PowerShell wrappers so that the compiler environment and dependency root are prepared consistently.

```powershell
# Fast development build and test suite
.\build.ps1

# Static, self-contained release build, tests, and release package
.\build.ps1 -Release

# Development build with live LLM and web-tool support
.\build.ps1 -Science

# Static release build with the optional science features
.\build.ps1 -Release -Science
```

The development executable is produced at `build\windows\bin\BioCAD.exe`. Static release presets produce an executable that does not require the Visual C++ redistributable beside it.

### CMake presets

Use the presets directly when integrating with your own Windows build environment:

```powershell
cmake --preset windows
cmake --build --preset windows
ctest --preset windows --output-on-failure
```

| Preset | Purpose |
| --- | --- |
| `windows` | Fast dynamic development build with the core application. |
| `windows-static` | Static-CRT, portable release-oriented build. |
| `windows-science` | Development build with curl-backed assistant and web-tool features. |
| `windows-science-static` | Static release build with the optional science features. |
| `windows-cuda` | Lean CUDA docking build for supported NVIDIA environments. |
| `windows-cuda-static` | CUDA build with a statically linked CUDA runtime for portable GPU execution. |

## Running a docking self-test

The executable exposes a headless acceptance path that provisions the required assets and runs one docking job:

```powershell
.\build\windows\bin\BioCAD.exe --selftest-dock --smiles "CC(N)Cc1ccccc1" --target DAT --compute cpu
```

`--compute` accepts `auto`, `cpu`, or `gpu`. A zero exit code means a real docking engine completed the run; exit code `2` means BioCAD returned its explicitly labeled descriptor estimate instead. First-run provisioning requires network access and stores downloaded runtime components under `%APPDATA%\BioCAD`.

## Local data and runtime layout

BioCAD keeps mutable state outside the installation directory:

```text
%APPDATA%\BioCAD\
├── biocad.db       # projects, molecules, receptors, runs, settings
├── artifacts\       # content-addressed structures, poses, and reports
├── runtime\         # provisioned engines, runtimes, and models
├── presets\         # shipped presets plus user overrides
├── logs\            # rotating logs and crash diagnostics
├── cache\           # download and retrieval cache
└── manifest.json    # provisioned component versions and checksums
```

API keys configured for supported assistant providers are protected with Windows DPAPI. Do not commit keys, generated artifacts, or local runtime data to source control.

## Repository guide

```text
src/
├── app/        Win32 application entry point
├── agent/      providers, tool registry, prompts, and web tools
├── chem/       molecular graph, descriptors, analysis, and 3D embedding
├── contracts/  stable service and backend interfaces
├── core/       configuration, paths, logging, manifests, and secrets
├── data/       domain models and preset-oriented data
├── modules/    application services and docking backends
├── render/     DirectX 11 device and molecular viewport
├── storage/    SQLite run and artifact storage
├── ui/         application shell, panels, and theme
└── workflow/   DAG execution and job scheduling

tests/          Catch2 unit and integration-style coverage
scripts/        Windows build, CI, packaging, signing, and capture helpers
```

For implementation details, browse the source tree and automated tests.

## Testing and CI

The test suite is built with Catch2 and registered through CTest. The GitHub Actions workflow validates the dynamic and static Windows presets, runs CTest, and packages the static release artifact.

```powershell
# Run a configured preset's tests
ctest --preset windows --output-on-failure

# Exercise the local CI workflow
.\scripts\ci.ps1
```

CUDA-capable docking should additionally be checked on a compatible NVIDIA machine with the `--selftest-dock --compute gpu` command above. GitHub-hosted runners do not provide a GPU.

## Contribution guidelines

1. Build and test with the supplied PowerShell scripts or the matching CMake preset.
2. Keep UI work non-blocking; move preparation, provisioning, and docking work onto worker services.
3. Preserve explicit result provenance—never present a fallback estimate as a real engine result.
4. Add deterministic tests for behavior changes and run the relevant CTest preset before opening a change.
5. Keep credentials, downloaded engines, local databases, and generated artifacts out of commits.

## License

No license file is currently included in this repository. Until one is added, do not assume permission to redistribute, modify, or use the project beyond the rights granted by its copyright holders.
