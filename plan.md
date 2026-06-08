# StimLab — Native Single-Binary Rewrite — Master Plan (self-contained)

> **How to use this document.** This is a complete, standalone brief for rewriting StimLab into a
> single native Windows application. It assumes **no prior context** — a fresh agent session can
> execute from this file alone. The legacy Python/React implementation lives in `backend/` and
> `frontend/` of this repo and is the **behavior-parity reference** (read it to learn *what* a
> feature does; we are **not** copying *how* it's done). Start at **§7 Execution model → Phase 0**.

---

## 1. What StimLab is (functional surface to preserve)

StimLab is a computational **drug-discovery suite** — an open alternative to Schrödinger
Maestro/Glide/Desmond. It prepares molecules and proteins, docks ligands into receptors, screens
properties (ADMET), and orchestrates these as reproducible workflows, with a built-in AI assistant
that can drive the whole app. Research focus: CNS stimulants and their targets (monoamine
transporters DAT/NET/SERT, TAAR1, aminergic GPCRs) — expressed as **configurable presets**, never
hardcoded.

Features that must survive the rewrite (all exist today in `backend/stimlab/` + `frontend/src/`):

| Area | What it does | Legacy reference |
|---|---|---|
| **Ligand prep** | SMILES → validated 3D molecule, conformers, protonation, format conversion | `backend/stimlab/modules/ligand_prep/` |
| **Protein prep** | PDB load/clean/protonate, define docking box (per-target presets) | `backend/stimlab/modules/protein_prep/` |
| **Docking** | dock ligand(s) into a receptor box → ranked scored poses; swappable engines | `backend/stimlab/modules/docking/` |
| **ADMET / metabolism** | RDKit-descriptor + SMARTS heuristics + data-driven screening of 16 metabolic endpoints; optional ML overlay | `backend/stimlab/modules/admet/` |
| **Library** | store/import ligands (SMILES/SDF) + receptors (PDB), dedup, "dock selected" | `backend/stimlab/services/library.py` |
| **Workflows** | DAG templates capturing prep→dock knobs, re-runnable; live DAG view | `backend/stimlab/workflow/`, `services/workflows.py` |
| **Runs / Results** | run history, per-run poses, readable Markdown summary | `backend/stimlab/services/summary.py`, `api/runs.py` |
| **Archive** | save run summaries, side-by-side compare matrix, download .md | `backend/stimlab/services/archive.py` |
| **Presets** | 29 CNS target presets (PDB ref + box) + ADMET panels, as YAML data | `backend/stimlab/presets/data/` |
| **AI agent** | multi-provider (Anthropic / OpenAI-compatible incl. DeepSeek / Ollama), tool-calling loop, autopilot + ask-first, ~10 tools, drives the UI (`navigate_ui`) | `backend/stimlab/agent/` |
| **Settings** | provider/key config, GPU mode, paths | `frontend/src/pages/SettingsPage.tsx` |

MD/FEP/QM/viz are **typed stubs** today and stay stubs in v1 (OpenMM-CUDA is a later phase).

---

## 2. Goal

Rewrite StimLab as a single self-provisioning native **x64 Windows** application:

- **C++ core** (no Python, no Docker, no localhost server) with a **Dear ImGui** desktop UI on
  **DirectX 11**, an in-app **3D molecular viewer**, **CUDA** hardware acceleration where it counts,
  and an expanded multi-provider **AI agent** with **web tools**.
- **One self-provisioning `.exe`**: first launch downloads/extracts its runtime (CUDA runtime libs,
  engine binaries, ML models) into `%APPDATA%/StimLab`. Everything — DB, artifacts, runtime, logs,
  presets — lives under `%APPDATA%/StimLab`.

### Locked decisions
- **Packaging:** one self-provisioning `.exe` (no installer required; self-heals its runtime).
- **Science:** bundle/compile proven native engines (Vina, AutoDock-GPU, RDKit, ONNX Runtime) —
  **do not reinvent the algorithms**.
- **Web tools:** keyless (no API key) — see §4 for the durable approach.
- **Render API:** DirectX 11.
- **Language:** C/C++ (C++20). This is settled; do not relitigate.
- Licenses are irrelevant (personal use) — static-link/bundle GPL/LGPL freely.

---

## 3. Architecture

Single process, three layers, all in-process:

```
 main thread │ UI — Dear ImGui (docking) on DirectX 11 + Win32
            │ shell/nav · feature panels · 3D mol viewport · assistant
            │      ▲ commands down / events up  (lock-free queues — UI NEVER blocks)
 worker pool │ Core/services — DAG workflow engine · job/thread pool ·
            │ storage (SQLite + artifact store) · modules · agent
            │      ▲ links / invokes
 native libs │ RDKit (linked) · ONNX Runtime (linked, prebuilt) ·
            │ Vina / AutoDock-GPU / gnina / obabel (subprocess binaries) ·
            │ libcurl · WebView2 · SQLite · yaml-cpp · nlohmann-json
```

**Hard rule:** the UI thread never blocks on science. Panels submit a **command** (e.g.
`SubmitDocking`) to the job system and render from the latest **snapshot/event**. Long work runs on
the worker pool; progress flows back as events. Make Phase-0 stubs **thick fakes** (return plausible
data) so UI work can proceed against something that actually renders.

---

## 4. Tech stack

### 4a. Critical principle — *link* vs *provision* (do not build what ships as a binary)

The #1 way to waste time/money here is building heavy libraries from source. Split deps:

**Linked in-process (the only hard linkage worth doing):**
- **RDKit** (C++) — molecule parsing (SMILES/PDB/SDF/Mol2), 3D embedding/conformers, descriptors,
  SMARTS (used by the viewer, ligand prep, and ADMET in the hot path). This is the gating linkage —
  prove it works first (see Phase 0.5).
- **ONNX Runtime** — **download Microsoft's official prebuilt Windows x64 GPU build** (CUDA EP) and
  link the prebuilt; **never vcpkg-build it from source**.

**Provisioned as prebuilt binaries, invoked as subprocesses** (matches how the legacy app worked):
- **AutoDock Vina 1.2** — ships official Windows binaries; default docking engine, subprocess.
- **Open Babel** — use **`obabel.exe` as a subprocess** for PDBQT/format conversion (linking
  libopenbabel on Windows is the shakiest dependency; avoid it).

**Real from-source build effort (no clean official Windows binaries — flag as the hard work):**
- **AutoDock-GPU (CUDA)** — GPU fast-screening engine; this is where "hardware-accelerated docking"
  actually comes from. Real Windows+CUDA build effort.
- **gnina (CUDA CNN)** — **stretch / post-v1.** No official Windows build. Fallback: export its CNN
  to ONNX and grid via libmolgrid, score through ONNX Runtime. **v1 ships without gnina** if the
  build doesn't converge — Vina + AutoDock-GPU already give real CPU + GPU docking.

**Dropped from the legacy stack:** `smina` (a pre-1.2 Vina fork; painful on Windows, ~redundant with
Vina 1.2 scoring). Vina 1.2 + AutoDock-GPU cover CPU + GPU.

### 4b. Everything else
- **Build:** MSVC v143 (VS 2022 Build Tools) + Windows 10/11 SDK · CMake ≥3.27 + Ninja · **vcpkg
  manifest mode** · CUDA Toolkit 12.4+. **Static-link the CRT (`/MT`)** so the bootstrapper never
  needs to install the VC++ redist. **Stand up a shared vcpkg binary cache** before fan-out
  (RDKit+Boost is an hours-long build — build once, restore everywhere).
- **UI:** Dear ImGui (docking branch) + ImPlot (charts) + ImGuizmo (viewer gizmos) + DX11/Win32 backends.
- **3D viewer:** custom DX11 renderer, GPU-instanced **impostor** spheres/cylinders (ball-and-stick);
  reference existing impostor shaders (Speck/3Dmol/PyMOL-style) rather than inventing them.
  Ribbon/cartoon/surface are a later pass.
- **Storage:** SQLite (amalgamation, C API, WAL mode) + content-addressed artifact store on disk.
- **HTTP/LLM:** libcurl (Schannel) with SSE streaming.
- **JSON/YAML:** nlohmann/json + yaml-cpp (reuse the existing preset YAML files as data).
- **DAG/jobs:** consider **Taskflow** (header-only, vcpkg) instead of hand-rolling a DAG executor +
  thread pool — it is almost exactly this problem.
- **GPU dashboard:** **NVML** (`nvml.dll`, ships with the NVIDIA driver) for device/VRAM/utilization.
- **Web tools (keyless, durable):** use **DuckDuckGo's HTML endpoint** or a **SearXNG** instance for
  the *search* step (scraping Google/Bing SERPs trips captcha/anti-bot fast, even via WebView2);
  reserve **WebView2 headless** for `web_fetch` of JS-heavy pages + HTML→text extraction.
- **Secrets:** API keys encrypted at rest via **Windows DPAPI**, never plaintext in the DB.
- **Hashing:** xxhash / BLAKE3 for cache keys + content addressing.
- **Logging/crash:** spdlog (rotating) + `MiniDumpWriteDump` on unhandled exceptions (native CUDA
  will crash — minidumps turn 2-hour diagnoses into 5-minute ones).
- **Test:** Catch2/doctest.

### Where hardware acceleration shows up
Docking (AutoDock-GPU CUDA; gnina CUDA stretch) · ADMET ML (ONNX Runtime CUDA EP) · rendering
(DX11 GPU instancing) · future MD/FEP (OpenMM CUDA). Detect NVIDIA GPU + driver at startup (NVML);
gracefully fall back to CPU (Vina, ONNX CPU EP) when absent. CPU/GPU is a setting, **auto** by default.

---

## 5. `%APPDATA%/StimLab` layout (all state here; nothing in Program Files)

```
%APPDATA%/StimLab/
  stimlab.db          SQLite: projects, ligands, receptors, runs, poses, library,
                      workflows, archive, settings, agent config (keys as DPAPI blob)
  artifacts/          content-addressed files: prepped structures, poses, reports (ab/cdef…)
  runtime/
    cuda/             provisioned CUDA 12 runtime DLLs (only if the system lacks them)
    engines/          vina.exe, autodock_gpu*, (gnina*), obabel.exe + data
    models/           ONNX ADMET models, (gnina CNN weights)
  presets/            shipped YAML presets + user overrides (29 CNS targets, ADMET panels)
  logs/               rotating logs + minidumps
  cache/              download/temp/scrape cache
  manifest.json       provisioned-component versions + checksums (self-heal source of truth)
```

---

## 6. Repository layout (new native tree)

```
/CMakeLists.txt            superbuild root
/CMakePresets.json        windows-cuda preset
/vcpkg.json               dependency manifest
/cmake/                   toolchain, CUDA, packaging helpers
/docs/
  ARCHITECTURE.md         condensed design
  contracts/             ONE markdown per frozen interface (shared spec, §8)
  workpackages/          ONE self-contained brief per WP (hand each to an agent)
/src/
  core/      paths, logging, errors, hashing, event bus, config, DPAPI secrets
  data/      domain types (Molecule, Ligand, Receptor, Pose, Run, Artifact, Provenance)
  storage/   sqlite wrapper, DAOs, artifact store, migrations
  workflow/  DAG engine, node/registry, job pool, progress
  modules/   module base + ligand_prep, protein_prep, docking/*, admet/*  (md/fep/qm/viz stubs)
  agent/     provider layer (anthropic/openai/ollama), tool registry, loop, web tools
  render/    dx11 device, imgui bootstrap, mol viewport renderer
  ui/        app shell, nav, panels (one file per feature page), assistant panel
  app/       WinMain, bootstrapper/provisioner, lifecycle
/third_party/  pinned engine build recipes (vina fetch, autodock-gpu, gnina)
/tests/        catch2 unit + golden tests (ported from legacy pytest assertions)
/packaging/    self-extracting payload assembly + release script
```

---

## 7. Execution model (built for many limited-context agents in parallel)

The strategy: **freeze interfaces once, prove one slice end-to-end, then fan out.** Each agent gets
only its WP brief + the 1-3 contract docs it touches + the `data/` header — never the whole repo.

### Phase 0 — Foundation / Contract Pack *(serialized; strong model; blocks everything)*
Produces a **compiling skeleton** where every interface exists as a header with a **thick fake** impl.
Deliverables:
1. Repo skeleton + superbuild `CMakeLists.txt` + `vcpkg.json` + `CMakePresets.json` + the **vcpkg
   binary cache** set up. Mark `rdkit` as the risky port; ONNX Runtime is a prebuilt download, not a port.
2. `src/core`: `AppPaths` (resolve/create `%APPDATA%/StimLab/*`), logging + minidump, `Error`/`Result<T>`,
   hashing, thread-safe typed `EventBus`, `Config`, DPAPI secret store.
3. `src/data`: domain structs (JSON-serializable) — the lingua franca.
4. **Frozen interface headers**, each with a `docs/contracts/*.md`: `IModule` + `ModuleContext`,
   `IDockingBackend`, `IStorage`/DAOs, `IArtifactStore`, `IJobSystem`, `ILlmProvider` + `ITool` +
   `IToolRegistry`, `IWebTools`, `IMolRenderer` + `IViewport`, and the `UiCommand`/`UiEvent` enums.
5. Thick fakes for every interface + a smoke-test `WinMain` (ImGui+DX11 window) + a Catch2 test.
   **Done-when:** `cmake --build` → `StimLab.exe` opens a window; `ctest` green.
6. Write `docs/ARCHITECTURE.md` + all `docs/workpackages/WP-*.md` briefs.

### Phase 0.5 — Walking skeleton *(serialized; strong model; do BEFORE fan-out)*
Build **one thin vertical slice on a single thread of work** that exercises every scary integration
point at once, proving the contracts against reality:
**paste SMILES → RDKit prep (linked) → Vina dock (subprocess) → poses in the 3D viewport → row in
SQLite + artifact store.** If RDKit won't link or the viewport won't draw, you find out *here*, not
across 8 parallel agents. This slice becomes the reference every later agent copies.

### Phase 1 — Parallel work packages (§8)
Fan out only after 0 + 0.5 are green. **Tier by difficulty:** strong model on the gnarly WPs
(CUDA builds, viewer shaders, agent loop); cheaper agents on the well-bounded mechanical WPs (UI
panels against fakes, DAO CRUD, tool-rebinding, golden tests).

### Phase 2 — Integration & packaging
A dedicated **integration owner** keeps a single trunk continuously green, swapping fakes for real
impls as they land. Every WP must ship a runnable `--selftest`/demo — "compiles" ≠ "integrates."
Then packaging (§8 WP-M).

---

## 8. Work packages (each is a self-contained agent brief)

Per-WP format: *Objective · Context budget (exact files to read) · Implements/consumes · Acceptance.*
Tier = suggested model strength.

### Wave 1 — depend only on the contract pack
- **WP-A · App shell + DX11/ImGui platform** *(strong).* Win32 window, DX11 device/swapchain, ImGui
  docking bootstrap, DPI/fonts/theme, dockspace, nav rail, persisted layout, the command/event pump
  bridging UI↔core. *Accept:* dockable empty panels per feature page, theme toggle, layout persists.
- **WP-B · Self-provisioning bootstrapper** *(strong).* Read `manifest.json`; ensure WebView2 runtime
  + (if no system CUDA) CUDA 12 runtime DLLs + engine binaries + ONNX models present & checksum-valid;
  download (libcurl) + extract (libarchive) into `runtime/`; **self-heal**: verify *completeness*, else
  wipe + reinstall (see §9). Progress via events. *Accept:* clean box provisions fully; a corrupted
  component re-provisions next launch.
- **WP-C · Storage layer** *(cheap-ok).* SQLite schema + migrations + DAOs for all entities;
  content-addressed `IArtifactStore`; settings + agent config (keys via DPAPI). *Accept:* CRUD
  round-trip tests; artifact dedup by hash; migration test.
- **WP-D · Workflow DAG + job system** *(medium).* Thread-pool `IJobSystem` (consider Taskflow); DAG
  with `cache_key = hash(module,version,inputs)`, content-cached resumable nodes; progress events;
  cancel tokens. *Accept:* a 3-node toy DAG runs, caches, resumes, cancels; concurrency stress test.
- **WP-N · Test harness + golden data** *(cheap-ok).* Catch2 scaffolding + golden assertions ported
  from the legacy pytest suite (`backend/tests/`). Gives every WP an acceptance gate.

### Wave 2 — depend on Wave-1 pieces
- **WP-F · Protein prep** *(medium; needs RDKit from 0.5).* PDB load/clean/protonate/add-H + box
  definition; reuse the 29 CNS target preset YAMLs (`backend/stimlab/presets/data/`) via yaml-cpp.
  *Accept:* prep a known transporter PDB → cleaned receptor + valid box.
- **WP-G · Docking engines** *(strong; split G1/G2/G3).*
  - **G1** Vina 1.2 subprocess `IDockingBackend` (default); SDF→PDBQT via `obabel.exe`.
  - **G2** AutoDock-GPU built with CUDA, subprocess backend (GPU fast screening).
  - **G3** gnina CNN — **stretch**, ONNX-export fallback; v1 may ship it disabled.
  *Accept (G1/G2):* dock known ligand/receptor → finite scores + ranked poses; GPU path runs on the
  NVIDIA box; scores within tolerance of reference.
- **WP-H · ADMET / metabolism** *(medium; needs RDKit).* Port the RDKit-descriptor + SMARTS heuristics
  and the data-driven `metabolism.yml` (16 endpoints) from `backend/stimlab/modules/admet/`; optional
  ADMET-AI overlay via ONNX Runtime. *Accept:* golden verdicts in §9.
- **WP-I · 3D molecular viewport** *(strong; needs RDKit + WP-A).* DX11 instanced impostor
  spheres/cylinders, element coloring, pick, camera (orbit/pan/zoom), **pose overlay**, receptor+ligand
  co-display, ImGui-docked. *Accept:* load receptor PDB + poses, cycle poses, pick atoms.
- **WP-J · Agent core** *(strong).* `ILlmProvider`: Anthropic, OpenAI-compatible (incl. DeepSeek),
  Ollama — libcurl + SSE streaming. Tool-calling loop, autopilot + ask-first, tool registry; rebind
  the legacy ~10 tools (list_backends, run docking, list/inspect runs+poses, summarize_run,
  archive_run, screen_metabolism, navigate_ui, presets) onto native services. *Accept:* mock-provider
  multi-tool turn test; live smoke vs one real key.
- **WP-K · Agent web tools** *(medium; needs WP-J).* `web_search` (DuckDuckGo HTML / SearXNG →
  structured hits), `web_fetch` (WebView2 headless for JS pages → HTML→text), `compare_results`
  synthesis tool; JSON schemas, rate-limit + cache in `cache/`. *Accept:* query returns parsed hits;
  fetched article yields clean text; agent calls them in a turn.
- **WP-L · Feature UI pages** *(cheap-ok against fakes; split L1/L2/L3).* ImGui panels backed by
  services: **L1** Dashboard + Docking + Runs/Results (PosesTable, RunSummaryCard, viewer embed);
  **L2** Library + Workflows (live DAG) + Archive (compare matrix); **L3** Metabolism + Presets/Targets
  + Settings (providers/keys, GPU mode, paths) + Assistant panel (streamed tool calls, autopilot/
  ask-first toggle, navigate_ui targets). *Accept:* each page drives its service end-to-end;
  navigate_ui focuses the right panel.

### Wave 3 — integration & packaging
- **WP-M · Packaging / release** *(strong).* Assemble the self-extracting payload + final
  `StimLab.exe`; embed base `manifest.json`; CI build on a Windows+CUDA runner; integration pass +
  perf check. *Accept:* fresh-VM run → app provisions and reaches the dashboard with no Docker/Python;
  UI holds ~60fps while a GPU dock runs.

---

## 9. Hard-won lessons to preserve (a cold session won't know these)

**Behavior-parity golden facts (use as test oracles):**
- Default ligand format is **SDF** for all engines. **Do not feed Meeko-style PDBQT to Vina** — it
  trips a `tree.h` error; convert SDF→PDBQT via `obabel`.
- ADMET verdicts: **amphetamine / methamphetamine / MDMA → WARN** (MAO + CYP2D6); **dopamine → WARN**
  (+ COMT, catechol); **acetaminophen → flags a reactive metabolite (NAPQI)**; **caffeine → INFO**.
- 29 CNS target presets exist (transporters, dopamine/serotonin receptors, TAAR1, adrenergic,
  histamine, muscarinic, opioid, sigma1, CB1, NMDA, enzymes MAO-A/B/COMT/AChE/CYP2D6, hERG safety),
  each with a real PDB ref + a box derived from the co-crystal ligand.
- ADMET heuristics are **not ML**: phenethylamine→MAO, catechol→COMT, ester/amide→hydrolysis,
  basic amine near arene→CYP2D6, structural alerts→bioactivation. Degrade gracefully if RDKit absent.

**Provisioning / CUDA self-heal (port this exactly into WP-B):**
- A provisioned runtime dir's *existence* is **not** proof of completeness. Track the required set of
  files/sonames; on launch verify **completeness**, and if incomplete **wipe + reinstall** rather than
  trusting the dir. (Legacy hit "libcufft.so.11 missing" because a stale CUDA cache was trusted.)
- gnina's CUDA needs are specific (libcudart, libcublas/Lt, libcufft, libnvToolsExt, libcusolver,
  libcusparse, libnvJitLink; cuDNN only for the GPU path). On Windows the equivalent CUDA 12 runtime
  DLLs must be provisioned the same way if gnina/AutoDock-GPU are dynamically linked.

**Windows text gotchas (relevant to any scripts):** keep launcher/CLI text **ASCII** — the Windows
console (cp1252) and PowerShell 5.1 (Windows-1252) mangle em-dashes and smart quotes.

---

## 10. Risk register (front-loaded)
- **RDKit linkage (gating).** Half the WPs depend on it. Prove it in Phase 0.5 before fan-out; vcpkg
  port is a long build → use the binary cache; allow a prebuilt-binary fallback.
- **AutoDock-GPU / gnina Windows+CUDA builds.** No clean official binaries. AutoDock-GPU is the real
  GPU-docking effort; gnina is stretch with an ONNX fallback and may ship disabled in v1.
- **ONNX Runtime GPU.** Use Microsoft's official prebuilt Windows GPU build; never build from source.
- **Keyless web search brittleness.** Prefer DuckDuckGo HTML / SearXNG over scraping Google/Bing;
  cache + throttle + degrade gracefully to "no results."
- **Parallel agents that compile but don't integrate.** Thick fakes + per-WP `--selftest` + an
  integration owner on a continuously-green trunk.
- **Unsigned self-downloading exe** trips SmartScreen/AV. Acceptable for personal use — expected.
- **3D viewer scope creep.** Ship ball-and-stick + poses in v1; ribbon/surface later.

---

## 11. Prerequisites for building agents
VS 2022 Build Tools (MSVC v143) + Windows 10/11 SDK · CMake ≥3.27 + Ninja · CUDA Toolkit 12.4+ ·
vcpkg (bootstrapped, with a shared binary cache) · Git · an NVIDIA-GPU machine for GPU-path
acceptance tests.

---

## 12. Verification (end-to-end, on the NVIDIA host)
1. **Build:** `cmake --preset windows-cuda && cmake --build --preset windows-cuda` → `StimLab.exe`; `ctest` green.
2. **Cold start:** clean profile → bootstrapper provisions `%APPDATA%/StimLab`; dashboard opens; no Docker/Python.
3. **Science slice:** import a stimulant SMILES → ligand prep → pick a CNS target preset → protein
   prep → dock (Vina CPU and AutoDock-GPU) → poses render in the viewport → archive + downloadable
   summary. Scores finite and within tolerance of golden references.
4. **ADMET:** screen amphetamine/MDMA/acetaminophen → verdicts match §9 goldens.
5. **Agent:** connect a provider → autopilot runs a docking job + a `web_search`/`web_fetch` and
   reports back; `navigate_ui` focuses panels.
6. **Perf:** UI holds ~60fps during a GPU dock; GPU dashboard (NVML) shows utilization.
7. **Self-heal:** corrupt a `runtime/` component → next launch re-provisions it.

---

## 13. Open items to confirm early
- Is **gnina CNN** required for v1, or acceptable as post-v1 (Vina + AutoDock-GPU ship real CPU+GPU docking now)?
- Keep **MD/FEP/QM/viz** as honest typed stubs in v1 (OpenMM-CUDA later)? (assumed yes)
- App identity/branding (window title, icon) — reuse "StimLab" or change?
