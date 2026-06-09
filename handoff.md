# StimLab — Session Handoff

> **Purpose of this file.** Bring a fresh session (human or agent) fully up to speed: what we're
> building, the hard rules, exactly what exists on disk right now, the verified environment, and the
> precise point to resume from. Read this first, then [plan.md](plan.md) (the master spec) and the
> approved build plan at `C:\Users\nigga\.claude\plans\i-want-you-to-deep-plum.md`.

Last updated: 2026-06-09 · Branch: `master` · **v1 FEATURE-COMPLETE** — all four user priorities landed.

> **STATUS: v1 COMPLETE. ctest 68/68 (windows + windows-static), science green too; all four priorities done.**
> The user's priority list is fully delivered:
> 1. **Real docking** — AutoDock Vina runs for real (WP-F receptor prep + box-from-co-crystal-ligand +
>    size-checked provisioning). Verified live: amphetamine→DAT −4.83, MDMA→SERT −4.26 kcal/mol.
> 2. **Workflows / DAG engine** (`src/workflow/`) — content-cached, resumable, cancellable prep→dock DAG
>    + live Workflows panel; off-thread. A review subagent caught + we fixed a critical notify-then-destroy UAF.
> 3. **Expansive agent tools (WP-K)** — service-action tools (analyze/screen_admet/dock/run_workflow/
>    list_runs/search_library/compare_compounds) + keyless web tools (web_search via DuckDuckGo HTML,
>    web_fetch) cached/rate-limited under `cache/web`, behind the science feature. Live web_search verified.
> 4. **Packaging (Phase E)** — `manifest.json` + self-heal (verify/heal corrupt components on launch) and a
>    **fully-static single-exe release** (`windows-static` preset, /MT, x64-windows-static): a 3.3 MB
>    StimLab.exe needing NO DLLs / NO VC++ redist. `scripts/package.ps1` zips it. Verified portable: the
>    exe run from a clean `%TEMP%` provisions + real-docks. Real docking works even in this curl-free build
>    (provisioning uses PowerShell); only the live LLM + web tools need the science build.
> See [docs/PHASE_E.md](docs/PHASE_E.md). Open follow-ups (post-v1): AutoDock-GPU/gnina (WP-G2/G3),
> WebView2 for JS-heavy `web_fetch`, prepare the other 25 receptors on demand, a science-static release.

> **STATUS: WORKFLOW / DAG ENGINE LANDED (WP-D, priority #2). ctest 61/61.** The re-runnable prep→dock
> pipelines are rebuilt as a content-cached, cancellable DAG with a live view:
> - **Engine** `src/workflow/` (hand-rolled, NO Taskflow dep): `JobSystem` (worker pool + cooperative
>   `CancelToken`), `Dag` (validate/topo), `DagExecutor` (parallel ready-queue scheduler). Content cache
>   key = `hash(module,version,params,{dep→dep key})` (transitive → unchanged graph = full cache-hit
>   resume; change one node → only it + downstream re-run). `INodeCache` with in-memory + `%APPDATA%/
>   StimLab/cache` disk backends.
> - **Pipeline** `src/modules/Pipelines.cpp`: `buildDockingPipeline(smiles,target,Services)` → 3 nodes
>   (ligand_prep → embed3D; receptor_prep → locate prepared PDBQT; dock → dockDetailed). External inputs
>   (receptor-prepared? vina present?) are fingerprinted into node params at build time so a provisioning
>   change correctly invalidates the cache.
> - **Live DAG view** = new **Workflows** panel: bezier node graph colored by per-node status
>   (pending/running/cached/done/failed/cancelled) + Run/Cancel + a node-output table. The run executes on
>   a worker thread via `AppShell::runWorkflow` (UI never blocks; mirrors the docking `dockFor` pattern).
> - An **adversarial concurrency review** (subagent) of the engine found a real critical use-after-free
>   (notify-then-destroy of the scheduler's `m`/`cv` by the last worker); fixed by joining worker futures
>   before `run()` returns. A throwing node fn is caught → Failed (no hang). Both fixes have tests.
> - Next per the user's priority order: **#3 expansive agent tools / web tools** (plan §8 WP-K). Packaging
>   (#4) stays deferred. See [docs/PHASE_D.md](docs/PHASE_D.md) for the docking detail.

> **STATUS: REAL DOCKING IS LIVE (WP-F + WP-G1), verified end-to-end. ctest 48/48.** Docking no longer
> falls back to the descriptor estimate when an engine is provisioned - it runs **AutoDock Vina for real**
> and renders ranked `real=true` poses. Landed this session:
> - **Receptor prep (WP-F)** `src/modules/docking/ReceptorPrep.*`: in-house **PDB → rigid receptor PDBQT**
>   (no RDKit / no hard OpenBabel link) - strips waters + additives + co-crystal ligand, AutoDock4 atom
>   types, heavy-atom (Vina scoring is charge-independent). Fetches the preset's PDB from **RCSB**, caches
>   `runtime/receptors/<id>.pdbqt`. Prefers `obabel.exe` (adds polar H) when it is dropped in.
> - **Box from the structure's own frame** (`receptorBoxFromPdb`): the docking box center comes from the
>   **co-crystal ligand centroid** in the fetched PDB's real coordinates (written as `REMARK STIMLAB_BOX`,
>   overrides the preset's literature center which is in a different frame). This is the load-bearing fix -
>   without it the box misses the receptor and Vina returns nothing.
> - **Provisioning** (`EngineLocator`): `ensureVina(true)` downloads `vina_1.2.5_win.exe` with a **size
>   check** (published 1,203,712 B; GitHub has no digest) + an **optional runtime SHA-256 pin** (env
>   `STIMLAB_VINA_SHA256` or `runtime/engines/vina.sha256`, no rebuild). `ensureObabel()` is locate-only.
> - **Off-thread** `docking::Provisioner` + `AppShell::dockFor()` (worker-thread, cached by molecule+target)
>   so the **UI never blocks**; Vina gets `--seed 1` for reproducibility. Docking panel has a "Provision
>   engine + receptors" button + live status; absent-engine path still degrades to the labeled estimate.
> - **Acceptance hook** `StimLab.exe --selftest-dock [--smiles S] [--target T]` (exit 0 = real dock).
>   **Verified live:** amphetamine→DAT = **−4.83 kcal/mol, 9 poses, real=true** (reproducible across runs);
>   MDMA→SERT = **−4.26, 9 poses, real**. GUI screenshot confirms the same in the Docking panel + 3D viewer.
> - Next per the user's priority order: **#2 Workflows / DAG engine** (plan §8 WP-D; Taskflow), then
>   **#3 expansive agent tools / web tools** (WP-K). Packaging stays deferred (do NOT start it).

> **STATUS: Phase B COMPLETE + ENRICHED, verified.** `./scripts/dev-build.ps1 -Test` (or
> `cmake --preset windows && cmake --build --preset windows`) produces `build/windows/bin/StimLab.exe`,
> which opens a clean three-pane DX11/ImGui GUI (Navigator | Workspace | Assistant) in **Segoe UI**.
> `ctest --preset windows` is green (**12/12**, incl. the §9 ADMET goldens). **14 feature panels** render
> against thick fakes over a **31-compound default library**; the assistant can focus + pulse-highlight
> any panel. Added this session:
> - **Absorption / PK** module (HIA, oral F%, Caco-2 permeability, logBB/CNS penetration, P-gp efflux).
> - **Analog Explorer** panel: tune a candidate derivative (property sliders + functional-group toggles)
>   and see live stability/absorption/ADMET, nearest existing sample, legal-analog score, and predicted
>   byproducts/interactions - analysis only, no synthesis.
> - **Compare** panel: grouped-bar + matrix comparison of up to three compounds.
> - Library expanded 16 -> 31 (MDA, methylone, MDPV, alpha-PVP, 4-FA, DMAA, bupropion, phenethylamine,
>   tyramine, theobromine, theophylline, norepinephrine, epinephrine, lisdexamfetamine, atomoxetine).
> - Segoe UI font + `scripts/dev-build.ps1` (fixes the vcvars VCPKG_ROOT override).
>
> **PHASE C (REAL CHEMISTRY) DONE — no more fakes in the app.** Key correction: **rdkit is NOT a vcpkg
> port** (verified on latest master; neither are openbabel/indigo/inchi), so the plan's "link RDKit via
> vcpkg" was impossible. Instead built an **in-house cheminformatics engine** (`src/chem`, zero external
> deps): a real SMILES parser -> molecular graph -> descriptors (Hill formula, MW, Lipinski HBD/HBA, Veber
> rotatable bonds, ring perception, **Ertl TPSA**, fraction-Csp3, Crippen-style logP), **Morgan/ECFP
> fingerprints + Tanimoto**, and graph-based **functional-group perception**. A **RealBackend**
> (`src/modules`) computes EVERY property from structure and drives all modules; the app links it (fakes
> retired from the app, kept only for legacy tests). **ctest = 25/25 green**, including formula-match over
> 14 structures, TPSA vs reference, and structure-derived ADMET goldens.
>
> **PHASE D DONE (this session).** ctest **38/38 green**. Landed:
> - **WP-1 · 3D embedding** (`src/chem/Embed3D.*`): distance-geometry conformers (covalent-radius / VSEPR
>   distance matrix -> classical metric-matrix MDS via **eigen3** -> bonded+angle+vdW steepest-descent,
>   explicit H), no RDKit. Deterministic, NaN-free, all 29 compounds.
> - **WP-2 · 3D molecular viewer** (`src/render/MolViewport.*`): instanced DX11 ball-and-stick / spacefill
>   (CPK colors, split-color cylinder bonds), off-screen RT -> ImGui image, orbit/zoom/pan, legend.
>   Embedded in the Structure Workbench (replaced the 2D schematic) and the Docking pose overlay.
> - **WP-3 · real docking** behind the frozen `contracts/IDockingBackend.h`: Vina/smina subprocess +
>   rigid-PDBQT writer (clean ROOT, avoids the tree.h gotcha) + REMARK-VINA-RESULT parser + 29 CNS
>   receptor presets + best-effort engine provisioner. **gnina disabled** (no Win/CUDA build). With no
>   engine on this host it degrades to the clearly-labeled descriptor estimate (poses still carry the
>   embedded ligand so the viewer has geometry).
> - **WP-4 (ALL 4):** full **Wildman-Crippen logP** atom typing (+regression test); a **Molecule Input**
>   panel (free-text SMILES -> full analysis + 3D); **SQLite** run-history persistence (live Runs panel,
>   `run_history` table under `%APPDATA%/StimLab/stimlab.db`); and the **real LLM agent loop** (below).
>
> **WP-4 Agent DONE (this session). ctest 42/42 green.** A real tool-calling assistant replaces the
> canned buttons:
> - Frozen contracts `contracts/ILlmProvider.h` + `contracts/IAgentTools.h` (ITool/IToolRegistry).
> - New `src/agent` static lib: a **MockProvider** (deterministic, offline default - keyword-routes a
>   panel + refuses synthesis asks), an **AnthropicProvider** (libcurl + **SSE** streaming of the Messages
>   API; all networking `#ifdef STIMLAB_HAVE_SCIENCE`), a threaded **Agent** tool-calling loop
>   (autopilot + ask-first), and the safety system prompt.
> - `AppShell` owns the agent + providers + registry + `Config`. A **thread-safe UI-action inbox**
>   marshals worker-thread tool calls back to the UI thread so tools never touch ImGui. Bound tools:
>   `highlight_panel`, `navigate_ui`, `list_panels`, `get_active_compound`, `what_can_stimlab_do`. The
>   assistant panel streams the live loop (canned quick-prompts kept, now routed through the agent).
> - Settings persist provider/model/mode + the **DPAPI-encrypted** API key via `Config` + `Secrets`.
> - **Both presets build**: default `windows` is curl-free (agent falls back to MockProvider); the live
>   Anthropic path builds under `windows-science`.
>
> **Build gotchas learned this session:**
> - **`curl[ssl]` = Schannel on Windows here.** The vcpkg curl 8.20 port has no `schannel` feature; its
>   `ssl` meta-feature pulls `sspi` (Schannel) on Windows, so `curl[ssl]` builds with **no OpenSSL**
>   (~26s). The "long OpenSSL build" warning in the old plan does not apply. (Also dropped the dead
>   `rdkit` from the `science` feature - it is not a port in this registry.)
> - Default model is `claude-opus-4-8`; the request omits `temperature`/`top_p`/`budget_tokens` (Opus 4.x
>   rejects sampling params with a 400). Model is a free-text Settings field (Haiku 4.5 is snappier).
> - A **live** smoke test needs a real key entered in Settings (not done here - no key on this host).
> - Legacy Python-backend `stimlab.db` can pre-exist in `%APPDATA%/StimLab` with a clashing `runs`
>   table - the native store uses `run_history` to avoid it.

---

## 1. What StimLab is

A native **CNS-stimulant computational drug-discovery / pharmacology-prediction suite** — an open,
single-`.exe` alternative to Schrödinger Maestro/Glide. It predicts **what a compound is and does**:
structure/properties, target **docking** (ligand→protein binding affinity at DAT/NET/SERT/TAAR1 etc.),
ADMET/metabolism, molecular **stability**, structural+pharmacological **similarity** to a known-substance
library, and **legal-analog "substantially similar" scorecards** for policy/forensic use — all behind a
clean GPU-rendered GUI with an AI agent that drives and explains the app.

**Stack (locked):** C++20 · Dear ImGui (docking) on **DirectX 11** + Win32 · RDKit (linked) · ONNX
Runtime (prebuilt GPU) · AutoDock **Vina**/AutoDock-GPU/gnina docking engines (real, from source) ·
SQLite + content-addressed artifact store · libcurl/WebView2 · single self-provisioning exe writing all
state to `%APPDATA%/StimLab`. Full architecture: [plan.md](plan.md) §3–§6.

---

## 2. ⚠️ Safety boundary — NON-NEGOTIABLE, carry this forward

**IN SCOPE — predict what a compound *is and does*:** structure/property analysis; docking =
ligand→protein **binding affinity** (pharmacology/activity); ADMET/metabolism incl. harmful metabolites
& drug–drug interactions; **absorption / pharmacokinetics** (HIA, oral F%, Caco-2 permeability, BBB/CNS
penetration, P-gp efflux — added this session to narrow candidates); **molecular-stability** scoring;
**similarity** to known substances; **legal status & analog scorecards**; library/workflows/runs/archive;
an AI assistant that explains, navigates, and highlights UI.

**OUT OF SCOPE — never build, even if asked:** synthesis routes/steps, reaction conditions, precursor
selection, "ease of manufacture"/synthesizability/**manufacturability** scoring, or any "how to make it"
guidance. Docking "binding energy" is **target binding (pharmacology)** and must never be repurposed as a
make-it signal. The **stability score replaces** the originally-requested manufacturability score. This
boundary was established with the user across the planning conversation and is final.

---

## 3. Authoritative documents

| Doc | What it is |
|---|---|
| [plan.md](plan.md) | **Master spec** (13 sections): functional surface, architecture, tech stack, `%APPDATA%` layout, repo tree, execution model, work packages (WP-A…WP-N), golden facts (§9), risk register (§10), verification (§12). The source of truth for *what* and *how*. |
| `C:\Users\nigga\.claude\plans\i-want-you-to-deep-plum.md` | **Approved build plan**: adopts plan.md + adds the new in-scope modules (stability/similarity/legal/structure-workbench/assistant-highlighting), the safety boundary, and the phase sequence (A→E). |
| **This file** | Live session state + resume point. |

---

## 4. Verified environment (probed this session — facts, not assumptions)

| Tool | Status | Notes |
|---|---|---|
| Working dir | `c:\Users\nigga\Desktop\Stimulant-Laboratory` | git repo (initialized this session) |
| OS / shell | Windows 10 Pro 19045 · PowerShell 7 · Bash available | not admin |
| **MSVC v143** | ✅ `cl` 19.44.35226 (x64) | VS **Community 2022** at `C:\Program Files\Microsoft Visual Studio\2022\Community`; Build Tools 2022 also present. **`cl`/`nvcc` are NOT on global PATH** — must enter the VS dev environment (vcvars64) before Ninja builds. |
| **CMake** | ✅ 4.3.1 | CMake **4.x** — old ports may trip policy minimums; preset sets `CMAKE_POLICY_VERSION_MINIMUM=3.5` as a guard. |
| **Ninja** | ✅ 1.13.2 | |
| **vcpkg** | ✅ 2026-05-27 (`d5b6777…`) at `C:\Users\nigga\vcpkg` | `VCPKG_ROOT` set & persisted (`setx`). Toolchain file present. Baseline pinned in `vcpkg.json`. |
| **git** | ✅ 2.53.0 | user `StimLab` / `erikmeeks2012@gmail.com` |
| **GPU** | ✅ RTX 3080 Ti, driver 596.36 | great for the eventual GPU docking path |
| **CUDA Toolkit** | ⚠️ **runtime only** | `C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.2` has `bin`+`lib` but **no `nvcc`**. Full toolkit must be installed before Phase D (AutoDock-GPU/ONNX-GPU). Note: v13.2 may mismatch ONNX Runtime's expected CUDA 12 — verify at Phase D. |
| Python/Node | Python 3.11.9, Node 24.14.1 present | **not used** by the native build (leftover from the rejected Python-stack option). |

**Critical build note:** every `cmake`/`ninja`/`ctest` invocation must run inside the MSVC environment.
Either run from a "Developer PowerShell for VS 2022", or wrap commands:
```powershell
cmd /c "\"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat\" && cmake --preset windows && cmake --build --preset windows && ctest --preset windows --output-on-failure"
```

---

## 5. Current on-disk state (exactly what exists)

**Present & done:**
- [plan.md](plan.md) — master spec (pre-existing).
- [CMakeLists.txt](CMakeLists.txt) — **root superbuild**, complete. Sets C++20, MSVC warnings, options
  (`STIMLAB_ENABLE_CUDA/SCIENCE/BUILD_TESTS`), finds vcpkg deps, aliases `yaml-cpp::yaml-cpp`, and
  `add_subdirectory()`s `src/core`, `src/data`, `src/contracts`, `src/fakes`, `src/render`, `src/ui`,
  `src/app`, `tests`. ⚠️ **It references those subdirs, but none exist yet** — so configure will fail
  until they're created (see §6).
- [CMakePresets.json](CMakePresets.json) — presets `windows` (light/skeleton deps), `windows-science`
  (adds RDKit/curl via vcpkg `science` feature), `windows-cuda` (adds CUDA). Ninja + vcpkg toolchain,
  `RelWithDebInfo`, triplet `x64-windows`.
- [vcpkg.json](vcpkg.json) — manifest. **Base deps** (fast): nlohmann-json, spdlog, fmt, sqlite3,
  yaml-cpp, imgui[docking-experimental,dx11-binding,win32-binding], implot, catch2. **`science` feature**
  (heavy, Phase C): rdkit, curl[ssl].
- [.gitignore](.gitignore), [.gitattributes](.gitattributes) — done.

**NOT yet created (this is the resume work — see §6):**
- `src/core/` — AppPaths, Log(+minidump), Error/`Result<T>`, Hash, EventBus, Config, Secrets(DPAPI) `+ CMakeLists.txt`
- `src/data/` — domain types (Molecule, Ligand, Receptor, Pose, Run, Artifact, Provenance) `+ CMakeLists.txt`
- `src/contracts/` — frozen interface headers (IModule, IDockingBackend, IStorage/IArtifactStore,
  IJobSystem, ILlmProvider/ITool/IToolRegistry, IWebTools, IMolRenderer/IViewport, UiCommand/UiEvent) `+ CMakeLists.txt`
- `src/fakes/` — thick fakes returning plausible data `+ CMakeLists.txt`
- `src/render/` — DX11 device + ImGui bootstrap `+ CMakeLists.txt`
- `src/ui/` — app shell, nav rail, dockspace, feature panels `+ CMakeLists.txt`
- `src/app/` — `WinMain.cpp`, lifecycle `+ CMakeLists.txt`
- `tests/` — Catch2 scaffolding + core tests `+ CMakeLists.txt`
- `docs/ARCHITECTURE.md`, `docs/contracts/*.md`, `docs/workpackages/WP-*.md`
- `scripts/dev-build.ps1` (vcvars wrapper), `README.md`

> Tip: to get configure passing incrementally, you can comment out `add_subdirectory()` lines in
> [CMakeLists.txt](CMakeLists.txt) for dirs not yet written, or (better) create every referenced subdir
> with at least a minimal `CMakeLists.txt` before the first configure.

---

## 6. Resume point — Phase C (Phase B is DONE)

**Phase B is complete and verified** (see the STATUS banner up top). Build with the new wrapper that
fixes the VCPKG_ROOT override gotcha: `./scripts/dev-build.ps1 -Test`.

### ⛔ Phase C blocker discovered this session (resolve FIRST)
`windows-science` fails at `vcpkg install` with **"the baseline does not contain an entry for port
rdkit"**, and `vcpkg search rdkit` returns nothing (while `boost-system` resolves). Findings:
- The user vcpkg at `C:\Users\nigga\vcpkg` is a **shallow clone** whose registry snapshot **omits the
  `rdkit` port** entirely (2833 ports on disk, no `ports/rdkit`). vcpkg itself prints "Run `git pull`".
- **`vcvars64.bat` overrides `VCPKG_ROOT`** to the VS-bundled vcpkg (`...\VC\vcpkg`), which ALSO lacks a
  resolvable rdkit baseline. The Phase-B build "worked" only because base ports exist in both.

**Fix options (pick one):**
1. **Update the user vcpkg registry** so rdkit appears, then pin the baseline:
   `git -C C:\Users\nigga\vcpkg pull` (or `git fetch --unshallow`), then set `builtin-baseline` in
   `vcpkg.json` to the new HEAD and run `./scripts/dev-build.ps1 windows-science`.
   Expect a **long first build** (RDKit + Boost = hours); a shared vcpkg **binary cache** is recommended.
2. **Overlay/prebuilt RDKit** (avoid the source build): point at a prebuilt RDKit and wire it via an
   overlay port / `find_package`. Faster to first-light but more manual.

Then proceed to the **Phase C walking skeleton**: SMILES → RDKit (linked) → real descriptors/properties
feeding the existing Structure/Stability/Absorption/ADMET/Similarity panels (which already render against
fakes), then Vina docking. Swap fakes for real impls behind the frozen `contracts/` interfaces.

Recommended order:
1. **`src/core`** — write the foundation (real, compiling impls). Windows-specific bits guarded by
   `_WIN32`: AppPaths via `%APPDATA%`, Log via spdlog + `MiniDumpWriteDump` (link `DbgHelp`), Secrets via
   DPAPI (`CryptProtectData`, link `Crypt32`). `Result<T>` (don't rely on C++23 `std::expected`).
2. **`src/data`** — JSON-serializable domain structs (nlohmann-json) — the lingua franca.
3. **`src/contracts`** — pure-virtual interface headers (header-only INTERFACE lib). Freeze these; every
   later WP codes against them. Mirror each with a `docs/contracts/*.md`.
4. **`src/fakes`** — minimal thick fakes (e.g. FakeLibrary/FakeDocking) so panels render plausible data.
5. **`src/render` + `src/ui` + `src/app`** — DX11 device + ImGui docking bootstrap + dockspace shell with
   a nav rail and empty feature panels (Dashboard, Structure Workbench, Docking, Runs, Library,
   Metabolism, Presets, Settings, Assistant) + theme. This is the visible deliverable.
6. **`tests`** — Catch2: AppPaths resolves, Hash stable, EventBus pub/sub, `Result<T>` semantics.
7. **Configure + build + test** from a VS dev shell (see §4). Fix dep target-name mismatches as they
   surface (sqlite3 → `SQLite::SQLite3`; yaml-cpp aliased; imgui → `imgui::imgui`; implot → `implot::implot`).
8. **Commit** the green skeleton (first commit). Then **Phase C**: flip to `windows-science` preset, wire
   SMILES→RDKit→Vina→3D poses→SQLite (the gating walking-skeleton slice).

After Phase C is green, **Phase D** fans out the work packages (plan.md §8) + the new modules, and
installs the **full CUDA Toolkit** for AutoDock-GPU. **Phase E** = integration + self-provisioning packaging.

---

## 7. Golden facts & gotchas to preserve (from plan.md §9 + learned this session)

- **Docking input:** default ligand format **SDF**; **do not feed Meeko-style PDBQT to Vina** (trips a
  `tree.h` error) — convert SDF→PDBQT via `obabel.exe`.
- **ADMET goldens (test oracles):** amphetamine/methamphetamine/MDMA → **WARN** (MAO+CYP2D6);
  dopamine → **WARN** (+COMT, catechol); acetaminophen → **reactive metabolite (NAPQI)** flag;
  caffeine → **INFO**. ADMET heuristics are **not ML** (RDKit descriptors + SMARTS).
- **29 CNS target presets** exist as YAML data (transporters, dopamine/serotonin receptors, TAAR1,
  adrenergic, histamine, muscarinic, opioid, sigma1, CB1, NMDA, enzymes MAO-A/B/COMT/AChE/CYP2D6, hERG),
  each with a real PDB ref + a box from the co-crystal ligand. Plus ship a curated known-stimulant
  reference library for the similarity/legal modules.
- **CUDA self-heal:** a provisioned runtime dir's *existence* is not proof of completeness — verify the
  required file set on launch; if incomplete, wipe + reinstall (legacy hit "libcufft missing" trusting a
  stale cache). Provision the matching CUDA 12 runtime DLLs if AutoDock-GPU/gnina are dynamically linked.
- **Windows text:** keep launcher/CLI text ASCII (cp1252 console mangles em-dashes/smart quotes).
- **Build env:** `cl`/`nvcc` not on PATH — always build from a vcvars/dev shell (§4).
- **CMake 4.x:** preset sets `CMAKE_POLICY_VERSION_MINIMUM=3.5`; a port that still fails may need a
  per-port override or a baseline bump.
- **RDKit via vcpkg is an hours-long build** (RDKit+Boost) — kept behind the `science` feature so the
  Phase-B skeleton stays fast. Consider a shared vcpkg binary cache before Phase-D fan-out.
- **Static CRT (`/MT`) deferred:** plan.md wants `/MT` for redist-free packaging, but the skeleton uses
  the default dynamic `x64-windows` (`/MD`) to minimize friction; revisit at Phase E (bundle redist or
  switch to `x64-windows-static`).

---

## 8. Open questions parked for later (plan.md §13)
- gnina CNN required for v1, or post-v1? (Vina + AutoDock-GPU already give CPU+GPU docking.)
- Keep MD/FEP/QM/viz as honest typed stubs in v1? (assumed yes)
- App identity/branding — keep "StimLab"? (assumed yes)
- ONNX Runtime GPU build vs installed CUDA **13.2** version compatibility — verify at Phase D.
