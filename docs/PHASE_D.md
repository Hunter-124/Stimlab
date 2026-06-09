# StimLab — Phase D brief (real docking + 3D viewer + improvements)

> Self-contained continuation brief. Phase A–C are DONE and verified. Read
> [handoff.md](../handoff.md) (STATUS banner + §6) and [plan.md](../plan.md) first.
> Safety boundary (handoff §2) is NON-NEGOTIABLE: predict what a compound IS/DOES;
> NO synthesis/route/precursor/manufacturability content. Stability replaces
> manufacturability. Docking "binding affinity" is pharmacology, never a make-it signal.

## Current state (don't redo)
- Build: `./scripts/dev-build.ps1 -Test` → `build/windows/bin/StimLab.exe`, **ctest 25/25 green**.
  (vcvars overrides `VCPKG_ROOT`; the wrapper forces it back to `C:\Users\nigga\vcpkg`.)
- **rdkit is NOT a vcpkg port** (verified latest master; neither openbabel/indigo/inchi). Do NOT try to
  link RDKit via vcpkg. We have our own engine in `src/chem` (SMILES→graph→descriptors, Ertl TPSA,
  Morgan/Tanimoto, functional groups). `src/modules/RealBackend` computes everything and drives the UI
  behind the frozen interfaces in `src/contracts/`.
- GUI: 3-pane DX11/ImGui (Navigator | Workspace | Assistant), 14 panels, assistant pulse-highlight,
  ~29-compound library, Segoe UI. `src/render/Dx11Device` owns the device/swapchain.
- `src/fakes` is retired from the app (kept only for legacy tests).

## WP-1 — 3D coordinate embedding  (`src/chem/Embed3D.{h,cpp}`)  [do first; both WP-2/WP-3 need it]
Generate 3D coordinates from the molecular graph (no RDKit). Pragmatic distance-geometry:
build a distance-bounds matrix from bonded (1-2), angle (1-3), and torsion (1-4) ideal distances using
covalent radii + VSEPR angles; metric-matrix embedding (eigendecompose, top-3 eigenvectors via `eigen3`
which IS in vcpkg) then a few steps of steepest-descent on a simple bonded+VdW force field. Add explicit
H positions for the viewer. Acceptance: a unit test asserting C–C bond lengths ≈1.5 Å and no atom
overlaps for amphetamine/caffeine; embed all 29 library compounds without NaNs.

## WP-2 — 3D molecular viewer  (`src/render/MolViewport.{h,cpp}` + `src/ui` integration)
DX11 renderer, GPU-instanced **impostor spheres** (atoms, CPK colors, VdW/ball radii) + **cylinders**
(bonds; split-color by element). Off-screen render target → ImGui image in the Structure Workbench
(replace the 2D `moleculeSchematic` placeholder) and the Docking panel (pose overlay). Orbit/pan/zoom
camera (drag + wheel), element legend, toggle ball-and-stick / spacefill. Reference existing impostor
shaders (3Dmol/Speck/PyMOL-style) — do not invent. Acceptance: selecting a compound shows a rotatable
ball-and-stick; ~60fps; atom pick highlights element.

## WP-3 — Real docking backends  (`src/contracts/IDockingBackend.h`, `src/modules/docking/*`, `third_party/`)
Replace the descriptor estimate in `RealDocking` with real engines behind a new `IDockingBackend`:
- **Vina** (default): provision official `vina.exe` (download + checksum into `%APPDATA%/StimLab/runtime/
  engines/`); ligand prep = WP-1 3D embed → write **PDBQT** (Gasteiger charges + atom types; write our own
  minimal PDBQT writer to avoid the OpenBabel dependency, OR provision `obabel.exe`); receptor presets =
  reuse the CNS target boxes (plan §9: 29 presets w/ PDB ref + box) under `src/presets` YAML; run Vina as
  a subprocess; parse `REMARK VINA RESULT` → scored `Pose`s. **Gotcha (plan §9):** default ligand format
  SDF; do NOT feed Meeko-style PDBQT to Vina (trips `tree.h`) — use our writer/obabel.
- **smina**: same interface, smina binary (custom scoring). 
- **gnina** (CUDA CNN): STRETCH — no official Windows build; fallback = export CNN to ONNX + ONNX Runtime
  (Microsoft prebuilt GPU build, NOT vcpkg-from-source). Ship disabled if it doesn't converge.
Wire `Settings` → engine選択 + GPU mode. Fall back to the descriptor estimate (clearly labeled) when no
engine is provisioned. Acceptance: dock a known ligand into DAT/SERT box → finite scores + ranked poses
rendered in the WP-2 viewport; engine selectable; absent-engine path degrades gracefully.

### WP-3 STATUS — REAL DOCKING LANDED + VERIFIED LIVE (ctest 48/48)
Docking now produces **real `real=true` engine poses** end-to-end. What was added this session:
- **WP-F receptor prep** (`src/modules/docking/ReceptorPrep.{h,cpp}`): an in-house **PDB → rigid receptor
  PDBQT** writer (no RDKit, no hard OpenBabel link). Strips waters + crystallization additives + the
  co-crystal ligand, infers elements, assigns AutoDock4 atom types (aromatic A / NA / OA / SA), writes a
  heavy-atom rigid receptor (Vina scoring is partial-charge-independent, so heavy-atom typing docks).
  `obabel.exe` is used instead when present (adds polar H). `ensureReceptor()` fetches the preset's PDB
  from **RCSB** (`files.rcsb.org`), prepares it, and caches `runtime/receptors/<id>.pdbqt`.
- **Box from the structure's OWN frame.** `receptorBoxFromPdb()` derives the box center from the
  **co-crystal ligand centroid** (largest non-water HETATM ≥ 8 atoms) in the fetched PDB's real
  coordinates, written as a `REMARK STIMLAB_BOX` and used to override the preset's literature-style center
  (which is in a *different* frame — e.g. DAT preset `(-2.1,12.4,-6.8)` vs 4M48 ligand `(-57.5,-12.6,49.0)`).
  Without this the box misses the receptor and Vina returns nothing. This is the load-bearing fix.
- **Engine provisioning** (`EngineLocator`): `ensureVina(true)` downloads the official `vina_1.2.5_win.exe`
  with a **size check** against the published 1,203,712 B (GitHub publishes no digest) and an **optional
  runtime SHA-256 pin** (env `STIMLAB_VINA_SHA256` or `runtime/engines/vina.sha256`, no rebuild needed).
  `ensureObabel()` is locate-only (OpenBabel ships as an installer; the built-in writer is the default).
- **Off-thread provisioning + docking.** `docking::Provisioner` runs the slow work on a worker thread;
  `AppShell::dockFor()` runs each dock on a worker and caches by molecule+target so the **UI never blocks**
  and a real dock never runs on the UI thread. Vina is given a fixed `--seed 1` so results are reproducible.
- **Headless acceptance hook:** `StimLab.exe --selftest-dock [--smiles S] [--target T]` provisions then
  docks through the real path; exit 0 = real dock. **Verified live on this host:** amphetamine→DAT =
  **−4.83 kcal/mol, 9 ranked poses, real=true** (reproducible); MDMA→SERT = **−4.26, 9 poses, real=true**.
  GUI Docking panel shows the same real poses + 3D viewer (screenshot-verified).

**Remaining to verify / future:** only the 4 headline receptors (DAT/NET/SERT/TAAR1) are provisioned up
front; the other 25 presets prepare on demand (same code path, not yet each individually exercised).
AutoDock-GPU (WP-G2) and gnina (WP-G3) are still future; obabel-based protonation is an optional quality
upgrade if `obabel.exe` is dropped into `runtime/engines`.

## WP-4 — improvements (pick up as capacity allows)
- **Full Wildman–Crippen logP** atom typing (replace the estimate in `chem/Descriptors.cpp::crippenLogP`);
  add a regression test vs reference logP for the library (±0.7).
- **SMILES/name input panel**: free-text candidate → `chem::parseSmiles` → full analysis (extends the
  Analog Explorer to arbitrary structures, not just slider-modeled ones).
- **Persistence**: use the already-linked `sqlite3` — store runs/archive/exported reports in
  `%APPDATA%/StimLab/stimlab.db` (schema in `src/storage`); make the Runs panel live.
- **Agent** [DONE]: `contracts/ILlmProvider.h` + `contracts/IAgentTools.h`; `src/agent` lib with a
  deterministic **MockProvider** (offline default) and an **AnthropicProvider** (libcurl + SSE Messages
  API, behind `STIMLAB_HAVE_SCIENCE`); a threaded tool-calling **Agent** (autopilot + ask-first). Tools
  `highlight_panel` / `navigate_ui` / `list_panels` / `get_active_compound` / `what_can_stimlab_do` bind
  to `AppShell` via a thread-safe UI-action inbox (worker threads never touch ImGui). The API key is
  persisted DPAPI-encrypted via `Settings`/`Secrets`. ctest 42/42. NOTE: vcpkg `curl[ssl]` resolves to
  **Schannel** (sspi) on Windows here - no OpenSSL build; and `rdkit` was removed from the `science`
  feature (not a port in this registry). The system prompt + tool set enforce the safety boundary:
  the agent refuses synthesis/route/precursor/manufacturability requests.

## Execution notes
- Freeze `IDockingBackend` + the `Embed3D` API FIRST, then WP-2 and WP-3 can run as parallel subagents
  (worktree isolation recommended: they both touch CMake/UI). Keep `ctest` at 25/25+ green throughout —
  add tests per WP. Build only via `./scripts/dev-build.ps1`.
- Verify the GUI by launching the exe and screenshotting (the prior session used a GDI window-capture
  PowerShell snippet since the computer-use resolver can't see an uninstalled exe).
- New vcpkg deps available and useful: `eigen3` (embedding), `curl` (agent, behind the `science` feature),
  `sqlite3` (already linked). Do NOT add rdkit/openbabel (not in this registry).
