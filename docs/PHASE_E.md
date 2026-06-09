# StimLab — Phase E + WP-K brief (agent tools, self-provisioning, packaging)

> Final v1 tracks. Read [handoff.md](../handoff.md) (STATUS banner + §2 safety boundary) and
> [docs/PHASE_D.md](PHASE_D.md) (real docking) first. Safety boundary is NON-NEGOTIABLE: predict what a
> compound IS/DOES; NO synthesis/route/precursor/manufacturability content anywhere.

This brief covers the two tracks that finished v1: **#3 expansive agent tools (WP-K)** and **#4 packaging /
self-provisioning (Phase E, WP-B + WP-M)**. ctest **68/68** on `windows` and `windows-static`; the science
build is green and the live `web_search` smoke test passes against DuckDuckGo.

## WP-K — expansive agent tools
The assistant can now DO things, not just navigate/explain. Tools are bound to the native services in
`AppShell::registerAgentServiceTools()` / `registerAgentWebTools()` (UI layer, so the agent lib stays
UI-free); the agent dispatches them via the frozen `IToolRegistry`.

- **Service-action tools** (always available, no network):
  - `analyze_compound` — identity + physchem (formula/MW/logP/TPSA/HBD/HBA) + stability + absorption +
    overall ADMET, for a library name/id OR a raw SMILES analyzed on the fly (`resolveAgentCompound`).
  - `screen_admet` — per-endpoint metabolism/ADMET verdicts.
  - `dock_compound` — dock into a CNS target → affinity + real/estimate + pose count.
  - `run_workflow` — kick the prep→dock DAG and focus the Workflows panel.
  - `list_runs`, `search_library`, `compare_compounds` (side-by-side props/stability/absorption/ADMET).
- **Web tools** (`src/agent/WebTools.{h,cpp}`, behind `STIMLAB_HAVE_SCIENCE`):
  - `web_search` — keyless DuckDuckGo HTML endpoint → structured hits (real URLs un-wrapped from the
    `/l/?uddg=` redirect). `web_fetch` — GET + HTML→text (scripts/styles/tags stripped, entities decoded).
  - Cached (12h) + rate-limited (≥700ms) under `%APPDATA%/StimLab/cache/web`. The HTML parsers are PURE and
    always compiled, so they're unit-tested from fixtures; the curl-free `windows` build reports the tools
    unavailable and the agent degrades. Returned content is flagged untrusted third-party text.

`AppShell::toolRegistry()` exposes the registry for enumeration / tests.

## WP-B — manifest.json + self-heal
`core/Manifest.{h,cpp}` records each provisioned component (engine binaries, prepared receptors) with
`path + size + FNV-1a content hash` into `%APPDATA%/StimLab/manifest.json`. The Provisioner writes it after
a successful provision; `verify()` re-checks existence + size + hash; `heal()` deletes corrupt files so the
next provision re-fetches them. `WinMain` self-heals on launch (logs N/M verified); Settings shows the
status + a "Verify + heal runtime" button. Download AUTHENTICITY stays with the pinned SHA-256 in
`EngineLocator` (env `STIMLAB_VINA_SHA256` or `runtime/engines/vina.sha256`); the manifest hash is fast
corruption/truncation detection. Verified live: a real provision writes vina.exe (1,203,712 B) + the four
headline receptors, each with size + hash.

## WP-M — packaging (fully-static single exe)
New `windows-static` configure preset: `x64-windows-static` triplet + static CRT (`/MT`, via
`CMAKE_MSVC_RUNTIME_LIBRARY`). The result is a **3.3 MB `StimLab.exe` that needs NO DLLs and NO VC++
redist** — a true single self-contained binary. The fast dynamic dev presets (`windows`/`windows-science`)
are untouched. `scripts/package.ps1` builds it, stages `dist/StimLab-<ver>-win-x64/` (exe + README; bundles
the fmt/spdlog/sqlite3 DLLs automatically when packaging a dynamic preset instead), and zips it (~1.7 MB).

**Verified:** the static build's ctest is 68/68; the standalone exe run from a clean `%TEMP%` dir (no DLLs,
no build tree) self-provisions and does a REAL Vina dock. First launch self-provisions the engine +
receptors into `%APPDATA%/StimLab` via the Docking panel's Provision button (needs network); until then
docking shows the clearly-labeled descriptor estimate. Real docking works even in this curl-free static
build because provisioning uses PowerShell `Invoke-WebRequest`, not libcurl — only the live LLM provider and
the web tools need the science build.

## Post-v1 enhancements — ALL FIVE LANDED (2026-06-09)
Each shipped on a continuously-green trunk (ctest 69/69 windows + windows-static), one commit per track.
See the handoff STATUS banner for the headline summary.

- **Receptors on demand** (`efd2c49`) — `AppShell::provisionTarget()` / `receptorReady()` provision ANY of the
  29 CNS presets (resolved by id OR display name) off the UI thread, manifest auto-updated; the Docking +
  Workflows panels show per-target readiness + a "Provision &lt;target&gt;" button for the 25 non-headline targets.
  `--selftest-dock --target X` now provisions the requested X. Verified: D2 −7.08, MAOB −7.50 kcal/mol real docks.
- **Science-static single exe** (`2a6ef65`) — `windows-science-static` preset (x64-windows-static + /MT + the
  science feature → static `curl[ssl]` = Schannel). ONE DLL-free ~3.85 MB exe carrying the live agent + web
  tools. `build.ps1 -Release -Science`; `package.ps1` advertises the bundled agent. Live web_search verified
  standalone from a clean dir.
- **CI + code-signing** (`a687b59`) — `scripts/ci.ps1` (local CI: build + ctest windows + windows-static +
  package) and `.github/workflows/ci.yml` (windows-latest, get-cmake + msvc-dev-cmd + vcpkg GH Actions cache).
  `scripts/sign.ps1` is an optional Authenticode hook (signs with STIMLAB_SIGN_PFX/_SHA1 via signtool; a clean
  exit-0 no-op otherwise). No git remote yet → the GH workflow is written-but-unrun; ci.ps1 is verified green.
- **WebView2 `web_fetch`** (`25fe293`) — optional headless WebView2 JS-render path
  (`src/agent/WebToolsRendered.cpp`: dedicated STA thread + COM + message pump → outerHTML → htmlToText),
  behind a new `STIMLAB_ENABLE_WEBVIEW2` (requires science). The `webview2` vcpkg port links the STATIC loader
  (`WebView2LoaderStatic.lib`) on x64-windows-static so the single exe stays DLL-free; runtime-detected
  (Evergreen Runtime) fallback to curl. `web_fetch` gained a `render` flag. Verified live: a data: URL whose
  visible text is JS-COMPUTED (`(6*7)`) → the marker is recovered (`stimlab_tests.exe "[webview2]"`).
- **GPU docking** (`3f3af7a`) — first-party **CUDA** rigid-body docking engine (`CudaBackend` + `CudaScore.cu`:
  AutoDock Vina inter-molecular scoring on the GPU) behind the lean `windows-cuda` preset (nvcc 13.3, sm_86,
  `-allow-unsupported-compiler`; MSVC C++ flags scoped to CXX so nvcc isn't handed `/permissive-`). The dead
  Compute selector is now wired (`docking::ComputeMode` Auto/GPU/CPU, persisted via Config; `realEngines()`
  orders by mode; honest GPU→Vina→estimate fallback). AutoDock-GPU/gnina native-Windows CUDA source builds
  were assessed INFEASIBLE (Unix Makefile + CUDA-13 drops pre-sm_75; gnina/Uni-Dock Linux-only) — the
  first-party CUDA engine is the honest replacement, clearly labeled rigid-body (flexible CPU Vina stays the
  accurate path). Verified live on the RTX 3080 Ti: amphetamine→DAT, 64,000 GPU poses, 9 real ranked poses,
  −3.66 kcal/mol, exit 0 (`--selftest-dock --compute gpu`).

### Still open (future, not blocking)
- Additional GPU engines behind the same `IDockingBackend` seam: Vina-GPU (OpenCL, prebuilt Windows exe,
  reuses the PDBQT+box contract) and a WSL2 subprocess to Linux Uni-Dock/AutoDock-GPU (the actual CUDA builds).
- Static cudart (`CUDA::cudart_static`) for a self-contained single-exe GPU build (a `windows-cuda-static`).
- A real git remote + a live Windows CI run; an actual code-signing certificate (unsigned trips SmartScreen).
