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

## Post-v1 follow-ups (deferred, not blocking)
- AutoDock-GPU (WP-G2) + gnina (WP-G3) CUDA engines; needs the full CUDA toolkit (nvcc).
- WebView2 headless for JS-heavy `web_fetch` pages (current fetch is a plain GET + HTML strip).
- Prepare the remaining 25 (non-headline) receptors — same code path, on demand.
- A science-static release (curl[ssl] static) bundling the live LLM + web tools into the single exe.
- CI on a Windows+CUDA runner; code-signing the exe (unsigned trips SmartScreen — expected for personal use).
