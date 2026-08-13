# External data sources and network endpoints

Every outbound endpoint in BioCAD was found by grepping this repository, and each row cites the
`file:line` where the URL is constructed. If an endpoint is not in this table, the application
does not contact it.

Two independent transports exist, and they matter for the table below:

- **libcurl**, compiled only when `BIOCAD_HAVE_SCIENCE` is defined (vcpkg feature `science`,
  CMake option `BIOCAD_ENABLE_SCIENCE`). Without it the agent's network code compiles to a stub
  (`src/agent/AnthropicProvider.cpp:8,293-296`, `src/agent/WebTools.cpp:18-20,300-306`).
- **PowerShell `Invoke-WebRequest`**, shelled out from the docking provisioning path. This is
  *not* gated on `BIOCAD_ENABLE_SCIENCE` - it is a `#if defined(_WIN32)` path with no libcurl
  dependency at all (`src/modules/docking/EngineLocator.cpp:229-260`,
  `src/modules/docking/ReceptorPrep.cpp:153-172`).

## Endpoints

| Source | Exact endpoint | When contacted | Needs `BIOCAD_ENABLE_SCIENCE` | Key required | Licence / terms of the data | Cacheable on disk? | Committable to this repo? |
|---|---|---|---|---|---|---|---|
| RCSB Protein Data Bank | `https://files.rcsb.org/download/<PDBID>.pdb` (`src/modules/docking/ReceptorPrep.cpp:420`) | Only from `ensureReceptor(..., allowDownload=true)`, i.e. when the user presses Provision or a target is prepared on demand; never on the cache-only hot docking path (`ReceptorPrep.cpp:397-400`, `src/modules/RealBackend.cpp:403-409`) | No (PowerShell transport) | No | PDB entries are released without restriction into the public domain; RCSB asks that the entry ID and its primary citation be quoted. Same statement as `NOTICE:56` | Yes - already cached to `runtime\receptors\<PDBID>.pdb` and the derived `<targetId>.pdbqt` | Not committed today. Structures are public domain, so committing a small fixture is permissible; the repo instead keeps synthetic PDB text in `tests/test_docking.cpp` |
| AutoDock Vina 1.2.5 Windows binary | `https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_win.exe` (`src/modules/docking/EngineLocator.cpp:44-45`) | Only from `ensureVina(allowDownload=true)`, i.e. the user pressed Provision. The startup probe is locate-only and touches no network (`EngineLocator.cpp:204-212`, `src/modules/docking/Provisioning.h:44-46`) | No (PowerShell transport) | No | AutoDock Vina is Apache-2.0 (Trott & Olson 2010, J Comput Chem 31:455-461); stated in `NOTICE:57`. BioCAD downloads it into the user's app-data directory and redistributes nothing | Yes - it *is* a cached artefact: `runtime\engines\vina.exe`, recorded in `manifest.json` | No. Binaries are never committed; provisioning is the distribution mechanism |
| Vina-GPU (OpenCL) binaries and kernel source | `https://raw.githubusercontent.com/DeltaGroupNJUPT/Vina-GPU/main/<file>` for the 17 files in `vinaGpuFiles()` (`src/modules/docking/EngineLocator.cpp:70-89`) | Only from `ensureVinaGpu(allowDownload=true)` (`EngineLocator.cpp:315-320`) | No (PowerShell transport) | No | A derivative of AutoDock Vina; the upstream repository's exact licence is **not verified in-repo** - `NOTICE:58` says to consult upstream before redistributing. Treat as unverified | Yes - `runtime\engines\vina-gpu\`, plus a locally compiled `Kernel2_Opt.bin` | No |
| Anthropic Messages API | `https://api.anthropic.com/v1/messages` (`src/agent/AnthropicProvider.cpp:243`) | Only when the user sends a message in the Assistant panel with the Anthropic provider selected and a key configured; `ready()` is false without a key (`AnthropicProvider.cpp:35-40`) | **Yes** - without it `send()` is a stub that returns an error (`AnthropicProvider.cpp:293-296`) | **Yes** - the user's own key, sent as `x-api-key` (`AnthropicProvider.cpp:240-241`), stored locally via DPAPI (`src/core/Secrets.cpp`) | Anthropic's commercial API terms; the user's own account and key. Model responses are the user's, subject to those terms | Not cached. Requests are streamed and nothing is written to the web cache | No |
| DuckDuckGo HTML search | `https://html.duckduckgo.com/html/?q=<escaped query>` (`src/agent/WebTools.cpp:338`) | Only when the assistant runs the `web_search` tool during a user-initiated turn | **Yes** | No | Used under DuckDuckGo's terms of use (`NOTICE:61`); the terms themselves are **not verified in-repo**. Result snippets are third-party content | Yes - parsed hits are cached as JSON under `cache\web\<hash>.json` for 12 hours (`WebTools.cpp:322-335`) | No. Third-party page text with unclear redistribution rights |
| Arbitrary user/agent-chosen web page (plain fetch) | Whatever URL is passed to `web_fetch`; issued via `curl_easy_setopt(c, CURLOPT_URL, url)` (`src/agent/WebTools.cpp:170`) | Only when the assistant runs the `web_fetch` tool during a user-initiated turn | **Yes** | No | Whatever the target site's terms are. Unknowable in advance; treat every fetched page as untrusted third-party text with no redistribution rights | Yes - extracted text cached as `cache\web\<hash>.txt` for 12 hours (`WebTools.cpp:373-379`) | No |
| Arbitrary web page, JavaScript-rendered | Same URL, loaded in a headless WebView2 instead of curl (`src/agent/WebToolsRendered.cpp:1-15`; invoked from `WebTools.cpp:381-384`) | Only when `web_fetch` is called with `renderJs` and the machine-level Evergreen WebView2 Runtime is present; otherwise it returns `nullopt` and falls back to curl | Compiled under the `webview2` feature; the calling tool still requires `science` | No | As above. The page's JavaScript runs in an isolated, invisible, host-object-free WebView2 | Yes - same 12-hour text cache, plus a WebView2 browser profile under `cache\web\wv2` (`WebToolsRendered.cpp:76-80`) | No |
| DigiCert timestamp authority | `http://timestamp.digicert.com`, overridable via `BIOCAD_SIGN_TS_URL` (`scripts/sign.ps1:37`) | **Build/release only.** Authenticode timestamping during signing. The application never contacts it | n/a | No | Release-signing infrastructure, not data | n/a | n/a |

Two entries that are deliberately *not* endpoints:

- **smina** is located on `PATH` or in `runtime\engines` and is **never downloaded**
  (`src/modules/docking/EngineLocator.cpp:56-57`). BioCAD only executes a binary the user supplies.
- **Open Babel** is locate-only by design; `ensureObabel()` contains no URL and says so
  (`EngineLocator.cpp:393-405`).

Data packs under `assets/packs/*.json` are static files shipped beside the executable and
loaded from disk (`src/packs/Pack.h:4-5`). Loading a pack performs no network access. The
identifiers packs carry (ChEMBL IDs, PubChem CIDs, InChIKeys, PDB IDs) are recorded with their
terms in `NOTICE:70-83`; none of them is fetched at runtime by any code in this tree.

## What is cached, and where

Everything persistent lives under one root, `%APPDATA%\BioCAD`, resolved by
`AppPaths::instance()` with the `BIOCAD_HOME` environment override
(`src/core/AppPaths.h:11-27`). Nothing is written to Program Files.

| Path | Contents | Origin |
|---|---|---|
| `runtime\engines\` | `vina.exe` (downloaded), any `smina.exe` / `obabel.exe` you dropped in, optional `vina.sha256` pin | `EngineLocator.cpp:129-131` |
| `runtime\engines\vina-gpu\` | `Vina-GPU.exe`, `Vina-GPU-K.exe`, the OpenCL kernel source, the example complex, and the locally compiled `Kernel2_Opt.bin` | `EngineLocator.cpp:133-135,298-390` |
| `runtime\receptors\` | Fetched `<PDBID>.pdb` source structures and prepared `<targetId>.pdbqt` receptors carrying the `REMARK BIOCAD_BOX` line | `ReceptorPrep.cpp:191-193,411-414` |
| `runtime\selftest-dock.txt` | Report from `--selftest-dock` | `src/app/WinMain.cpp:130` |
| `manifest.json` | Path, size and FNV-1a hash of every provisioned engine and receptor; the self-heal source of truth | `src/modules/docking/Provisioning.cpp:15-30`, `src/core/Manifest.h:22-40` |
| `cache\` | Docking scratch: `lig_<target>.pdbqt` and `dock_<target>.pdbqt` per run, plus the workflow DAG's node cache | `src/modules/docking/Backends.cpp:96-102`, `src/workflow/Dag.cpp:107-110` |
| `cache\web\` | 12-hour agent web cache: `<hash>.json` search results, `<hash>.txt` page text | `src/agent/WebTools.cpp:202-207,322,373` |
| `cache\web\wv2\` | WebView2 browser profile for the JS-rendered fetch path | `src/agent/WebToolsRendered.cpp:76-80` |
| `packs\` | Your own data packs; override built-ins by pack id | `src/packs/Pack.h:4-5` |
| `presets\`, `artifacts\`, `logs\biocad.log`, `biocad.db`, `config.json` | Application state, run history, settings | `src/core/AppPaths.h:20-27` |

### Purging

All of it is plain files under one directory; deletion is the supported purge.

```powershell
# Everything network-derived, keeping settings and run history:
Remove-Item -Recurse -Force "$env:APPDATA\BioCAD\cache"              # web + docking scratch
Remove-Item -Recurse -Force "$env:APPDATA\BioCAD\runtime\receptors"  # fetched PDBs + prepared PDBQTs
Remove-Item -Recurse -Force "$env:APPDATA\BioCAD\runtime\engines"    # provisioned binaries
Remove-Item -Force         "$env:APPDATA\BioCAD\manifest.json"       # provisioning record

# Everything, including settings, API key and run history:
Remove-Item -Recurse -Force "$env:APPDATA\BioCAD"
```

Notes on purging:

- Delete `manifest.json` alongside `runtime\`, or leave it and let `selfHealManifest()` reconcile
  on the next launch - it deletes components whose size or hash no longer matches and reports the
  missing ones so provisioning re-fetches them (`Provisioning.cpp:32-37`, `Manifest.h:52-63`).
- The web cache is time-bounded anyway: entries older than 12 hours are ignored and refetched
  (`WebTools.cpp:323,374`).
- The stored Anthropic API key lives in `config.json`, DPAPI-encrypted (`src/core/Secrets.cpp`);
  deleting that file removes it.
- The migration path already does one targeted purge for you: moving a pre-rename
  `%APPDATA%\StimLab` root deletes `runtime\receptors` wholesale, because those PDBQTs carry the
  old box marker (`src/core/AppPaths.cpp:72-73`).

## Nothing is contacted without a user action

- Startup runs a **locate-only** probe. `Provisioner::start(allowDownload=false, ...)` checks the
  filesystem and `PATH` and returns; every download call site is guarded by `allowDownload`
  (`src/modules/docking/Provisioning.h:44-46`, `EngineLocator.cpp:204-212,315-320`,
  `ReceptorPrep.cpp:397-400`). The UI states this and offers an explicit Provision button
  (`src/ui/Panels.cpp:1002-1004`).
- The docking hot path is cache-only. `RealDocking::dockDetailed()` calls
  `locatePreparedReceptor()`, never `ensureReceptor()`, so running a dock never fetches anything
  (`src/modules/RealBackend.cpp:397-409`).
- The assistant's `web_search` / `web_fetch` tools only run inside a turn the user started, and
  only in a `science` build; otherwise they report unavailable and the agent degrades
  (`src/agent/WebTools.cpp:300-314`).
- The LLM provider needs a key the user entered; with no key `ready()` is false and nothing is
  sent (`src/agent/AnthropicProvider.cpp:22-40`). Without the `science` feature the default
  offline MockProvider is used and there is no transport at all
  (`src/agent/CMakeLists.txt:2-5`).
- Outbound requests on the curl path are rate-limited process-wide to one per 700 ms
  (`WebTools.cpp:149-159`) and results are cached for 12 hours, so a chatty agent cannot hammer a
  third-party host.
