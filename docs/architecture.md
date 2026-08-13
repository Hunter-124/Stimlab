# Architecture

BioCAD is one Win32 process, one executable, and thirteen static libraries. There is no
service, no IPC, and no scripting runtime. This document describes the seams that matter when
you change something.

## Target graph

Read off the `target_link_libraries` calls in each `src/*/CMakeLists.txt`; arrows point from a
target to what it links.

```text
                        BioCAD (src/app, WIN32 exe)
                              |
        +---------------------+----------------------+
        |                     |                      |
    biocad_ui            biocad_render         biocad_modules
        |                     |                      |
   +----+----+----+-----+     +-- biocad_chem   +----+----+-----+------+
   |    |    |    |     |         d3d11 dxgi    |    |    |    |      |
   |    |    |    |     |         d3dcompiler   |    |    |    |      |
   v    v    v    v     v         windowscodecs v    v    v    v      v
 render chem modules workflow    spdlog     packs storage workflow  chem
        agent                                  |     |
                                               |     |
        +--------------------------------------+-----+
        |
   biocad_contracts (INTERFACE) --> biocad_data, biocad_chem
        |
   biocad_core  --> nlohmann_json, spdlog, fmt, DbgHelp, Crypt32, Shell32, Ole32
   biocad_data  --> nlohmann_json
   biocad_chem  --> Eigen3
```

| Target | Sources | Responsibility |
| --- | --- | --- |
| `biocad_core` | `src/core` | Paths, config, logging, minidumps, DPAPI secrets, manifest, generated `Version.h` |
| `biocad_data` | `src/data` | Domain DTOs, `Provenance`, `Quantity`, JSON (de)serialisation |
| `biocad_chem` | `src/chem` | SMILES parsing, descriptors, analysis, 3D embedding, the ADMET/PK model (header-only) |
| `biocad_contracts` | `src/contracts` | Header-only frozen interfaces plus the `Services` struct |
| `biocad_packs` | `src/packs` | Versioned JSON data packs: schema, loader, merge and error reporting |
| `biocad_storage` | `src/storage` | SQLite run store |
| `biocad_workflow` | `src/workflow` | Job system, cancel tokens, the content-cached DAG |
| `biocad_agent` | `src/agent` | LLM providers, tool registry, system prompt, web tools |
| `biocad_modules` | `src/modules` | Real backend services, docking backends, provisioning, receptor prep |
| `biocad_render` | `src/render` | DX11 device, WIC back-buffer capture, molecular viewport |
| `biocad_ui` | `src/ui` | Theme, `AppShell`, panels |
| `BioCAD` | `src/app` | `wWinMain`, CLI parsing, window, main loop |
| `biocad_tests` | `tests` | One Catch2 binary, registered via `catch_discover_tests` |

`biocad_contracts` is an INTERFACE target: the contracts are headers, and everything that
implements or consumes them inherits `biocad_data` and `biocad_chem` through it.

## The Services seam

`src/contracts/Services.h` is the whole dependency-injection mechanism:

```cpp
struct Services {
    ILibrary*           library    = nullptr;
    IStabilityModule*   stability  = nullptr;
    IAdmetModule*       admet      = nullptr;
    IAbsorptionModule*  absorption = nullptr;
    ISimilarityModule*  similarity = nullptr;
    ILegalModule*       legal      = nullptr;
    IDockingModule*     docking    = nullptr;
    IRunStore*          runs       = nullptr;

    [[nodiscard]] bool valid() const;
};
```

It is a non-owning bundle of interface pointers, populated by `RealBackend::services()`
(`src/modules/RealBackend.cpp`). The UI codes only against the interfaces, so a module can be
replaced without touching a panel.

**There is deliberately no fake backend.** An earlier revision carried one - a second
implementation of every contract, roughly a thousand lines of parallel science - and the test
suite exercised it. That is worse than having no tests: it validates a double while the product
ships the original, and the two drift silently. `tests/test_backend.cpp` now runs
`RealBackend`. Hermeticity comes from the design instead of from a double: the catalog is data
on disk, every property is computed from SMILES, and the docking hot path is cache-only so it
never downloads and never spawns an unprovisioned engine.

This is **compile-time** dependency injection. There is no plugin loader, no registry, no
dynamic discovery, and adding one is not an improvement anybody has asked for: the binary is
meant to be a single self-contained exe.

### Checklist: adding a module kind

There is no shortcut, and skipping a step fails in a different place each time. In order:

1. **Contract.** Add the interface to `src/contracts/IModules.h` (or a sibling header, as
   `IDockingBackend.h` does). Pure virtual, virtual destructor, no implementation.
2. **`Services` member.** Add the pointer to `src/contracts/Services.h` **and** extend
   `Services::valid()`. A member left out of `valid()` is a null-pointer crash waiting for the
   first panel that uses it.
3. **Real implementation.** Add the class in `src/modules/RealBackend.cpp`, hold it in
   `RealBackend::Impl`, and assign it in `RealBackend::services()`.
4. **Determinism.** The implementation must be deterministic for a given input: no clock, no
   randomness, and no network on a read path. Do NOT add a second implementation to make it
   testable - if the real one cannot be tested, that is a defect in the real one.
5. **Panel row.** Add a `PanelInfo{id, title, group}` entry to the table in
   `src/ui/AppShell.cpp`. Panel **ids are persisted** in `imgui.ini` and are used by the agent
   tools and `--shot-panel`, so an id is an API: never rename one.
6. **Router branch.** Add the `else if` in the panel router in `src/ui/AppShell.cpp` that calls
   your `panels::` function, and declare that function in `src/ui/Panels.h`.
7. **Agent tool (optional).** Register a `FunctionTool` with a JSON schema in
   `AppShell::registerAgentServiceTools()`. Tool ids are brand-free and stable.
8. **Tests.** Add cases in `tests/test_backend.cpp` against the real implementation. If it
   needs an engine or a network, the module needs a cache-only path, not a stand-in.

Every derived number the new module emits must carry a `Provenance`; see
[provenance.md](provenance.md).

## Workflow DAG

`src/workflow/Dag.h` defines a re-runnable, content-cached, cancellable graph. The important
property is that **node bodies are code, not documents**:

```cpp
using NodeInputs = std::map<std::string, std::string>;
using NodeFn = std::function<NodeResult(const NodeInputs&, const CancelToken&)>;
```

A `Node` carries `id`, `module`, `version`, `params`, `deps` and its `NodeFn`. The cache key is
`hash(module, version, params, {dep -> dep cache key})` and is therefore transitive:

- Re-running an unchanged graph is a full cache hit and completes instantly.
- Changing one node's `params` or `version` re-runs that node and everything downstream of it,
  and nothing else.

What is cached is **only the node's `output` string**, behind an `INodeCache` (in-memory, or
on disk under `%APPDATA%\BioCAD\cache`). That string is simultaneously the value handed to
dependents and the cached payload, which is why node functions must be deterministic in
`(inputs, params)`. Nothing else about a node - no side effect it performed, no file it wrote -
is captured by the cache, so a node that writes an artifact must key that artifact by content
itself.

`NodeStatus` transitions (`Pending`, `Running`, `Cached`, `Done`, `Failed`, `Cancelled`,
`Skipped`) fire an optional progress callback, which is what the live DAG view in the Workflows
panel renders.

## Job system and thread discipline

`src/workflow/JobSystem.h` is a deliberately tiny fixed-size thread pool: `std::thread`, a
queue, and `std::future`. The project hand-rolls it rather than taking Taskflow because the
value is in the DAG layer's caching and cancellation, not in the scheduler.

- Default width is `hardware_concurrency() - 1`, clamped to at least 1, leaving a core for the
  UI thread.
- `~JobSystem()` drains every queued task before joining, so no handed-out future is ever
  broken.
- `CancelToken` copies share one `std::atomic<bool>`. The UI flips it; long-running node
  functions poll `cancelled()` and bail cooperatively. Nothing is killed from outside.

The discipline is:

| Runs on the UI thread | Runs on workers |
| --- | --- |
| `wWinMain` message pump, ImGui frame, all panel drawing | DAG node functions |
| Reading module results that are already resolved | Receptor preparation and PDB fetch |
| Flipping a `CancelToken`, starting a run | Engine provisioning and download |
| SQLite reads/writes through `IRunStore` | Docking subprocess invocation |

Results come back as node `output` strings and through progress callbacks; the UI polls its
own copies rather than being called into from a worker.

## Persistence

`src/storage/RunStore.cpp` owns exactly **one** table, `run_history`, with six columns, opened
in WAL mode. Startup executes only `CREATE TABLE IF NOT EXISTS`. There is deliberately:

- **no migration mechanism** and no `PRAGMA user_version`. Adding a column is therefore a
  breaking change that needs a real migration story first, not a quiet schema edit;
- **no use of the name `runs`**. A legacy `biocad.db` left in `%APPDATA%` by an older Python
  backend has a `runs` table with an incompatible column layout, where
  `CREATE TABLE IF NOT EXISTS` would silently no-op and every `INSERT` would then fail on the
  column mismatch. `run_history` sidesteps that entirely.

Everything else on disk lives under `AppPaths` (`src/core/AppPaths.h`): `artifacts/`,
`runtime/`, `packs/`, `presets/`, `logs/`, `cache/`, `config.json`, `manifest.json`.

## The one-shot appdata migration

`AppPaths::migrateLegacyRoot()` (`src/core/AppPaths.cpp`) runs as the first statement of
`ensureLayout()`. If the BioCAD root does not exist but a sibling `StimLab` root does, it
renames the whole tree once, then:

- deletes `runtime/receptors` outright, because every cached PDBQT was written with the old
  `REMARK STIMLAB_BOX` marker and the parser now accepts only `REMARK BIOCAD_BOX` - there is
  deliberately no dual-accept path;
- renames `stimlab.db` to `biocad.db` and `logs/stimlab.log` to `logs/biocad.log`.

A failure logs a warning and returns false; the app then starts with a fresh root rather than
aborting. This is the entire compatibility story: no dual-read, no fallback resolver, no
permanent alias layer. The behaviour is pinned by
`TEST_CASE("AppPaths migrates a legacy StimLab root exactly once")` in `tests/test_core.cpp`.

The point of the migration is that provisioned docking engines are gigabytes; re-downloading
them because of a rename would be an unforced insult to the user.

## Feature gating

Three CMake options, each producing a `BIOCAD_HAVE_*` compile definition:

| Option | Macro | Effect |
| --- | --- | --- |
| `BIOCAD_ENABLE_SCIENCE` | `BIOCAD_HAVE_SCIENCE` | Links curl: the live Anthropic provider and the web tools |
| `BIOCAD_ENABLE_CUDA` | `BIOCAD_HAVE_CUDA` | Compiles the first-party CUDA docking backend (`sm_86` by default) |
| `BIOCAD_ENABLE_WEBVIEW2` | `BIOCAD_HAVE_WEBVIEW2` | Headless JS-rendered `web_fetch`; requires `BIOCAD_ENABLE_SCIENCE` |

The rule for every one of them: **the feature degrades cleanly when it is off.** A build
without science still runs the offline assistant, still navigates, still explains; it simply
has no network transport. `AnthropicProvider::transportAvailable()` returning false is a
first-class state the Settings panel renders, not an error.

`project(BioCAD VERSION x.y.z)` in the top-level `CMakeLists.txt` is the single source of the
version: `configure_file` generates `core/Version.h`, the PowerShell scripts parse the same
line via `scripts/version.ps1`, and CI asserts that `vcpkg.json` agrees.
