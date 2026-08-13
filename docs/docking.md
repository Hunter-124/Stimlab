# Docking

BioCAD does not implement a docking search. It locates, provisions and drives real
third-party engines (AutoDock Vina, smina, Vina-GPU), prepares receptors itself, parses the
engine output, and - when no engine can actually run - falls back to an explicitly labelled
descriptor estimate that is not a docked score.

Everything below is traceable to `src/modules/docking/`, `src/contracts/IDockingBackend.h`,
`src/modules/RealBackend.cpp` and `tests/test_docking.cpp`.

## Engine locate and provision order

### What is searched for

`docking::Engine` (`src/modules/docking/EngineLocator.h:23`) has four members: `Vina`, `Smina`,
`Obabel`, `VinaGpu`. Candidate file names per engine are listed most-specific-first in
`candidateNames()` (`EngineLocator.cpp:51-66`):

| Engine | Candidate names |
|---|---|
| Vina | `vina.exe`, `vina_1.2.5_win.exe`, `vina_1.2.5.exe`, `vina_1.2.3_win.exe`, `vina_win.exe`, `vina` |
| Smina | `smina.exe`, `smina_win.exe`, `smina.static.exe`, `smina` |
| Obabel | `obabel.exe`, `obabel` |
| VinaGpu | `Vina-GPU.exe`, `vina-gpu.exe`, `Vina-GPU` |

### Where it looks, in order

`locateEngine()` (`EngineLocator.cpp:160-196`):

1. `%APPDATA%\BioCAD\runtime\engines` (`enginesDir()`, `EngineLocator.cpp:129-131`). Vina-GPU is
   searched in its own subfolder `runtime\engines\vina-gpu` (`vinaGpuDir()`,
   `EngineLocator.cpp:133-135`) because it needs companion files beside the exe.
2. A loose-prefix scan of that same directory: any `vina*.exe` / `smina*.exe` / `obabel*.exe`
   dropped in by hand is accepted (`EngineLocator.cpp:172-190`). This scan is deliberately
   skipped for Vina-GPU, because it would otherwise match the `Vina-GPU-K.exe` kernel compiler
   that sits beside the docking binary.
3. `PATH`, split on `;` (`pathDirs()`, `EngineLocator.cpp:88-113`; `locateEngine`
   `EngineLocator.cpp:191-196`).

There is no other search location. Nothing is looked for in Program Files or the registry.

### Provisioning

`ensureVina(allowDownload)` (`EngineLocator.cpp:200-296`):

- If the engine is already located, it returns `fetched=true` immediately and touches no
  network.
- With `allowDownload=false` (the default, and what the app does at startup) it is a pure
  locate probe.
- With `allowDownload=true` it shells out to one `powershell -NoProfile -NonInteractive`
  `Invoke-WebRequest` (TLS 1.2 forced) for
  `https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_win.exe`
  (`EngineLocator.cpp:44-45`) into `runtime\engines\vina.exe`. Exit codes are meaningful:
  `0` ok, `2` transfer error, `3` SHA-256 mismatch, `4` implausibly small (an HTML error page).
  A file below `kVinaMinSizeBytes = 200000` (`EngineLocator.cpp:49`) is deleted so a later
  locate cannot trust it.

Integrity has two independent levels (`EngineLocator.cpp:38-49, 137-158`):

| Check | Source | Always on? |
|---|---|---|
| Size sanity | `kVinaSizeBytes = 1203712` published size; hard floor 200000 B | Yes |
| SHA-256 | first non-empty of env `BIOCAD_VINA_SHA256`, `runtime\engines\vina.sha256`, the compile-time pin `kVinaSha256` (empty by default) | Only when configured |

The runtime resolution order means a user can pin the hash without a rebuild.
`expectedVinaSha256()` is exposed so the UI can say whether the binary will be
cryptographically verified or only size-checked; `tests/test_docking.cpp:367` asserts that when
a pin exists it is well-formed 64-hex.

`ensureVinaGpu(allowDownload)` (`EngineLocator.cpp:298-390`) fetches only the ~3.7 MB of files
needed to run and compile the kernel from
`https://raw.githubusercontent.com/DeltaGroupNJUPT/Vina-GPU/main/` (`EngineLocator.cpp:70-71`;
file list `EngineLocator.cpp:75-89`) - never the 37 MB GUI in that repo. It then runs
`Vina-GPU-K.exe` once on the bundled `2bm2` example to compile a `Kernel2_Opt.bin` for *this*
machine's GPU, because OpenCL binaries are device specific. It reports `fetched=true` only when
both `Vina-GPU.exe` and a compiled `Kernel2_Opt.bin` are present; exit code `5` distinguishes
"downloaded but no usable OpenCL device".

`ensureObabel()` (`EngineLocator.cpp:393-405`) is **locate-only by design**: Open Babel ships as
an installer, not a drop-in exe, so it is never auto-downloaded. Drop `obabel.exe` into
`runtime\engines` to enable the higher-quality receptor prep path.

All of this runs off the UI thread. `docking::Provisioner` (`src/modules/docking/Provisioning.h:38-70`,
`Provisioning.cpp:53-92`) owns a worker thread, publishes a mutex-guarded status string plus
atomic counters (`vinaReady`, `obabelReady`, `receptorsReady/Total`), and runs three steps:
engines, then the headline receptors, then `writeRuntimeManifest()`. Nothing blocks startup; any
failure leaves a human-readable note and the app stays in the descriptor-estimate fallback.

### Manifest self-heal

`writeRuntimeManifest()` (`Provisioning.cpp:15-30`) records every located engine binary
(`engine:vina`, `engine:obabel`) and every `*.pdbqt` in the receptors directory
(`receptor:<stem>`) into `%APPDATA%\BioCAD\manifest.json` with its path, size and content hash.

`selfHealManifest()` (`Provisioning.cpp:32-37`) is the launch path: load the manifest, `heal()`,
then `verify()`. `Manifest::heal()` deletes every component whose file is present but whose size
or hash no longer matches, so the next provisioning run re-fetches it; missing files are simply
reported and left to provisioning (`src/core/Manifest.h:52-63`). The hash is a fast
non-cryptographic FNV-1a and is explicitly only a corruption/truncation check - download
*authenticity* is the pinned SHA-256 above, and the header says so
(`src/core/Manifest.h:10-12`).

## Receptor preparation

`ensureReceptor(target, allowDownload)` (`ReceptorPrep.cpp:396-500`) is the orchestrator:

1. **Cache wins.** `locatePreparedReceptor(id)` (`ReceptorPrep.cpp:371-394`) returns
   `runtime\receptors\<id>.pdbqt` if it exists and is non-empty. No network on this path, which
   is what the hot docking path uses (`RealBackend.cpp:403-409`).
2. **Fetch.** Otherwise, and only with `allowDownload`, the preset's PDB is downloaded from
   `https://files.rcsb.org/download/<PDBID>.pdb` (`ReceptorPrep.cpp:420`) via the same
   PowerShell `Invoke-WebRequest` helper (`ReceptorPrep.cpp:153-172`), cached beside the PDBQT.
   A response of 256 bytes or less is treated as an error page and deleted.
3. **Box derivation**, once, from the fetched text (`receptorBoxFromPdb`,
   `ReceptorPrep.cpp:316-368`): the centroid of the largest non-water HETATM residue with at
   least 8 atoms - the co-crystal ligand - or, failing that, the centroid of all protein atoms.
   `boxSource` records which.
4. **Prep, obabel path (preferred).** If `obabel.exe` was located, `obabel <pdb> -xr -O <pdbqt>`
   is run with a 120 s timeout and the box remark is prepended to its output
   (`ReceptorPrep.cpp:446-470`). This adds polar hydrogens.
5. **Prep, built-in path (always available).** `pdbToRigidReceptor()`
   (`ReceptorPrep.cpp:195-313`).

### What the built-in writer does

- First model only (`MODEL`/`ENDMDL` handling, `ReceptorPrep.cpp:211-220`); one conformation
  only (altLoc blank, `A` or `1`, `ReceptorPrep.cpp:235`).
- Waters and additives dropped (`isWater`, `ReceptorPrep.cpp:64-67`), and so is the co-crystal
  ligand: an empty `keepHetero` set yields a clean apo-like pocket. A named cofactor can be kept
  (`keepHetero = {"FAD"}`), covered by `tests/test_docking.cpp:333`.
- **Modified amino acids are kept, not stripped** (`aminoAcidParent`, `ReceptorPrep.cpp:76-89`):
  MSE, FME, MHO, SEP, TPO, PTR, CSO, CSD, OCS, CME, KCX, MLY, M3L, ALY, HYP, HIC, PCA. These are
  recorded as HETATM but are covalently part of the chain; dropping them tears a gap in the
  backbone and silently deletes pocket residues. They are emitted as `ATOM` under their standard
  parent name; MSE's `SE` becomes methionine's `SD`, and selenium (Z=34) is re-typed as sulfur
  because AutoDock4 has no Se type (`ReceptorPrep.cpp:255-268`). `tests/test_docking.cpp:281`
  guards this.
- Hydrogens are dropped: the built-in path is heavy-atom (`ReceptorPrep.cpp:267`).
- Elements are inferred from PDB columns 77-78 when valid, else from the first alphabetic
  character of the atom name, with a carbon-like default that is never 0 or NaN
  (`inferElement`, `ReceptorPrep.cpp:91-118`).
- Aromatic ring atoms of PHE/TYR/TRP/HIS (and HID/HIE/HIP) are recognised from an explicit
  residue+atom-name table (`isAromaticRingAtom`, `ReceptorPrep.cpp:38-61`) so their carbons get
  the AutoDock4 `A` type and the imidazole/indole nitrogens get `NA`.
- Output is a fixed-width PDBQT `ATOM` line: standard PDB columns 1-66, then the PDBQT trailing
  charge (`%10.3f`) and a right-justified 2-character atom type (`ReceptorPrep.cpp:283-286`).
  There is no torsion tree - a receptor is rigid by definition.

### Atom typing

Shared with the ligand writer: `autodockAtomType(z, aromatic, polarH)`
(`PdbqtWriter.cpp:115-130`) - `HD` for polar hydrogen, `A` for aromatic carbon else `C`, `NA`
for aromatic nitrogen else `N`, `OA`, `SA`/`S`, `F`, `P`, `Cl`, `Br`, `I`. The ligand writer
emits only polar hydrogens (bonded to N/O/S), matching AutoDock's united-atom convention
(`isPolarHydrogen`, `PdbqtWriter.cpp:66-76`).

### Why approximate charges are acceptable

Receptor charges come from `approxReceptorCharge()` (`ReceptorPrep.cpp:120-138`): a flat
per-element table (O -0.30, aromatic N -0.20, N -0.30, S -0.10, aromatic C 0.00, C 0.08,
P 0.30, polar H 0.20). Ligand charges come from `approxPartialCharge()`, an
electronegativity-difference smear over the bond list, not a converged Gasteiger iteration
(`PdbqtWriter.h:36-39`, electronegativity table `PdbqtWriter.cpp:18-33`).

This is honest rather than lazy: **AutoDock Vina's empirical scoring function is
partial-charge-independent**, so the charge column is informational for AutoDock4-style
consumers and does not enter the Vina score. The generated receptor says so in its own header
remarks (`ReceptorPrep.cpp:296-298`):

```text
REMARK  BioCAD rigid receptor (built-in heavy-atom prep; AutoDock4 types)
REMARK  binding-affinity / target-engagement prediction only; not a synthesis artifact
REMARK  Vina scoring is partial-charge-independent; charges are approximate
```

A converged protein PEOE solve is out of scope without RDKit, and pretending otherwise would be
worse than stating the approximation.

## The `REMARK BIOCAD_BOX` marker

**Written by** `boxRemark()` (`ReceptorPrep.cpp:182-187`), which formats exactly:

```text
REMARK BIOCAD_BOX %.3f %.3f %.3f
```

and is prepended to the prepared PDBQT on both the obabel path (`ReceptorPrep.cpp:461-466`) and
the built-in path (`ReceptorPrep.cpp:485-488`).

**Parsed by** `locatePreparedReceptor()` (`ReceptorPrep.cpp:379-388`), which scans the leading
remark block for the literal prefix `REMARK BIOCAD_BOX`, reads three doubles, and stops at the
first `ATOM`/`HETATM` line.

**Why it exists.** The pack preset carries a literature-style box centre, which is frequently in
a different coordinate frame from the structure that was actually downloaded. `RealDocking::dockDetailed()`
overrides the preset centre with the recovered structure-frame centre when one is present
(`RealBackend.cpp:397-409`), which is what makes the search box actually overlap the receptor.

**Why there is no dual-accept for the old `REMARK STIMLAB_BOX`.** A dual-accept parser would be
dead code by construction: `AppPaths::migrateLegacyRoot()` deletes the entire cached receptor
tree during the one-shot rename migration -
`std::filesystem::remove_all(root_ / "runtime" / "receptors")` with the comment "Cached receptors
carry the old REMARK STIMLAB_BOX marker and must be re-prepared" (`src/core/AppPaths.cpp:72-73`).
After migration no file on disk can carry the old marker, so accepting it would only add a
second parse path with no reachable input. One token, one writer, one reader.

## The convergent search

A single Vina run is a stochastic Monte Carlo search: run it twice with different seeds and you
get different numbers. `runVinaSearch()` (`Backends.cpp:81-201`), shared by the CPU and GPU Vina
backends, therefore does not report one run.

- The ligand is written as a **torsionally flexible** PDBQT (`writeFlexiblePdbqt`,
  `PdbqtWriter.h:60-68`) so the engine may bend and rotate it. A rigid PDBQT is also prepared;
  if the flexible torsion tree is rejected on the first run, the search silently retries rigid
  (`Backends.cpp:161-165`) so a flex-writer bug can never lose a dock that would otherwise work.
- Up to `kMaxRuns = 5` **independent** runs, each with a fresh `--seed` (`Backends.cpp:119`,
  seed `i+1`, `Backends.cpp:160`).
- On CPU, escalating `--exhaustiveness` on the schedule `{16, 24, 32, 32, 32}`
  (`Backends.cpp:148`). Vina-GPU has no `--exhaustiveness`; thoroughness comes from
  `--thread 8000` lanes plus independent seeds, and its box edges are clamped to < 30 A per the
  engine's own limitation (`Backends.cpp:113,125-128`).
- Stop early when two consecutive runs agree within `kTol = 0.2` kcal/mol; otherwise stop at the
  run cap or once ~150 s of wall clock has elapsed, but always do at least two runs
  (`Backends.cpp:145-175`).
- Each run's output file is deleted before the run so a stale output can never be misparsed as
  this run's poses (`Backends.cpp:137`). Each subprocess has a 300 s hard timeout and is
  terminated rather than allowed to wedge the worker (`runProcess`, `Backends.cpp:36-62`).

### What the four UI fields mean

`DockJobResult` (`src/contracts/IDockingBackend.h:59-79`), rendered in the docking panel
(`src/ui/Panels.cpp:1047-1105`):

| Field | Meaning to a reader |
|---|---|
| `searchRuns` | How many independent searches (fresh seed) were aggregated. The reported pose is the best over all of them. UI row "Independent runs". |
| `affinitySpread` | max minus min of the per-run best affinity, in kcal/mol. This is the search's own reproducibility, **not** an experimental error bar - it says nothing about whether the scoring function is right. It is passed as the `error` term of the rendered `Quantity` (`Panels.cpp:1058`). |
| `converged` | Two consecutive runs agreed within 0.2 kcal/mol. The panel shows CONFIDENCE HIGH when true, MODERATE when the run/time budget was hit instead (`Panels.cpp:1047-1049`). "Converged" means the *search* stabilised, nothing more. |
| `torsions` | Active rotatable bonds the ligand was docked with (TORSDOF). `0` means the flexible tree was rejected and the rigid fallback ran - the ligand was docked as a frozen conformer, which is a real limitation of that result. |

The engine's own log line records all of it verbatim, e.g.
`AutoDock Vina: 3 run(s), best -9.31 kcal/mol, spread 0.14; converged; 4 torsions (flexible).`
(`Backends.cpp:191-196`).

### Engine selection order

`RealDocking::realEngines()` (`RealBackend.cpp:460-486`) builds the try-order from the
process-wide `ComputeMode`:

| Mode | Order |
|---|---|
| `Auto` | Vina, smina, Vina-GPU, CUDA (CUDA only under `BIOCAD_HAVE_CUDA`) |
| `Cpu` | Vina, smina |
| `Gpu` | Vina-GPU, then CUDA |

The first backend that is `available()` and returns a non-empty result with
`provenance == Model` wins (`RealBackend.cpp:420-424`); otherwise the descriptor estimate
carries the result. The first-party CUDA backend is honest about being a **rigid-body** grid
search that does not vary torsions, so flexible CPU Vina stays the accurate path
(`CudaBackend.h:8-12`).

## What a Vina score is and is not

Read this before quoting a number out of this app.

- It **is** the value of an empirical scoring function, in kcal/mol, for one pose of one ligand
  in one rigid receptor inside one box.
- Its **reported standard error is 2.85 kcal/mol** (Trott & Olson 2010, PMC3041641), cited in
  `src/chem/AdmetModel.h:170-173`. At 298 K that is roughly a factor of 123 in Kd. Any
  nanomolar figure derived from a docked score is manufactured precision.
- It is **`Provenance::Model`**, not `Measured` and not `Predicted`
  (`IDockingBackend.h:59-63`). A pose is a constructed artefact. Nothing was measured, and no
  published affinity model with a benchmark error ran.
- It **ranks poses far better than it predicts affinity**. Comparing poses of one ligand in one
  box is the operation it was fitted for. Comparing scores across different targets, different
  receptor preparations, different box volumes or different rotatable-bond counts is not.
- The app **refuses to convert it into a Kd**. `kdFromDeltaG` / `deltaGFromKd` exist in
  `src/chem/AdmetModel.h:179-186` with the constant `R = 1.987204259e-3 kcal/(mol*K)`, and the
  header states they are "DELIBERATELY NOT WIRED TO DOCKING SCORES" and take a *measured* dG.
  The docking panel shows kcal/mol with the `model` badge and never a nM affinity
  (`Panels.cpp:1053-1061`).
- `affinitySpread` is search noise, not accuracy. A spread of 0.05 kcal/mol over five runs and a
  standard error of 2.85 kcal/mol are entirely compatible statements.

## The descriptor-estimate fallback

`EstimateBackend::dock()` (`Backends.cpp:204-231`) fires whenever no real engine produced a
dock. In practice that is: no engine binary located; a binary present but no prepared receptor
for the target; the engine ran but parsed no poses; or the ligand SMILES failed to parse
(`RealBackend.cpp:411-427`, `Backends.cpp:238-262`).

It computes six ranked pseudo-poses from `-(5.0 + logP*0.8 + MW*0.004)`, stepped by 0.35, using
the in-house `chem::crippenLogP` and `chem::molecularWeight`, and attaches the embedded
conformer purely so the 3D viewer has geometry. It is **not** a dock and it did not consult the
receptor at all.

It is labelled as such everywhere:

| Surface | Label |
|---|---|
| `DockJobResult::engine` | the literal string `descriptor-estimate` (`Backends.cpp:207`), asserted by `tests/test_docking.cpp:237` |
| `DockJobResult::provenance` | `Provenance::Heuristic` (`Backends.cpp:208`), so `fromEngine()` is false |
| `DockJobResult::log` | "No docking engine provisioned; affinity is a structure-descriptor estimate (logP/MW), NOT a docked score..." (`Backends.cpp:225-229`) |
| Docking panel | `drawQuantity` with **no unit** and source "descriptor estimate - rank ordering only", rendered in the Heuristic colour (`Panels.cpp:1056-1061`) |
| Legacy summary string | "Estimated affinity ... (descriptor-estimate - structure-descriptor model, not a docked score)." (`RealBackend.cpp:438-441`) |
| Headless selftest | `--selftest-dock` exits `0` for a real engine dock and `2` for the descriptor fallback (`src/app/WinMain.cpp:54-55`) |

The unit-lessness is enforced, not merely conventional: `makeQuantity` throws when a
`Provenance::Heuristic` value carries a unit, and `tests/test_provenance.cpp:21-27` asserts that
a heuristic `-7.4 kcal/mol` is rejected while a `Provenance::Model` one is accepted. A heuristic
score with a kcal/mol label is unrepresentable in this codebase.
