# StimLab — Session Card: "do everything" (remaining post-v1 future work)

> Paste-to-start brief for a FRESH session. **StimLab** is a native Windows (C++20, DirectX 11 + Dear ImGui)
> CNS-stimulant computational pharmacology suite with an **in-house cheminformatics engine (NO RDKit** —
> rdkit/openbabel/indigo/inchi are NOT vcpkg ports here; do not try to link them). Repo:
> `C:\Users\nigga\Desktop\Stimulant-Laboratory` (git, branch `master`).

## 0. READ FIRST (in order)
`handoff.md` (STATUS banner + §2 safety boundary + §4 verified env + §7 gotchas), `docs/PHASE_E.md` (the
"Post-v1 enhancements — ALL FIVE LANDED" section), and the auto-memory: `stimlab-post-v1`,
`stimlab-build-and-deps`, `stimlab-real-docking`, `stimlab-v1-complete`.

## 1. SAFETY BOUNDARY — NON-NEGOTIABLE
Predict what a compound IS and DOES (structure, pharmacology, ADMET, stability, similarity, legal-analog,
docking = target engagement). **NEVER** produce synthesis routes, reaction conditions, precursor selection, or
manufacturability content. Docking "binding affinity" is a pharmacology signal, never a make-it signal. The
agent system prompt + tool set must keep enforcing this.

## 2. CURRENT STATE (all green, committed on `master`)
v1 feature-complete + all 5 optional post-v1 tracks landed (newest first):
`233b6b9` docs · `3f3af7a` Track 1 CUDA GPU docking · `25fe293` Track 2 WebView2 web_fetch ·
`a687b59` Track 5 CI+signing · `2a6ef65` Track 4 science-static exe · `efd2c49` Track 3 receptors-on-demand.
ctest **69/69** on windows + windows-static; science, science-static, windows-cuda all green. Tree clean.

## 3. BUILD / VERIFY (hard-won — obey exactly)
- Build ONLY via the **PowerShell tool** (NOT the Bash tool — it mangles the `.ps1` path; a `127` exit means it
  ran the script as bash). `cl`/`nvcc` are NOT on PATH; the wrappers enter vcvars + pin `VCPKG_ROOT=C:\Users\nigga\vcpkg`.
- `.\scripts\dev-build.ps1 <preset> -Test`  |  `.\build.ps1 [-Release] [-Science]`  |  `.\scripts\ci.ps1 [-Science]`.
- Presets: `windows` (fast, curl-free) · `windows-science` (+curl[ssl]=Schannel, +webview2) · `windows-static`
  (static single exe) · `windows-science-static` (static + live agent/web tools) · `windows-cuda` (lean +
  first-party CUDA GPU docking; nvcc 13.3, sm_86). Keep ctest GREEN; commit per item.
- Headless: `StimLab.exe --selftest-dock [--smiles S] [--target T] [--compute auto|gpu|cpu]` (exit 0 = real
  engine dock; 2 = labeled estimate). Hidden live tests (run explicitly): `stimlab_tests.exe "[live]"`
  (web_search), `"[webview2]"` (JS render). GUI: `.\scripts\capture-window.ps1 -Exe <exe> -Out <ABSOLUTE>.png
  -Panel <Name>` (+ `STIMLAB_TARGET` env). Commit trailer:
  `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

## 4. ENVIRONMENT (verified)
RTX 3080 Ti (sm_86), driver 596.36. **CUDA Toolkit is now FULL: nvcc 13.3.33** at
`C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.3\bin` (v13.2 stays runtime-only). MSVC v143 (cl 19.44),
CMake 4.x, Ninja, vcpkg at `C:\Users\nigga\vcpkg`. **Remote: `Hunter-124/Stimulant-Laboratory` (private)** via the
cached GCM credential. WSL2 not used (user declined; Vina-GPU covers native GPU docking).

## 5. THE WORK — "do everything"  (one commit per item, green trunk; do A→B first, they're self-contained)

### A. Static cudart → self-contained GPU exe  [self-contained; do first]
> ✅ **DONE** (`913e09d`): `windows-cuda-static` preset + `CUDA::cudart_static` (triplet-selected). 3.21 MB
> exe, no cudart DLL import, clean-dir `--compute gpu` real dock, ctest 70/70.

`windows-cuda` links DYNAMIC `CUDA::cudart` (needs `cudart64_*.dll` on PATH). Add a `windows-cuda-static` preset
(x64-windows-static + `/MT` + `STIMLAB_ENABLE_CUDA`, point `CMAKE_CUDA_COMPILER` at nvcc 13.3, arch 86) and link
**`CUDA::cudart_static`** in `src/modules/CMakeLists.txt` under `STIMLAB_ENABLE_CUDA`. Verify the exe runs
standalone from a clean dir with NO cudart DLL beside it, and `--selftest-dock --compute gpu` gives a real GPU
dock. Watch the `/MT` ↔ CUDA host-runtime match.

### B. Vina-GPU (OpenCL) docking backend  [self-contained; genuinely FLEXIBLE GPU docking — high value]
> ✅ **DONE** (`3145966` + UI `8832466`): `Engine::VinaGpu` + `ensureVinaGpu()` (fetches ~3.7 MB, compiles a
> GPU-matched kernel) + `VinaGpuBackend` (reuses `parseVinaPdbqt`, box <30, `--seed 1`) + `realEngines()` Gpu
> ordering + Settings provision button. Verified amphetamine→DAT real dock −4.80 kcal/mol. NOTE: rigid-ligand
> like the CPU path (the in-house PDBQT writer is rigid); the win over the CUDA grid is the real Vina MC+BFGS
> search, not torsions. Truly-flexible ligand docking stays Track F.

Add a 2nd GPU engine behind the frozen `IDockingBackend` (`src/contracts/IDockingBackend.h`). The
**DeltaGroupNJUPT/Vina-GPU** repo ships PREBUILT Windows `Vina-GPU.exe` + `Vina-GPU-K.exe` (OpenCL on NVIDIA via
the driver ICD — no CUDA build), taking Vina-style PDBQT receptor/ligand + `center_x/y/z` + `size_x/y/z` — the
exact existing `DockBox`+PDBQT contract. Mirror `VinaBackend` (`src/modules/docking/Backends.cpp`): add
`Engine::VinaGpu` + candidate names to `EngineLocator`, an `ensureVinaGpu()` provisioner (FOLDER-based: the exe
+ `Kernel2_Opt.bin`; run `Vina-GPU-K.exe` once to compile the kernel), a `VinaGpuBackend` (write `config.txt`,
run, **reuse `parseVinaPdbqt`**), and add it to `RealDocking::realEngines()` ordered by compute mode. Label
"GPU (OpenCL, Vina-GPU)". This is FLEXIBLE docking (more accurate than the rigid first-party CUDA engine).
Verify a real dock vs DAT.

### C. WSL2 → Uni-Dock / AutoDock-GPU subprocess  [PRECONDITION: WSL2 + a distro — confirm with user]
> ⏭ **SUPERSEDED BY B — NO WSL** (user declined WSL). Vina-GPU (Track B) delivers the real Vina GPU search
> natively on Windows on this same NVIDIA GPU, so the WSL2/Uni-Dock path is intentionally not pursued; a
> native-Windows AutoDock-GPU CUDA source build stays infeasible. (Torsional flexibility = Track F.)

Real CUDA *flexible* GPU docking via Linux binaries. A `WslDockBackend : IDockingBackend` that shells
`wsl.exe -e <linux binary>` to **Uni-Dock** (CUDA Vina re-impl; same PDBQT + center/size; compute ≥7.0 ✓) or
AutoDock-GPU `adgpu_*_cuda12`, with `/mnt/c` ↔ `\\wsl$` path translation. WSL2 GPU passthrough is real (do NOT
install a Linux GPU driver — the Windows driver is stubbed as `libcuda.so`). Confirm `wsl --status` first; offer
`wsl --install` if absent. (Do NOT attempt the AutoDock-GPU native-Windows CUDA *source* build — assessed
infeasible: Unix Makefile + CUDA 13 drops pre-sm_75.)

### D. Git remote + live CI run  [PRECONDITION: a GitHub repo — `gh` or user-provided remote]
> ✅ **DONE — CI GREEN**: private repo `Hunter-124/Stimulant-Laboratory` created + `master` pushed (cached GCM
> credential). CI fix: removed `vcpkg.json` `builtin-baseline` (a LOCAL-only vcpkg commit the runner couldn't
> resolve → `vcpkg install` failed). Run #2 SUCCESS (22 min): windows + windows-static build+ctest, release zip
> packaged + uploaded as artifact `StimLab-win-x64`, code-sign skipped (no cert).

The repo has NO remote, so `.github/workflows/ci.yml` is written-but-unrun. Create/connect a GitHub repo
(`gh repo create` or a user remote), push `master`, confirm the Actions run goes green on windows-latest (it
uses the preinstalled vcpkg + x-gha binary cache; sanity-check the `vcpkg.json` builtin-baseline against the
runner). The job uploads the release zip artifact.

### E. Code-signing  [PRECONDITION: a real Authenticode cert — USER provides]
> ⏭ **DECIDED — INTENTIONALLY UNSIGNED** (user: educational/personal use, doesn't want a cert). `sign.ps1`
> stays a clean exit-0 no-op; the packaged README.txt now states the build is unsigned (SmartScreen expected).
> No cert acquired. Nothing further to do unless a cert appears.

`scripts/sign.ps1` is ready (signtool via `STIMLAB_SIGN_PFX`/`_PASSWORD` or `STIMLAB_SIGN_SHA1`; clean no-op
without a cert). Ask the user for a cert; if none, leave it unsigned (SmartScreen warning expected) + documented.
Do NOT acquire a cert on the user's behalf.

### F. (Stretch) Flexible/torsional first-party CUDA docking
`CudaBackend`/`CudaScore.cu` is RIGID-body. A rotatable-bond torsional search, or a CUDA rescoring layer over
Vina poses, would raise fidelity. Large — only after A–E.

## 6. USER-DEPENDENT BLOCKERS (raise in the first message)
C needs WSL2 + a distro; D needs a GitHub repo/remote; E needs a code-signing cert. A, B, F are self-contained —
start there while the user resolves C/D/E.

## 7. KEY GOTCHAS
- nvcc not on the vcvars PATH → the windows-cuda preset sets `CMAKE_CUDA_COMPILER` explicitly +
  `CMAKE_CUDA_ARCHITECTURES=86` + `CMAKE_CUDA_FLAGS=-allow-unsupported-compiler`; MSVC C++ flags are scoped to
  `$<COMPILE_LANGUAGE:CXX>` so nvcc isn't handed `/permissive-`.
- WebView2 static loader (`WebView2LoaderStatic.lib`) keeps the exe DLL-free; the rendered path needs the
  machine's Evergreen WebView2 Runtime, else `web_fetch` degrades to curl.
- `curl[ssl]` = Schannel here (no OpenSSL). New vcpkg adds: `eigen3`, `webview2`, `wil`.
- A real engine dock must NEVER run on the UI thread (use the off-thread `Provisioner`/`dockFor` pattern).
- Late-session GUI screenshots can come back blank if the display sleeps — that is ENVIRONMENTAL (re-capture
  or verify headlessly), not a code regression.
