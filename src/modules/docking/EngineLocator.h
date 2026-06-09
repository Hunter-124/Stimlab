// modules/docking/EngineLocator.h - find (or best-effort provision) docking exes.
//
// Locates AutoDock Vina (vina.exe / vina_*.exe), smina (smina.exe / smina_*.exe)
// and Open Babel (obabel.exe) under %APPDATA%/StimLab/runtime/engines and then on
// PATH. Also exposes an OPTIONAL, non-blocking provisioner that downloads the
// official Vina Windows binary via a PowerShell Invoke-WebRequest + Get-FileHash
// SHA-256 verification into the engines dir. Provisioning never blocks startup and
// degrades silently to "unavailable" on any network/AV/hash failure (no libcurl).
//
// SAFETY SCOPE: this only locates third-party docking BINARIES used for a binding-
// affinity prediction. It contains no synthesis content.
#pragma once

#include <filesystem>
#include <optional>
#include <string>

namespace stimlab::docking {

// Which engine to look for / provision.
//   Vina/Smina - CPU AutoDock Vina / smina (flexible-ligand, the accurate path).
//   Obabel     - Open Babel (receptor protonation helper).
//   VinaGpu    - Vina-GPU: an OpenCL GPU re-implementation of Vina (flexible-ligand,
//                NVIDIA/AMD via the driver ICD). Lives in its OWN subfolder because it
//                needs companion files beside the exe (a GPU-matched Kernel2_Opt.bin +
//                the OpenCL kernel source the -K compiler reads).
enum class Engine { Vina, Smina, Obabel, VinaGpu };

// User compute-mode preference for docking-engine selection:
//   Auto - best available: the accurate CPU engine (Vina) first, the GPU engine as a
//          fallback when no CPU engine is present, then the labeled estimate.
//   Gpu  - force the CUDA GPU engine (then fall back to the labeled estimate).
//   Cpu  - CPU engines only (skip the GPU engine).
// Process-wide (a single global user setting); the docking module reads it when
// ordering engines. Defaults to Auto. Thread-safe.
enum class ComputeMode { Auto, Gpu, Cpu };
void                      setComputeMode(ComputeMode m);
[[nodiscard]] ComputeMode computeMode();

// The directory engine binaries live in: %APPDATA%/StimLab/runtime/engines.
std::filesystem::path enginesDir();

// The Vina-GPU OpenCL engine's own subfolder: %APPDATA%/StimLab/runtime/engines/vina-gpu.
// It needs companion files beside the exe (Kernel2_Opt.bin + the OpenCL kernel source),
// so it gets a folder of its own rather than sitting loose in engines/.
std::filesystem::path vinaGpuDir();

// Locate an engine binary. Search order: runtime/engines (incl. vina_1.2.5_win.exe
// style names), then PATH. Returns the absolute path, or nullopt if not found.
std::optional<std::filesystem::path> locateEngine(Engine e);

// Convenience: is this engine available right now?
bool engineAvailable(Engine e);

// Result of a provisioning attempt (best-effort, never throws).
struct ProvisionResult {
    bool        fetched = false;   // true iff a verified binary now exists on disk
    std::string path;              // absolute path to the binary (if fetched)
    std::string note;             // human-readable status / failure reason
};

// Best-effort: ensure vina.exe exists under runtime/engines. If already present,
// returns fetched=true immediately. Otherwise tries ONE PowerShell download with a
// size sanity-check + optional SHA-256 verify; on any failure returns fetched=false
// with a note. `allowDownload` must be true to touch the network (default false =
// locate-only). The expected SHA-256 is resolved at runtime (first non-empty of:
// env STIMLAB_VINA_SHA256, a `vina.sha256` file in the engines dir, the compile-time
// pin) so a user can pin integrity WITHOUT a rebuild; the published file size is
// always checked so an HTML error page can never masquerade as the binary.
ProvisionResult ensureVina(bool allowDownload = false);

// Best-effort: ensure the Vina-GPU (OpenCL) engine is provisioned under
// runtime/engines/vina-gpu. Downloads only the pieces needed (the prebuilt Vina-GPU.exe
// + Vina-GPU-K.exe + the small OpenCL kernel source + one example complex -- NOT the
// 37 MB GUI in the upstream repo), then runs Vina-GPU-K.exe ONCE on the example to
// compile a `Kernel2_Opt.bin` matched to THIS machine's GPU (OpenCL binaries are device
// specific). `allowDownload` must be true to touch the network (default false =
// locate-only). Reports fetched=true only when BOTH Vina-GPU.exe and a compiled
// Kernel2_Opt.bin are present. Degrades to fetched=false (never throws) on any
// download / no-OpenCL-device failure, so the GPU engine simply stays unavailable.
ProvisionResult ensureVinaGpu(bool allowDownload = false);

// Best-effort: locate obabel.exe (runtime/engines or PATH). Receptor prep prefers
// it when present (adds polar H), else uses the built-in heavy-atom writer. This is
// LOCATE-ONLY by design: OpenBabel ships as an installer, not a drop-in exe, so we
// never auto-download it - drop obabel.exe into runtime/engines to enable it.
ProvisionResult ensureObabel();

// The expected SHA-256 for vina.exe actually in force right now (env / file / pin),
// or empty if none is configured. Exposed so the UI can show whether the provisioned
// binary will be cryptographically verified or only size-checked.
std::string expectedVinaSha256();

}  // namespace stimlab::docking
