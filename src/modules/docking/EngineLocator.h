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
enum class Engine { Vina, Smina, Obabel };

// The directory engine binaries live in: %APPDATA%/StimLab/runtime/engines.
std::filesystem::path enginesDir();

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
// returns fetched=true immediately. Otherwise tries ONE PowerShell download +
// SHA-256 verify; on any failure returns fetched=false with a note. `allowDownload`
// must be true to touch the network at all (default false = locate-only).
ProvisionResult ensureVina(bool allowDownload = false);

}  // namespace stimlab::docking
