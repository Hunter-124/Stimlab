#include "modules/docking/EngineLocator.h"

#include <array>
#include <atomic>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <string>
#include <system_error>
#include <vector>

#include "core/AppPaths.h"

#if defined(_WIN32)
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#endif

namespace biocad::docking {

// Process-wide compute-mode preference (a single global user setting). Thread-safe.
namespace {
std::atomic<int> g_computeMode{static_cast<int>(ComputeMode::Auto)};
}  // namespace
void setComputeMode(ComputeMode m) {
    g_computeMode.store(static_cast<int>(m), std::memory_order_relaxed);
}
ComputeMode computeMode() {
    return static_cast<ComputeMode>(g_computeMode.load(std::memory_order_relaxed));
}

namespace {

namespace fs = std::filesystem;

// Official AutoDock Vina 1.2.5 Windows release asset. GitHub publishes no per-asset
// digest for this release, so kVinaSha256 is left empty by default; integrity can be
// pinned at runtime (env BIOCAD_VINA_SHA256 or a `vina.sha256` file in the engines
// dir) WITHOUT a rebuild. The published byte size IS known and is always enforced as
// a sanity check, so an HTML error page / truncated transfer is rejected even when no
// hash is configured. Provisioning stays best-effort and non-blocking either way.
constexpr const char* kVinaUrl =
    "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_win.exe";
constexpr const char* kVinaSha256 = "";          // optional compile-time pin
constexpr long long   kVinaSizeBytes = 1203712;  // published size of vina_1.2.5_win.exe
constexpr long long   kVinaMinSizeBytes = 200000;  // reject anything implausibly small

// Candidate file names per engine (most specific first).
std::vector<std::string> candidateNames(Engine e) {
    switch (e) {
        case Engine::Vina:
            return {"vina.exe", "vina_1.2.5_win.exe", "vina_1.2.5.exe",
                    "vina_1.2.3_win.exe", "vina_win.exe", "vina"};
        case Engine::Smina:
            return {"smina.exe", "smina_win.exe", "smina.static.exe", "smina"};
        case Engine::Obabel:
            return {"obabel.exe", "obabel"};
        case Engine::VinaGpu:
            return {"Vina-GPU.exe", "vina-gpu.exe", "Vina-GPU"};
    }
    return {};
}

// Vina-GPU (OpenCL) provisioning source: the prebuilt Windows binaries + OpenCL kernel
// source committed in the DeltaGroupNJUPT/Vina-GPU repo (no GitHub release; raw files).
// We fetch ONLY the pieces needed to run + compile the kernel (~3.7 MB), never the 37 MB
// GUI that dominates the repo. The kernel binary is compiled locally for this GPU.
constexpr const char* kVinaGpuRawBase =
    "https://raw.githubusercontent.com/DeltaGroupNJUPT/Vina-GPU/main/";
constexpr long long kVinaGpuMinExeBytes = 200000;  // the real exe is ~1.14 MB
// Files to fetch, in repo-relative form (forward slashes; rewritten to the local tree).
const std::vector<std::string>& vinaGpuFiles() {
    static const std::vector<std::string> files = {
        "Vina-GPU.exe", "Vina-GPU-K.exe",
        "OpenCL/inc/commonMacros.h", "OpenCL/inc/kernel2.h", "OpenCL/inc/kernel_string.h",
        "OpenCL/inc/wrapcl.h", "OpenCL/src/kernels/code_head.cpp",
        "OpenCL/src/kernels/kernel1.cl", "OpenCL/src/kernels/kernel2.cl",
        "OpenCL/src/kernels/matrix.cpp", "OpenCL/src/kernels/mutate_conf.cpp",
        "OpenCL/src/kernels/quasi_newton.cpp", "OpenCL/src/wrapcl.cpp",
        "input_file_example/2bm2_config.txt", "input_file_example/2bm2_ligand.pdbqt",
        "input_file_example/2bm2_protein.pdbqt"};
    return files;
}

// Split a PATH-style environment string on ';' (Windows) into directories.
std::vector<fs::path> pathDirs() {
    std::vector<fs::path> dirs;
#if defined(_WIN32)
    if (const char* path = std::getenv("PATH"); path && *path) {
        std::string cur;
        for (const char* p = path; ; ++p) {
            if (*p == ';' || *p == '\0') {
                if (!cur.empty()) {
                    std::error_code ec;
                    fs::path d(cur);
                    if (fs::exists(d, ec)) dirs.push_back(d);
                }
                cur.clear();
                if (*p == '\0') break;
            } else {
                cur.push_back(*p);
            }
        }
    }
#endif
    return dirs;
}

std::string trimHash(std::string s) {
    // Keep only hex digits (Get-FileHash emits uppercase; a file may have a newline).
    std::string out;
    for (char c : s) if (std::isxdigit(static_cast<unsigned char>(c))) out.push_back(c);
    return out;
}

std::string readFirstLine(const fs::path& p) {
    std::error_code ec;
    if (!fs::exists(p, ec)) return {};
    std::ifstream in(p);
    std::string line;
    std::getline(in, line);
    return line;
}

}  // namespace

fs::path enginesDir() {
    return AppPaths::instance().runtime() / "engines";
}

fs::path vinaGpuDir() {
    return enginesDir() / "vina-gpu";
}

std::string expectedVinaSha256() {
    // Runtime override precedence: env var, then a vina.sha256 file beside the binary,
    // then the compile-time pin. All normalised to bare uppercase hex.
    if (const char* env = std::getenv("BIOCAD_VINA_SHA256"); env && *env) {
        std::string h = trimHash(env);
        if (h.size() == 64) {
            for (auto& c : h) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            return h;
        }
    }
    std::string fileHash = trimHash(readFirstLine(enginesDir() / "vina.sha256"));
    if (fileHash.size() == 64) {
        for (auto& c : fileHash) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        return fileHash;
    }
    std::string pin = trimHash(kVinaSha256);
    if (pin.size() == 64) {
        for (auto& c : pin) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        return pin;
    }
    return {};
}

std::optional<fs::path> locateEngine(Engine e) {
    std::error_code ec;
    const auto names = candidateNames(e);

    // 1) runtime/engines (preferred, app-provisioned). Vina-GPU lives in its own
    // subfolder (it carries companion files) so it is searched there.
    const fs::path dir = (e == Engine::VinaGpu) ? vinaGpuDir() : enginesDir();
    for (const auto& n : names) {
        const fs::path cand = dir / n;
        if (fs::exists(cand, ec) && fs::is_regular_file(cand, ec)) return cand;
    }
    // Also accept any vina_*/smina_* file dropped into engines without an exact name.
    // Skipped for Vina-GPU: its name is fixed and a loose-prefix scan would wrongly
    // match the Vina-GPU-K.exe kernel compiler that sits beside the docking exe.
    if (e != Engine::VinaGpu && fs::exists(dir, ec) && fs::is_directory(dir, ec)) {
        const std::string stem = (e == Engine::Vina) ? "vina" :
                                 (e == Engine::Smina) ? "smina" : "obabel";
        for (fs::directory_iterator it(dir, ec), end; it != end && !ec; it.increment(ec)) {
            if (!it->is_regular_file(ec)) continue;
            std::string fn = it->path().filename().string();
            std::string low = fn;
            for (auto& c : low) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            if (low.rfind(stem, 0) == 0 &&
                (low.size() >= 4 && low.substr(low.size() - 4) == ".exe"))
                return it->path();
        }
    }

    // 2) PATH.
    for (const auto& d : pathDirs()) {
        for (const auto& n : names) {
            const fs::path cand = d / n;
            if (fs::exists(cand, ec) && fs::is_regular_file(cand, ec)) return cand;
        }
    }
    return std::nullopt;
}

bool engineAvailable(Engine e) { return locateEngine(e).has_value(); }

ProvisionResult ensureVina(bool allowDownload) {
    ProvisionResult r;

    // Already provisioned? Done.
    if (auto p = locateEngine(Engine::Vina)) {
        r.fetched = true;
        r.path = p->string();
        r.note = "vina already present at " + r.path;
        return r;
    }
    if (!allowDownload) {
        r.fetched = false;
        r.note = "vina not found; download not requested (locate-only).";
        return r;
    }

    std::error_code ec;
    const fs::path dir = enginesDir();
    fs::create_directories(dir, ec);
    if (ec) {
        r.note = "could not create engines dir: " + ec.message();
        return r;
    }
    const fs::path target = dir / "vina.exe";

#if defined(_WIN32)
    // Best-effort PowerShell download with a size sanity-check + optional SHA-256
    // verify. Quoted so spaces in the profile path are safe; -ErrorAction Stop turns
    // transfer failures into a non-zero exit we can detect. We NEVER block app
    // startup on this (callers run it off the UI thread / opportunistically), and any
    // failure leaves us in the locate-only state with a note. Exit codes: 0 = ok,
    // 2 = transfer error, 3 = hash mismatch, 4 = implausibly small (error page).
    const std::string expected = expectedVinaSha256();
    std::string ps;
    ps += "$ErrorActionPreference='Stop';";
    ps += "try{";
    ps += "[Net.ServicePointManager]::SecurityProtocol=[Net.SecurityProtocolType]::Tls12;";
    ps += "Invoke-WebRequest -UseBasicParsing -Uri '";
    ps += kVinaUrl;
    ps += "' -OutFile '";
    ps += target.string();
    ps += "';";
    // Reject an HTML error page / truncated transfer masquerading as the binary.
    ps += "if((Get-Item '";
    ps += target.string();
    ps += "').Length -lt ";
    ps += std::to_string(kVinaMinSizeBytes);
    ps += "){Remove-Item -Force '";
    ps += target.string();
    ps += "';exit 4};";
    if (expected.size() == 64) {
        ps += "$h=(Get-FileHash -Algorithm SHA256 '";
        ps += target.string();
        ps += "').Hash;";
        ps += "if($h -ne '";
        ps += expected;
        ps += "'){Remove-Item -Force '";
        ps += target.string();
        ps += "';exit 3};";
    }
    ps += "exit 0";
    ps += "}catch{exit 2}";

    std::string cmd = "powershell -NoProfile -NonInteractive -Command \"" + ps + "\"";
    const int code = std::system(cmd.c_str());

    if (code == 0 && fs::exists(target, ec) && fs::file_size(target, ec) >= kVinaMinSizeBytes) {
        r.fetched = true;
        r.path = target.string();
        const long long sz = static_cast<long long>(fs::file_size(target, ec));
        const bool sizeMatch = (sz == kVinaSizeBytes);
        if (expected.size() == 64)
            r.note = "downloaded + SHA-256 verified vina.exe";
        else if (sizeMatch)
            r.note = "downloaded vina.exe (size matches published " +
                     std::to_string(kVinaSizeBytes) + " B; no pinned SHA-256 - "
                     "set BIOCAD_VINA_SHA256 or runtime/engines/vina.sha256 to verify)";
        else
            r.note = "downloaded vina.exe (" + std::to_string(sz) +
                     " B; UNVERIFIED - no pinned SHA-256 and size differs from published)";
        return r;
    }
    r.fetched = false;
    r.note = code == 3   ? "vina download SHA-256 mismatch; discarded (check your pinned hash)."
             : code == 4 ? "vina download too small (likely an error page); discarded."
                         : "vina download failed (exit " + std::to_string(code) +
                               "); staying in descriptor-estimate fallback.";
    // Remove a partial/zero-byte file so a later locate() does not trust it.
    if (fs::exists(target, ec) && fs::file_size(target, ec) < kVinaMinSizeBytes)
        fs::remove(target, ec);
    return r;
#else
    (void)target;
    r.note = "provisioning only implemented on Windows.";
    return r;
#endif
}

ProvisionResult ensureVinaGpu(bool allowDownload) {
    ProvisionResult r;
    std::error_code ec;
    const fs::path dir = vinaGpuDir();
    const fs::path exe = dir / "Vina-GPU.exe";
    const fs::path kernel = dir / "Kernel2_Opt.bin";

    const bool haveExe = fs::exists(exe, ec) && fs::file_size(exe, ec) >= kVinaGpuMinExeBytes;
    const bool haveKernel = fs::exists(kernel, ec);

    // Already fully provisioned (exe + a kernel compiled for this GPU)? Done.
    if (haveExe && haveKernel) {
        r.fetched = true;
        r.path = exe.string();
        r.note = "Vina-GPU already present (exe + compiled kernel) at " + dir.string();
        return r;
    }
    if (!allowDownload) {
        r.fetched = false;
        r.note = haveExe ? "Vina-GPU exe present but its Kernel2_Opt.bin is not compiled yet; "
                           "provision (with download) to build the kernel for this GPU."
                         : "Vina-GPU (OpenCL) not provisioned; download not requested (locate-only).";
        return r;
    }

#if defined(_WIN32)
    fs::create_directories(dir, ec);
    if (ec) {
        r.note = "could not create vina-gpu dir: " + ec.message();
        return r;
    }

    // Fetch only the needed files (NOT the 37 MB GUI), then compile the kernel for THIS
    // GPU by running Vina-GPU-K.exe once on the bundled example complex. OpenCL binaries
    // are device-specific, so we never trust a prebuilt kernel - we build our own.
    // Exit codes: 0 = exe + kernel present, 5 = exe present but no kernel produced
    // (no usable OpenCL device), 4 = exe missing after download, 2 = download error.
    std::string fileList;
    for (const auto& f : vinaGpuFiles()) {
        if (!fileList.empty()) fileList += ",";
        fileList += "'" + f + "'";
    }
    std::string ps;
    ps += "$ErrorActionPreference='Stop';";
    ps += "try{";
    ps += "[Net.ServicePointManager]::SecurityProtocol=[Net.SecurityProtocolType]::Tls12;";
    ps += "$base='"; ps += kVinaGpuRawBase; ps += "';";
    ps += "$dir='"; ps += dir.string(); ps += "';";
    ps += "$files=@("; ps += fileList; ps += ");";
    ps += "foreach($f in $files){";
    ps += "$dest=Join-Path $dir ($f -replace '/','\\');";
    ps += "$pdir=Split-Path -Parent $dest;";
    ps += "if(-not (Test-Path $pdir)){New-Item -ItemType Directory -Force -Path $pdir|Out-Null};";
    ps += "Invoke-WebRequest -UseBasicParsing -Uri ($base+$f) -OutFile $dest;";
    ps += "}";
    ps += "$exe=Join-Path $dir 'Vina-GPU.exe';";
    ps += "if(-not (Test-Path $exe) -or (Get-Item $exe).Length -lt ";
    ps += std::to_string(kVinaGpuMinExeBytes);
    ps += "){exit 4};";
    ps += "$k=Join-Path $dir 'Vina-GPU-K.exe';";
    ps += "$a='--receptor input_file_example/2bm2_protein.pdbqt --ligand input_file_example/2bm2_ligand.pdbqt"
          " --out _kernelgen_out.pdbqt --center_x 40.415 --center_y 110.986 --center_z 82.673"
          " --size_x 25 --size_y 25 --size_z 25 --thread 1000';";
    ps += "$p=Start-Process -FilePath $k -ArgumentList $a -WorkingDirectory $dir -PassThru -NoNewWindow"
          " -RedirectStandardOutput (Join-Path $dir '_kgen.log') -RedirectStandardError (Join-Path $dir '_kgen.err');";
    ps += "$null=$p.WaitForExit(300000);";
    ps += "if(-not $p.HasExited){try{$p.Kill()}catch{}};";
    ps += "if(Test-Path (Join-Path $dir 'Kernel2_Opt.bin')){exit 0}else{exit 5};";
    ps += "}catch{exit 2}";

    std::string cmd = "powershell -NoProfile -NonInteractive -Command \"" + ps + "\"";
    const int code = std::system(cmd.c_str());

    const bool nowExe = fs::exists(exe, ec) && fs::file_size(exe, ec) >= kVinaGpuMinExeBytes;
    const bool nowKernel = fs::exists(kernel, ec);
    if (code == 0 && nowExe && nowKernel) {
        r.fetched = true;
        r.path = exe.string();
        r.note = "downloaded Vina-GPU + compiled Kernel2_Opt.bin for this GPU";
        return r;
    }
    r.fetched = false;
    r.path = nowExe ? exe.string() : "";
    r.note = code == 5 ? "Vina-GPU downloaded but the kernel compile produced no Kernel2_Opt.bin "
                         "(no usable OpenCL device?); GPU-OpenCL docking stays unavailable."
           : code == 4 ? "Vina-GPU download incomplete (exe missing or too small); discarded."
                       : "Vina-GPU provisioning failed (exit " + std::to_string(code) +
                             "); GPU-OpenCL docking stays unavailable.";
    return r;
#else
    r.note = "Vina-GPU provisioning only implemented on Windows.";
    return r;
#endif
}

ProvisionResult ensureObabel() {
    ProvisionResult r;
    if (auto p = locateEngine(Engine::Obabel)) {
        r.fetched = true;
        r.path = p->string();
        r.note = "obabel located at " + r.path + " (used for receptor protonation).";
    } else {
        r.fetched = false;
        r.note = "obabel.exe not found; receptor prep uses the built-in heavy-atom writer. "
                 "Drop obabel.exe into runtime/engines for full protonation.";
    }
    return r;
}

}  // namespace biocad::docking
