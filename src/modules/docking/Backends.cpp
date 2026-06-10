#include "modules/docking/Backends.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

#include "chem/Embed3D.h"

#include "chem/Descriptors.h"
#include "core/AppPaths.h"
#include "modules/docking/PdbqtWriter.h"
#include "modules/docking/VinaParser.h"

#ifdef _WIN32
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#endif

namespace stimlab::docking {
namespace {

namespace fs = std::filesystem;
namespace chem = stimlab::chem;

// Run a process to completion with no console window; return its exit code (or -1
// if it could not be launched). Windows-only; the fallback path never calls this.
int runProcess(const std::wstring& cmdline, const fs::path& workdir, unsigned timeoutMs = 180000) {
#ifdef _WIN32
    STARTUPINFOW si{};
    si.cb = sizeof(si);
    PROCESS_INFORMATION pi{};
    std::wstring buf = cmdline;  // CreateProcessW may write to the command buffer
    const std::wstring wd = workdir.wstring();
    BOOL ok = CreateProcessW(nullptr, buf.data(), nullptr, nullptr, FALSE, CREATE_NO_WINDOW,
                             nullptr, wd.empty() ? nullptr : wd.c_str(), &si, &pi);
    if (!ok) return -1;
    if (WaitForSingleObject(pi.hProcess, timeoutMs) == WAIT_TIMEOUT) {
        TerminateProcess(pi.hProcess, 1);            // a hung search must not wedge the worker
        WaitForSingleObject(pi.hProcess, 5000);
        CloseHandle(pi.hProcess);
        CloseHandle(pi.hThread);
        return -2;  // distinct from a clean non-zero exit
    }
    DWORD code = 1;
    GetExitCodeProcess(pi.hProcess, &code);
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);
    return static_cast<int>(code);
#else
    (void)cmdline;
    (void)workdir;
    (void)timeoutMs;
    return -1;
#endif
}

std::string readFile(const fs::path& p) {
    std::ifstream in(p, std::ios::binary);
    if (!in) return {};
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

// Convergent FLEXIBLE-ligand search shared by the CPU and GPU Vina backends. It
// repeats an INDEPENDENT search (fresh --seed, escalating --exhaustiveness on CPU)
// with a torsionally flexible ligand and keeps going until the best affinity
// stabilises within a tolerance across two consecutive runs, or a run/time budget
// is reached - so a single noisy search can't produce a low-confidence number.
// Falls back to a rigid ligand if the flexible torsion tree is rejected on the
// first run, so a flex-writer bug can never lose a dock that would otherwise work.
DockJobResult runVinaSearch(const std::string& engineName, const fs::path& bin,
                            const fs::path& workdir, bool gpu, const std::string& tag,
                            const chem::Molecule& graph, const chem::Conformer& ligand3d,
                            const ReceptorTarget& target) {
    DockJobResult r;
    r.engine = engineName;
    r.real = false;
    r.targetId = target.id;

    const PdbqtLigand flex = writeFlexiblePdbqt(graph, ligand3d);
    const PdbqtLigand rigid = writeRigidPdbqt(graph, ligand3d);
    if (flex.atomCount == 0 && rigid.atomCount == 0) {
        r.log = "Ligand PDBQT preparation produced no atoms.";
        return r;
    }

    const fs::path dir = AppPaths::instance().cache();
    std::error_code ec;
    fs::create_directories(dir, ec);
    const fs::path ligPath = dir / ("lig_" + tag + ".pdbqt");
    const fs::path outPath = dir / ("dock_" + tag + ".pdbqt");

    bool useFlexible = flex.atomCount > 0;
    int  torsions = useFlexible ? flex.torsions : 0;
    const std::vector<int>* activeOrder = useFlexible ? &flex.serialToConf : &rigid.serialToConf;
    auto writeLig = [&](const std::string& text) {
        std::ofstream o(ligPath, std::ios::binary);
        o << text;
    };
    writeLig(useFlexible ? flex.text : rigid.text);

    auto num = [](double v) { return std::to_wstring(v); };
    auto clampSize = [](double s) { return s > 28.0 ? 28.0 : (s < 1.0 ? 1.0 : s); };  // Vina-GPU < 30 A
    auto buildCmd = [&](int seed, int exh) -> std::wstring {
        std::wstring c = L"\"" + bin.wstring() + L"\"" +
            L" --receptor \"" + fs::path(target.receptorPath).wstring() + L"\"" +
            L" --ligand \"" + ligPath.wstring() + L"\"" +
            L" --out \"" + outPath.wstring() + L"\"" +
            L" --seed " + std::to_wstring(seed) +
            L" --center_x " + num(target.box.cx) + L" --center_y " + num(target.box.cy) +
            L" --center_z " + num(target.box.cz);
        if (gpu) {
            // Vina-GPU has no --exhaustiveness; thoroughness comes from --thread lanes +
            // independent seeds. Box is clamped < 30 A per the engine's limitation.
            c += L" --thread 8000";
            c += L" --size_x " + num(clampSize(target.box.sx)) +
                 L" --size_y " + num(clampSize(target.box.sy)) +
                 L" --size_z " + num(clampSize(target.box.sz));
        } else {
            c += L" --exhaustiveness " + std::to_wstring(exh);
            c += L" --size_x " + num(target.box.sx) + L" --size_y " + num(target.box.sy) +
                 L" --size_z " + num(target.box.sz);
        }
        return c;
    };
    auto runOnce = [&](int seed, int exh) -> std::vector<DockPose> {
        fs::remove(outPath, ec);  // never let a stale output misparse as this run's poses
        const int code = runProcess(buildCmd(seed, exh), workdir, 300000);
        if (code != 0) return {};
        return parseVinaPdbqt(readFile(outPath), ligand3d, activeOrder);
    };

    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    constexpr int  kMaxRuns = 5;
    constexpr long long kBudgetMs = 150000;  // stop STARTING new runs past ~2.5 min
    constexpr double kTol = 0.2;             // kcal/mol: best affinity considered "stable"
    const int exhSchedule[5] = {16, 24, 32, 32, 32};

    std::vector<DockPose> best;
    double bestAff = 1e9;
    std::vector<double> runBests;
    bool converged = false;

    for (int i = 0; i < kMaxRuns; ++i) {
        const long long elapsed =
            std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - t0).count();
        if (i >= 2 && elapsed > kBudgetMs) break;  // always do >= 2 runs, then respect the budget

        auto poses = runOnce(i + 1, exhSchedule[i]);
        if (poses.empty() && useFlexible && i == 0) {
            useFlexible = false; torsions = 0;     // flexible tree rejected -> rigid fallback
            activeOrder = &rigid.serialToConf;
            writeLig(rigid.text);
            poses = runOnce(i + 1, exhSchedule[i]);
        }
        if (poses.empty()) continue;               // failed run; keep trying within budget

        const double a = poses.front().affinityKcalPerMol;
        runBests.push_back(a);
        if (a < bestAff) { bestAff = a; best = poses; }
        if (runBests.size() >= 2 &&
            std::abs(runBests[runBests.size() - 1] - runBests[runBests.size() - 2]) <= kTol) {
            converged = true;
            break;
        }
    }

    if (best.empty()) {
        r.log = "Ran " + engineName + " but parsed no poses from its output.";
        return r;
    }
    double mn = runBests.front(), mx = runBests.front();
    for (double v : runBests) { mn = std::min(mn, v); mx = std::max(mx, v); }
    r.poses = std::move(best);
    r.real = true;
    r.searchRuns = static_cast<int>(runBests.size());
    r.affinitySpread = mx - mn;
    r.converged = converged;
    r.torsions = torsions;
    char note[256];
    std::snprintf(note, sizeof note,
                  "%s: %d run(s), best %.2f kcal/mol, spread %.2f; %s; %d torsions (%s).",
                  engineName.c_str(), r.searchRuns, bestAff, r.affinitySpread,
                  converged ? "converged" : "budget/limit reached", torsions,
                  useFlexible ? "flexible" : "rigid fallback");
    r.log = note;
    return r;
}

}  // namespace

// ----------------------------------------------------------- EstimateBackend
DockJobResult EstimateBackend::dock(const chem::Molecule& graph, const chem::Conformer& ligand3d,
                                    const ReceptorTarget& target) const {
    DockJobResult r;
    r.engine = "descriptor-estimate";
    r.real = false;
    r.targetId = target.id;

    // Same descriptor relationship the Phase-C estimate used (logP + size), now
    // computed from the parsed graph rather than a stored field.
    const double logP = chem::crippenLogP(graph);
    const double mw = chem::molecularWeight(graph);
    const double base = -(5.0 + logP * 0.8 + mw * 0.004);
    for (int i = 0; i < 6; ++i) {
        DockPose p;
        p.rank = i + 1;
        p.affinityKcalPerMol = base + i * 0.35;
        p.rmsdLb = i * 0.9;
        p.rmsdUb = i * 0.9 + 1.2;
        p.ligand = ligand3d;  // the embedded ligand, so the viewer has geometry
        r.poses.push_back(std::move(p));
    }
    r.log = "No docking engine provisioned; affinity is a structure-descriptor estimate "
            "(logP/MW), NOT a docked score. Provision vina.exe + a prepared receptor under "
            "runtime/engines to enable real docking.";
    return r;
}

// --------------------------------------------------------------- VinaBackend
std::string VinaBackend::id() const { return engine_ == Engine::Smina ? "smina" : "vina"; }
std::string VinaBackend::displayName() const {
    return engine_ == Engine::Smina ? "smina" : "AutoDock Vina";
}
bool VinaBackend::available() const { return engineAvailable(engine_); }

DockJobResult VinaBackend::dock(const chem::Molecule& graph, const chem::Conformer& ligand3d,
                                const ReceptorTarget& target) const {
    DockJobResult r;
    r.engine = displayName();
    r.real = false;
    r.targetId = target.id;

    const auto bin = locateEngine(engine_);
    if (!bin) {
        r.log = displayName() + " binary not found under runtime/engines or PATH.";
        return r;
    }
    if (target.receptorPath.empty() || !fs::exists(target.receptorPath)) {
        r.log = "No prepared receptor for " + target.name +
                " (drop a PDBQT under runtime/engines or presets/receptors); cannot dock.";
        return r;
    }

    // Flexible-ligand, convergent search (bends/rotates the ligand; repeats until the
    // best affinity stabilises). Falls back to a rigid ligand internally if needed.
    return runVinaSearch(displayName(), *bin, AppPaths::instance().cache(), /*gpu=*/false,
                         target.id, graph, ligand3d, target);
}

// ------------------------------------------------------------ VinaGpuBackend
bool VinaGpuBackend::available() const {
    const auto bin = locateEngine(Engine::VinaGpu);
    if (!bin) return false;
    // Vina-GPU.exe loads a precompiled OpenCL kernel; without it the engine can't run.
    std::error_code ec;
    return fs::exists(bin->parent_path() / "Kernel2_Opt.bin", ec);
}

DockJobResult VinaGpuBackend::dock(const chem::Molecule& graph, const chem::Conformer& ligand3d,
                                   const ReceptorTarget& target) const {
    DockJobResult r;
    r.engine = displayName();
    r.real = false;
    r.targetId = target.id;

    const auto bin = locateEngine(Engine::VinaGpu);
    if (!bin) {
        r.log = "Vina-GPU.exe not found under runtime/engines/vina-gpu or PATH.";
        return r;
    }
    const fs::path engineDir = bin->parent_path();
    std::error_code ec;
    if (!fs::exists(engineDir / "Kernel2_Opt.bin", ec)) {
        r.log = "Vina-GPU present but its Kernel2_Opt.bin is not compiled; provision it to "
                "build the OpenCL kernel for this GPU.";
        return r;
    }
    if (target.receptorPath.empty() || !fs::exists(target.receptorPath)) {
        r.log = "No prepared receptor for " + target.name + "; cannot GPU-dock.";
        return r;
    }

    // Flexible-ligand, convergent search on the GPU (runs in the engine's own folder so
    // it finds Kernel2_Opt.bin). Same torsion-tree ligand + multi-seed convergence as
    // the CPU path; --thread lanes replace --exhaustiveness for the GPU engine.
    return runVinaSearch(displayName(), *bin, engineDir, /*gpu=*/true, "gpu_" + target.id,
                         graph, ligand3d, target);
}

}  // namespace stimlab::docking
