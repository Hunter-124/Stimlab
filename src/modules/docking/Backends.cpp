#include "modules/docking/Backends.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <system_error>

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
int runProcess(const std::wstring& cmdline, const fs::path& workdir) {
#ifdef _WIN32
    STARTUPINFOW si{};
    si.cb = sizeof(si);
    PROCESS_INFORMATION pi{};
    std::wstring buf = cmdline;  // CreateProcessW may write to the command buffer
    const std::wstring wd = workdir.wstring();
    BOOL ok = CreateProcessW(nullptr, buf.data(), nullptr, nullptr, FALSE, CREATE_NO_WINDOW,
                             nullptr, wd.empty() ? nullptr : wd.c_str(), &si, &pi);
    if (!ok) return -1;
    WaitForSingleObject(pi.hProcess, 180000);  // 3-minute cap
    DWORD code = 1;
    GetExitCodeProcess(pi.hProcess, &code);
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);
    return static_cast<int>(code);
#else
    (void)cmdline;
    (void)workdir;
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

    // Ligand prep with OUR rigid PDBQT writer (a clean ROOT avoids the Meeko/tree.h
    // gotcha that flexible PDBQTs trip in Vina).
    const PdbqtLigand lig = writeRigidPdbqt(graph, ligand3d);
    if (lig.atomCount == 0) {
        r.log = "Ligand PDBQT preparation produced no atoms.";
        return r;
    }

    const fs::path dir = AppPaths::instance().cache();
    std::error_code ec;
    fs::create_directories(dir, ec);
    const fs::path ligPath = dir / ("lig_" + target.id + ".pdbqt");
    const fs::path outPath = dir / ("dock_" + target.id + ".pdbqt");
    {
        std::ofstream o(ligPath, std::ios::binary);
        o << lig.text;
    }

    auto num = [](double v) { return std::to_wstring(v); };
    const std::wstring cmd =
        L"\"" + bin->wstring() + L"\"" + L" --receptor \"" + fs::path(target.receptorPath).wstring() +
        L"\"" + L" --ligand \"" + ligPath.wstring() + L"\"" + L" --out \"" + outPath.wstring() + L"\"" +
        L" --center_x " + num(target.box.cx) + L" --center_y " + num(target.box.cy) +
        L" --center_z " + num(target.box.cz) + L" --size_x " + num(target.box.sx) +
        L" --size_y " + num(target.box.sy) + L" --size_z " + num(target.box.sz);

    const int code = runProcess(cmd, dir);
    if (code != 0) {
        r.log = displayName() + " subprocess exited with code " + std::to_string(code) + ".";
        return r;
    }

    auto poses = parseVinaPdbqt(readFile(outPath), ligand3d);
    if (poses.empty()) {
        r.log = "Ran " + displayName() + " but parsed no poses from its output.";
        return r;
    }
    r.poses = std::move(poses);
    r.real = true;
    r.log = "Docked with " + displayName() + " (" + std::to_string(r.poses.size()) + " poses).";
    return r;
}

}  // namespace stimlab::docking
