// app/WinMain.cpp - Win32 + DirectX 11 + Dear ImGui entry point.
// Creates the window, brings up the device and ImGui, wires the RealBackend into
// the AppShell, and runs the frame loop until the user exits.
#ifndef WIN32_LEAN_AND_MEAN
#  define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <shellapi.h>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include <imgui.h>
#include <implot.h>
#include <imgui_impl_dx11.h>
#include <imgui_impl_win32.h>

#include "contracts/Services.h"
#include "core/AppPaths.h"
#include "core/Config.h"
#include "core/Log.h"
#include "data/Domain.h"
#include "modules/RealBackend.h"
#include "modules/docking/EngineLocator.h"
#include "modules/docking/Presets.h"
#include "modules/docking/Provisioning.h"
#include "render/Dx11Device.h"
#include "ui/AppShell.h"
#include "ui/Theme.h"

// App icon resource id (see src/app/BioCAD.rc, compiled into the exe).
#define IDI_APPICON 101

// Forward decl from imgui_impl_win32 (handles input/message translation).
extern IMGUI_IMPL_API LRESULT ImGui_ImplWin32_WndProcHandler(HWND, UINT, WPARAM, LPARAM);

namespace {
biocad::Dx11Device* g_device = nullptr;

std::string wToUtf8(const std::wstring& w) {
    if (w.empty()) return {};
    int n = WideCharToMultiByte(CP_UTF8, 0, w.c_str(), -1, nullptr, 0, nullptr, nullptr);
    std::string s(n > 0 ? n - 1 : 0, '\0');
    if (n > 1) WideCharToMultiByte(CP_UTF8, 0, w.c_str(), -1, s.data(), n, nullptr, nullptr);
    return s;
}

// Headless docking selftest / provisioning hook (no GUI). Provisions vina + the
// headline receptors, docks one ligand into one target through the REAL backend
// path, and writes a report to stdout (parent console) AND runtime/selftest-dock.txt.
// Exit code 0 = a real engine dock succeeded; 2 = fell back to the descriptor
// estimate; other = setup error. This is the WP-3 acceptance gate made runnable.
int runHeadlessDock(const std::string& smiles, const std::string& target, const std::string& compute) {
    using namespace biocad;
    AppPaths::instance().ensureLayout();
    log::init();

    const bool con = AttachConsole(ATTACH_PARENT_PROCESS) != 0 || AllocConsole() != 0;
    FILE* fp = nullptr;
    if (con) freopen_s(&fp, "CONOUT$", "w", stdout);

    std::ostringstream rep;
    auto emit = [&](const std::string& s) {
        rep << s << "\n";
        std::printf("%s\n", s.c_str());
        std::fflush(stdout);
        spdlog::info("[selftest-dock] {}", s);
    };

    emit("== BioCAD docking selftest (real-engine acceptance) ==");
    emit("ligand : " + smiles);
    emit("target : " + target);

    // Optional compute-mode override (auto/gpu/cpu) so the GPU engine is headlessly
    // testable: --selftest-dock --compute gpu exercises the first-party CUDA backend.
    docking::ComputeMode mode = docking::ComputeMode::Auto;
    if (compute == "gpu") mode = docking::ComputeMode::Gpu;
    else if (compute == "cpu") mode = docking::ComputeMode::Cpu;
    docking::setComputeMode(mode);
    emit("compute: " + compute);

    docking::Provisioner prov;
    // Provision exactly the requested target (resolved by id or display name) so a
    // non-headline preset is genuinely prepared before we dock it; fall back to the
    // headline set only when the target string matches no preset.
    std::vector<ReceptorTarget> toProvision;
    if (const auto* p = docking::findPreset(target)) toProvision.push_back(*p);
    else toProvision = docking::headlinePresets();
    prov.start(/*allowDownload=*/true, std::move(toProvision));
    while (prov.running()) Sleep(250);
    emit("provision : " + prov.status());
    emit("vina      : " + std::string(prov.vinaReady() ? "ready" : "absent") +
         "   receptors: " + std::to_string(prov.receptorsReady()) + "/" +
         std::to_string(prov.receptorsTotal()));

    // In GPU mode, also provision the Vina-GPU (OpenCL) engine so its real-engine path is
    // exercised headlessly. It is tried before the first-party CUDA engine, so when it
    // provisions successfully a `--compute gpu` dock runs on Vina-GPU; otherwise the dock
    // falls through to CUDA (if built) or the labeled estimate.
    if (compute == "gpu") {
        emit("vina-gpu  : provisioning Vina-GPU (OpenCL)...");
        const auto vg = docking::ensureVinaGpu(/*allowDownload=*/true);
        emit("vina-gpu  : " + std::string(vg.fetched ? "ready - " : "unavailable - ") + vg.note);
    }

    RealBackend backend;
    Services s = backend.services();
    Molecule lig;
    lig.id = "__cli__";
    lig.name = "cli-ligand";
    lig.smiles = smiles;

    const DockJobResult d = s.docking->dockDetailed(lig, target);
    emit("engine    : " + d.engine);
    emit(std::string("REAL dock : ") + (d.fromEngine() ? "YES (engine dock)" : "no (descriptor estimate)"));
    emit("poses     : " + std::to_string(d.poses.size()));
    if (!d.poses.empty()) {
        char b[64];
        std::snprintf(b, sizeof(b), "%.2f", d.poses.front().affinityKcalPerMol);
        emit(std::string("best aff  : ") + b + " kcal/mol");
    }
    emit("log       : " + d.log);
    emit(d.fromEngine() ? "RESULT    : PASS - real ranked poses produced."
                : "RESULT    : fell back to labeled estimate (see log).");

    {
        std::ofstream o(AppPaths::instance().runtime() / "selftest-dock.txt", std::ios::binary);
        o << rep.str();
    }
    if (fp) std::fclose(fp);
    if (con) FreeConsole();
    log::shutdown();
    return d.fromEngine() ? 0 : 2;
}
}  // namespace

static LRESULT WINAPI WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    if (ImGui_ImplWin32_WndProcHandler(hwnd, msg, wParam, lParam)) return true;
    switch (msg) {
        case WM_SIZE:
            if (g_device && wParam != SIZE_MINIMIZED)
                g_device->resize(static_cast<unsigned>(LOWORD(lParam)),
                                 static_cast<unsigned>(HIWORD(lParam)));
            return 0;
        case WM_SYSCOMMAND:
            if ((wParam & 0xFFF0) == SC_KEYMENU) return 0;  // disable ALT app menu
            break;
        case WM_DESTROY:
            PostQuitMessage(0);
            return 0;
        default:
            break;
    }
    return DefWindowProcW(hwnd, msg, wParam, lParam);
}

namespace {

// Deterministic-capture options. --shot renders a fixed number of warm-up frames so
// async provisioning, layout and animation settle, then writes PNGs and exits. This
// replaces desktop screen-scraping, which cannot work on a runner with no
// interactive desktop session.
struct ShotOptions {
    std::filesystem::path out;      // empty => capture disabled
    std::string panel;              // panel id forced before the first frame
    int frames = 1;                 // >1 writes <stem>-0000.png, <stem>-0001.png, ...
    int warmup = 120;               // frames rendered before the first capture
    int width = 1600;
    int height = 1000;

    [[nodiscard]] bool enabled() const { return !out.empty(); }

    // Path for capture index i: the single-shot case keeps the exact requested name
    // so callers can predict it.
    [[nodiscard]] std::filesystem::path pathFor(int i) const {
        if (frames <= 1) return out;
        char suffix[16];
        std::snprintf(suffix, sizeof suffix, "-%04d", i);
        auto p = out;
        p.replace_filename(out.stem().string() + suffix + out.extension().string());
        return p;
    }
};

}  // namespace

int WINAPI wWinMain(HINSTANCE hInstance, HINSTANCE, PWSTR, int nCmdShow) {
    // Headless modes (no GUI): "--selftest-dock [--smiles S] [--target T]" provisions
    // the engine + receptors and docks one ligand through the real backend path.
    ShotOptions shot;
    {
        int argc = 0;
        LPWSTR* argv = CommandLineToArgvW(GetCommandLineW(), &argc);
        bool selftest = false;
        std::string smiles = "Cn1cnc2c1c(=O)n(C)c(=O)n2C";  // caffeine
        std::string target = "DAT";
        std::string compute = "auto";  // auto | gpu | cpu
        for (int i = 1; i < argc; ++i) {
            const std::wstring a = argv[i];
            if (a == L"--selftest-dock") selftest = true;
            else if (a == L"--smiles" && i + 1 < argc) smiles = wToUtf8(argv[++i]);
            else if (a == L"--target" && i + 1 < argc) target = wToUtf8(argv[++i]);
            else if (a == L"--compute" && i + 1 < argc) compute = wToUtf8(argv[++i]);
            else if (a == L"--shot" && i + 1 < argc) shot.out = std::filesystem::path(argv[++i]);
            else if (a == L"--shot-panel" && i + 1 < argc) shot.panel = wToUtf8(argv[++i]);
            else if (a == L"--shot-frames" && i + 1 < argc) shot.frames = std::stoi(wToUtf8(argv[++i]));
            else if (a == L"--shot-warmup" && i + 1 < argc) shot.warmup = std::stoi(wToUtf8(argv[++i]));
            else if (a == L"--shot-size" && i + 2 < argc) {
                shot.width = std::stoi(wToUtf8(argv[++i]));
                shot.height = std::stoi(wToUtf8(argv[++i]));
            }
        }
        if (argv) LocalFree(argv);
        if (selftest) return runHeadlessDock(smiles, target, compute);
        if (shot.frames < 1) shot.frames = 1;
        if (shot.warmup < 0) shot.warmup = 0;
    }

    biocad::AppPaths::instance().ensureLayout();
    biocad::log::init();
    spdlog::info("BioCAD starting (real chem engine + 3D viewer + docking + workflows + agent)");

    // Self-heal the provisioned runtime: delete any corrupt engine/receptor so it is
    // re-fetched on the next provision (manifest.json is the source of truth).
    if (const auto st = biocad::docking::selfHealManifest(); st.total > 0) {
        spdlog::info("Runtime manifest: {}/{} components verified ({} missing, {} corrupt-healed)",
                     st.present, st.total, st.missing.size(), st.corrupt.size());
    }

    ImGui_ImplWin32_EnableDpiAwareness();

    WNDCLASSEXW wc{};
    wc.cbSize = sizeof(wc);
    wc.style = CS_HREDRAW | CS_VREDRAW;
    wc.lpfnWndProc = WndProc;
    wc.hInstance = hInstance;
    wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
    // Brand icon (title bar + taskbar). The same resource is the exe's Explorer icon.
    wc.hIcon = static_cast<HICON>(LoadImageW(hInstance, MAKEINTRESOURCEW(IDI_APPICON), IMAGE_ICON,
                                             0, 0, LR_DEFAULTSIZE | LR_SHARED));
    wc.hIconSm = static_cast<HICON>(LoadImageW(hInstance, MAKEINTRESOURCEW(IDI_APPICON), IMAGE_ICON,
                                               GetSystemMetrics(SM_CXSMICON),
                                               GetSystemMetrics(SM_CYSMICON), LR_SHARED));
    wc.lpszClassName = L"BioCADWindow";
    RegisterClassExW(&wc);

    // A capture run uses a fixed CLIENT area so every image in docs/media has the
    // same dimensions regardless of DPI, theme metrics or the runner's desktop.
    RECT wantClient{0, 0, shot.enabled() ? shot.width : 1480, shot.enabled() ? shot.height : 920};
    AdjustWindowRect(&wantClient, WS_OVERLAPPEDWINDOW, FALSE);
    HWND hwnd = CreateWindowExW(0, wc.lpszClassName, L"BioCAD",
                                WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT,
                                wantClient.right - wantClient.left,
                                wantClient.bottom - wantClient.top,
                                nullptr, nullptr, hInstance, nullptr);
    if (!hwnd) {
        spdlog::error("CreateWindow failed");
        return 1;
    }

    biocad::Dx11Device device;
    if (!device.init(hwnd)) {
        spdlog::error("DirectX 11 device init failed");
        MessageBoxW(hwnd, L"Failed to initialize DirectX 11.", L"BioCAD", MB_OK | MB_ICONERROR);
        return 1;
    }
    g_device = &device;

    // WIC (the --shot PNG encoder) needs an initialised apartment. Harmless when no
    // capture is requested, and RPC_E_CHANGED_MODE just means someone got here first.
    const HRESULT comInit = CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED);

    // A capture run keeps the window off-screen-ish but realised: SW_SHOWNA avoids
    // stealing focus on a developer desktop while still giving DXGI a valid target.
    ShowWindow(hwnd, shot.enabled() ? SW_SHOWNA : nCmdShow);
    UpdateWindow(hwnd);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;

    static std::string iniPath = (biocad::AppPaths::instance().root() / "imgui.ini").string();
    io.IniFilename = iniPath.c_str();

    // Crisp UI fonts: Segoe UI (body) + Segoe UI Semibold (semi), both scalable.
    // The first font added becomes the ImGui default face.
    {
        const char* segoeBody = "C:\\Windows\\Fonts\\segoeui.ttf";
        const char* segoeSemi = "C:\\Windows\\Fonts\\seguisb.ttf";

        ImFont* bodyPtr = nullptr;
        ImFont* semiPtr = nullptr;

        if (std::filesystem::exists(segoeBody)) {
            ImFontConfig cfg;
            cfg.OversampleH = 2;
            cfg.OversampleV = 2;
            bodyPtr = io.Fonts->AddFontFromFileTTF(segoeBody, 17.0f, &cfg);
        }

        if (std::filesystem::exists(segoeSemi)) {
            ImFontConfig cfg;
            cfg.OversampleH = 2;
            cfg.OversampleV = 2;
            semiPtr = io.Fonts->AddFontFromFileTTF(segoeSemi, 17.0f, &cfg);
        }

        biocad::theme::fonts().body = bodyPtr;
        biocad::theme::fonts().semi = (semiPtr != nullptr) ? semiPtr : bodyPtr;
    }

    biocad::theme::apply();
    ImGui_ImplWin32_Init(hwnd);
    ImGui_ImplDX11_Init(device.device(), device.context());

    // Persisted settings (incl. the DPAPI-encrypted agent API key). Declared
    // before the shell so it outlives the shell's agent during teardown.
    biocad::Config config(biocad::AppPaths::instance().config());
    config.load();

    biocad::RealBackend backend;
    biocad::AppShell shell(backend.services());
    // WP-2: give the shell the live DX11 device so the 3D molecular viewport can
    // render into an off-screen target shown by ImGui in the Structure/Docking panels.
    shell.setRenderDevice(device.device(), device.context());
    // WP-4 (agent): hand the shell the config so the assistant loop can read the
    // provider/key/model and the Settings panel can persist them.
    shell.setConfig(&config);

    // --shot-panel wins over BIOCAD_PANEL: an explicit flag beats ambient state.
    // Optional: BIOCAD_PANEL / BIOCAD_MOLECULE select the initial panel + compound
    // (handy for screenshots/automation; the default panel is Dashboard).
    if (char buf[128]; GetEnvironmentVariableA("BIOCAD_PANEL", buf, sizeof(buf)) > 0)
        shell.state().activePanel = buf;
    if (!shot.panel.empty()) shell.state().activePanel = shot.panel;
    if (char buf[128]; GetEnvironmentVariableA("BIOCAD_MOLECULE", buf, sizeof(buf)) > 0)
        shell.state().selectedMolecule = buf;
    // BIOCAD_TARGET selects the initial docking target (preset id or display name) so a
    // capture can show the on-demand "Provision <target>" affordance for an unprepared one.
    if (char buf[128]; GetEnvironmentVariableA("BIOCAD_TARGET", buf, sizeof(buf)) > 0) {
        if (const auto* p = biocad::docking::findPreset(buf)) shell.state().dockTarget = p->name;
        else shell.state().dockTarget = buf;
    }

    // Automation hook: deep-linking to the Workflows panel kicks one pipeline run so a
    // capture shows the live DAG executing (harmless in normal use - this env path is
    // only exercised by the screenshot tooling).
    if (shell.state().activePanel == "Workflows") {
        const biocad::Molecule mm = shell.currentMolecule();
        shell.runWorkflow(mm.smiles, "DAT", mm.name + " -> DAT");
    }

    spdlog::info("UI ready - {} compounds in library (real chem engine)",
                 backend.services().library->count());

    const float clear[4] = {0.06f, 0.08f, 0.11f, 1.0f};
    bool running = true;
    int frameIndex = 0;      // frames rendered so far
    int captured = 0;        // PNGs written so far
    bool captureFailed = false;
    while (running) {
        MSG msg;
        while (PeekMessageW(&msg, nullptr, 0, 0, PM_REMOVE)) {
            TranslateMessage(&msg);
            DispatchMessageW(&msg);
            if (msg.message == WM_QUIT) running = false;
        }
        if (!running) break;

        if (IsIconic(hwnd)) {
            Sleep(16);
            continue;
        }

        ImGui_ImplDX11_NewFrame();
        ImGui_ImplWin32_NewFrame();
        ImGui::NewFrame();

        shell.draw();

        ImGui::Render();
        device.beginFrame(clear);
        ImGui_ImplDX11_RenderDrawData(ImGui::GetDrawData());

        // Capture BEFORE present: the swap chain discards on Present, so a
        // post-present read-back would sample undefined memory.
        if (shot.enabled() && frameIndex >= shot.warmup && captured < shot.frames) {
            if (!device.captureBackBufferPng(shot.pathFor(captured))) {
                captureFailed = true;
                running = false;
            }
            ++captured;
        }

        device.present(!shot.enabled());  // a capture run must not wait on vsync
        ++frameIndex;

        if (shot.enabled() && captured >= shot.frames) running = false;
        if (shell.state().quitRequested) running = false;
    }

    spdlog::info("BioCAD shutting down");
    ImGui_ImplDX11_Shutdown();
    ImGui_ImplWin32_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    g_device = nullptr;
    device.shutdown();
    DestroyWindow(hwnd);
    UnregisterClassW(wc.lpszClassName, hInstance);
    if (SUCCEEDED(comInit)) CoUninitialize();
    biocad::log::shutdown();
    // Exit 3 means "a capture was requested and did not produce every PNG": CI must
    // fail loudly rather than publish a docs artifact with missing or black images.
    if (shot.enabled() && (captureFailed || captured < shot.frames)) return 3;
    return 0;
}
