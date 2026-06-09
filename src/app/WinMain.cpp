// app/WinMain.cpp - Win32 + DirectX 11 + Dear ImGui entry point.
// Creates the window, brings up the device and ImGui, wires the FakeBackend into
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
#include "modules/docking/Presets.h"
#include "modules/docking/Provisioning.h"
#include "render/Dx11Device.h"
#include "ui/AppShell.h"
#include "ui/Theme.h"

// Forward decl from imgui_impl_win32 (handles input/message translation).
extern IMGUI_IMPL_API LRESULT ImGui_ImplWin32_WndProcHandler(HWND, UINT, WPARAM, LPARAM);

namespace {
stimlab::Dx11Device* g_device = nullptr;

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
int runHeadlessDock(const std::string& smiles, const std::string& target) {
    using namespace stimlab;
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

    emit("== StimLab docking selftest (real-engine acceptance) ==");
    emit("ligand : " + smiles);
    emit("target : " + target);

    docking::Provisioner prov;
    prov.start(/*allowDownload=*/true, docking::headlinePresets());
    while (prov.running()) Sleep(250);
    emit("provision : " + prov.status());
    emit("vina      : " + std::string(prov.vinaReady() ? "ready" : "absent") +
         "   receptors: " + std::to_string(prov.receptorsReady()) + "/" +
         std::to_string(prov.receptorsTotal()));

    RealBackend backend;
    Services s = backend.services();
    Molecule lig;
    lig.id = "__cli__";
    lig.name = "cli-ligand";
    lig.smiles = smiles;

    const DockJobResult d = s.docking->dockDetailed(lig, target);
    emit("engine    : " + d.engine);
    emit(std::string("REAL dock : ") + (d.real ? "YES (engine dock)" : "no (descriptor estimate)"));
    emit("poses     : " + std::to_string(d.poses.size()));
    if (!d.poses.empty()) {
        char b[64];
        std::snprintf(b, sizeof(b), "%.2f", d.poses.front().affinityKcalPerMol);
        emit(std::string("best aff  : ") + b + " kcal/mol");
    }
    emit("log       : " + d.log);
    emit(d.real ? "RESULT    : PASS - real ranked poses produced."
                : "RESULT    : fell back to labeled estimate (see log).");

    {
        std::ofstream o(AppPaths::instance().runtime() / "selftest-dock.txt", std::ios::binary);
        o << rep.str();
    }
    if (fp) std::fclose(fp);
    if (con) FreeConsole();
    log::shutdown();
    return d.real ? 0 : 2;
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

int WINAPI wWinMain(HINSTANCE hInstance, HINSTANCE, PWSTR, int nCmdShow) {
    // Headless modes (no GUI): "--selftest-dock [--smiles S] [--target T]" provisions
    // the engine + receptors and docks one ligand through the real backend path.
    {
        int argc = 0;
        LPWSTR* argv = CommandLineToArgvW(GetCommandLineW(), &argc);
        bool selftest = false;
        std::string smiles = "CC(N)Cc1ccccc1";  // amphetamine
        std::string target = "DAT";
        for (int i = 1; i < argc; ++i) {
            const std::wstring a = argv[i];
            if (a == L"--selftest-dock") selftest = true;
            else if (a == L"--smiles" && i + 1 < argc) smiles = wToUtf8(argv[++i]);
            else if (a == L"--target" && i + 1 < argc) target = wToUtf8(argv[++i]);
        }
        if (argv) LocalFree(argv);
        if (selftest) return runHeadlessDock(smiles, target);
    }

    stimlab::AppPaths::instance().ensureLayout();
    stimlab::log::init();
    spdlog::info("StimLab starting (Phase D: real chem engine + 3D viewer + docking)");

    ImGui_ImplWin32_EnableDpiAwareness();

    WNDCLASSEXW wc{};
    wc.cbSize = sizeof(wc);
    wc.style = CS_HREDRAW | CS_VREDRAW;
    wc.lpfnWndProc = WndProc;
    wc.hInstance = hInstance;
    wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wc.lpszClassName = L"StimLabWindow";
    RegisterClassExW(&wc);

    HWND hwnd = CreateWindowExW(0, wc.lpszClassName, L"StimLab - CNS-stimulant analysis suite",
                                WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT, 1480, 920,
                                nullptr, nullptr, hInstance, nullptr);
    if (!hwnd) {
        spdlog::error("CreateWindow failed");
        return 1;
    }

    stimlab::Dx11Device device;
    if (!device.init(hwnd)) {
        spdlog::error("DirectX 11 device init failed");
        MessageBoxW(hwnd, L"Failed to initialize DirectX 11.", L"StimLab", MB_OK | MB_ICONERROR);
        return 1;
    }
    g_device = &device;

    ShowWindow(hwnd, nCmdShow);
    UpdateWindow(hwnd);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;

    static std::string iniPath = (stimlab::AppPaths::instance().root() / "imgui.ini").string();
    io.IniFilename = iniPath.c_str();

    // Crisp UI font: prefer Segoe UI (ships on Windows); fall back to ImGui default.
    {
        const char* segoe = "C:\\Windows\\Fonts\\segoeui.ttf";
        if (std::filesystem::exists(segoe)) {
            ImFontConfig cfg;
            cfg.OversampleH = 2;
            cfg.OversampleV = 2;
            io.Fonts->AddFontFromFileTTF(segoe, 18.0f, &cfg);
        }
    }

    stimlab::theme::apply();
    ImGui_ImplWin32_Init(hwnd);
    ImGui_ImplDX11_Init(device.device(), device.context());

    // Persisted settings (incl. the DPAPI-encrypted agent API key). Declared
    // before the shell so it outlives the shell's agent during teardown.
    stimlab::Config config(stimlab::AppPaths::instance().config());
    config.load();

    stimlab::RealBackend backend;
    stimlab::AppShell shell(backend.services());
    // WP-2: give the shell the live DX11 device so the 3D molecular viewport can
    // render into an off-screen target shown by ImGui in the Structure/Docking panels.
    shell.setRenderDevice(device.device(), device.context());
    // WP-4 (agent): hand the shell the config so the assistant loop can read the
    // provider/key/model and the Settings panel can persist them.
    shell.setConfig(&config);

    // Optional: STIMLAB_PANEL / STIMLAB_MOLECULE select the initial panel + compound
    // (handy for screenshots/automation; defaults are Dashboard/amphetamine).
    if (char buf[128]; GetEnvironmentVariableA("STIMLAB_PANEL", buf, sizeof(buf)) > 0)
        shell.state().activePanel = buf;
    if (char buf[128]; GetEnvironmentVariableA("STIMLAB_MOLECULE", buf, sizeof(buf)) > 0)
        shell.state().selectedMolecule = buf;

    spdlog::info("UI ready - {} compounds in library (real chem engine)",
                 backend.services().library->count());

    const float clear[4] = {0.06f, 0.08f, 0.11f, 1.0f};
    bool running = true;
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
        device.present(true);

        if (shell.state().quitRequested) running = false;
    }

    spdlog::info("StimLab shutting down");
    ImGui_ImplDX11_Shutdown();
    ImGui_ImplWin32_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    g_device = nullptr;
    device.shutdown();
    DestroyWindow(hwnd);
    UnregisterClassW(wc.lpszClassName, hInstance);
    stimlab::log::shutdown();
    return 0;
}
