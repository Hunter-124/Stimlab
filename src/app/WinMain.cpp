// app/WinMain.cpp - Win32 + DirectX 11 + Dear ImGui entry point.
// Creates the window, brings up the device and ImGui, wires the FakeBackend into
// the AppShell, and runs the frame loop until the user exits.
#ifndef WIN32_LEAN_AND_MEAN
#  define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>

#include <filesystem>
#include <string>

#include <imgui.h>
#include <implot.h>
#include <imgui_impl_dx11.h>
#include <imgui_impl_win32.h>

#include "core/AppPaths.h"
#include "core/Log.h"
#include "modules/RealBackend.h"
#include "render/Dx11Device.h"
#include "ui/AppShell.h"
#include "ui/Theme.h"

// Forward decl from imgui_impl_win32 (handles input/message translation).
extern IMGUI_IMPL_API LRESULT ImGui_ImplWin32_WndProcHandler(HWND, UINT, WPARAM, LPARAM);

namespace {
stimlab::Dx11Device* g_device = nullptr;
}

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
    stimlab::AppPaths::instance().ensureLayout();
    stimlab::log::init();
    spdlog::info("StimLab starting (Phase B skeleton)");

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

    stimlab::RealBackend backend;
    stimlab::AppShell shell(backend.services());
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
