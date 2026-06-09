// agent/WebToolsRendered.cpp - optional WebView2 headless JS-render fetch (post-v1).
//
// webFetchRenderedHtml() loads a URL in a HEADLESS WebView2 (a full, invisible Chromium),
// lets the page's JavaScript run, then returns the rendered DOM HTML so the caller can
// extract text from JS-built pages a plain curl GET would miss. It is compiled in only
// when STIMLAB_ENABLE_WEBVIEW2 (STIMLAB_HAVE_WEBVIEW2); otherwise this file provides a
// stub returning std::nullopt so web_fetch transparently keeps using the curl path.
//
// It needs the machine-level Evergreen WebView2 Runtime at RUNTIME (preinstalled on
// Win11 and most Win10). When that Runtime is absent - or anything fails or times out -
// the function returns nullopt and web_fetch falls back to curl, so behaviour never
// regresses below the existing baseline.
//
// SAFETY SCOPE: this only RENDERS public web pages to read their text - the same scope as
// the curl path. The page JS runs in an isolated, invisible, host-object-free WebView2;
// the returned HTML is treated as untrusted third-party text exactly as before. No
// synthesis/route/precursor content is produced or solicited here.
#include "agent/WebTools.h"

#ifdef STIMLAB_HAVE_WEBVIEW2

#include <chrono>
#include <filesystem>
#include <mutex>
#include <string>
#include <thread>

#include <windows.h>
#include <wrl.h>
#include <WebView2.h>

#include <nlohmann/json.hpp>

#include "core/AppPaths.h"

namespace stimlab::agent {
namespace {

namespace fs = std::filesystem;
using Microsoft::WRL::Callback;
using Microsoft::WRL::ComPtr;

std::wstring toWide(const std::string& s) {
    if (s.empty()) return {};
    const int n = MultiByteToWideChar(CP_UTF8, 0, s.c_str(), static_cast<int>(s.size()), nullptr, 0);
    std::wstring w(static_cast<size_t>(n), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, s.c_str(), static_cast<int>(s.size()), w.data(), n);
    return w;
}

std::string toUtf8(PCWSTR w) {
    if (!w) return {};
    const int n = WideCharToMultiByte(CP_UTF8, 0, w, -1, nullptr, 0, nullptr, nullptr);
    std::string s(n > 0 ? static_cast<size_t>(n - 1) : 0, '\0');
    if (n > 1) WideCharToMultiByte(CP_UTF8, 0, w, -1, s.data(), n, nullptr, nullptr);
    return s;
}

// A reused, never-shown host window for the (invisible) WebView2 controller. WebView2
// needs a real parent HWND; we keep it zero-size and never call ShowWindow.
const wchar_t* kHostClass = L"StimLabWebView2Host";
void ensureHostClass() {
    static std::once_flag once;
    std::call_once(once, [] {
        WNDCLASSEXW wc{};
        wc.cbSize = sizeof(wc);
        wc.lpfnWndProc = DefWindowProcW;
        wc.hInstance = GetModuleHandleW(nullptr);
        wc.lpszClassName = kHostClass;
        RegisterClassExW(&wc);  // ERROR_CLASS_ALREADY_EXISTS on re-entry is fine
    });
}

// WebView2 needs a writable user-data folder (it spawns msedgewebview2.exe with a
// profile). Keep it under the same web cache root the curl tools use.
std::wstring userDataDir() {
    const fs::path d = AppPaths::instance().cache() / "web" / "wv2";
    std::error_code ec;
    fs::create_directories(d, ec);
    return d.wstring();
}

bool runtimeInstalled() {
    LPWSTR version = nullptr;
    const HRESULT hr = GetAvailableCoreWebView2BrowserVersionString(nullptr, &version);
    const bool ok = SUCCEEDED(hr) && version && version[0] != L'\0';
    if (version) CoTaskMemFree(version);
    return ok;
}

}  // namespace

bool webViewRenderAvailable() {
    // Compiled in (STIMLAB_HAVE_WEBVIEW2) AND the Evergreen Runtime is present.
    return runtimeInstalled();
}

std::optional<std::string> webFetchRenderedHtml(const std::string& url, int timeoutMs) {
    std::optional<std::string> result;

    // All WebView2 work runs on a dedicated STA thread that owns its own COM apartment
    // and message pump; the call is synchronous from the caller's perspective (join()).
    std::thread worker([&] {
        if (FAILED(CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED))) return;
        struct CoGuard {
            ~CoGuard() { CoUninitialize(); }
        } coGuard;  // declared first -> CoUninitialize runs LAST, after the ComPtrs below

        if (!runtimeInstalled()) return;  // -> nullopt -> caller uses curl

        ensureHostClass();
        HWND host = CreateWindowExW(0, kHostClass, L"", WS_OVERLAPPED, 0, 0, 0, 0, nullptr,
                                    nullptr, GetModuleHandleW(nullptr), nullptr);
        if (!host) return;

        ComPtr<ICoreWebView2Controller> controller;
        ComPtr<ICoreWebView2>           webview;
        bool          done = false;
        const std::wstring wurl = toWide(url);
        const std::wstring wudf = userDataDir();

        const HRESULT hrEnv = CreateCoreWebView2EnvironmentWithOptions(
            nullptr, wudf.c_str(), nullptr,
            Callback<ICoreWebView2CreateCoreWebView2EnvironmentCompletedHandler>(
                [&](HRESULT ec, ICoreWebView2Environment* env) -> HRESULT {
                    if (FAILED(ec) || !env) { done = true; PostQuitMessage(0); return S_OK; }
                    env->CreateCoreWebView2Controller(
                        host,
                        Callback<ICoreWebView2CreateCoreWebView2ControllerCompletedHandler>(
                            [&](HRESULT ec2, ICoreWebView2Controller* ctrl) -> HRESULT {
                                if (FAILED(ec2) || !ctrl) { done = true; PostQuitMessage(0); return S_OK; }
                                controller = ctrl;
                                controller->put_IsVisible(FALSE);
                                controller->get_CoreWebView2(&webview);
                                if (!webview) { done = true; PostQuitMessage(0); return S_OK; }
                                EventRegistrationToken navToken{};
                                webview->add_NavigationCompleted(
                                    Callback<ICoreWebView2NavigationCompletedEventHandler>(
                                        [&](ICoreWebView2*, ICoreWebView2NavigationCompletedEventArgs*) -> HRESULT {
                                            webview->ExecuteScript(
                                                L"document.documentElement.outerHTML",
                                                Callback<ICoreWebView2ExecuteScriptCompletedHandler>(
                                                    [&](HRESULT ec3, PCWSTR json) -> HRESULT {
                                                        if (SUCCEEDED(ec3) && json) {
                                                            try {
                                                                result = nlohmann::json::parse(toUtf8(json))
                                                                             .get<std::string>();
                                                            } catch (...) {
                                                            }
                                                        }
                                                        done = true;
                                                        PostQuitMessage(0);
                                                        return S_OK;
                                                    })
                                                    .Get());
                                            return S_OK;
                                        })
                                        .Get(),
                                    &navToken);
                                webview->Navigate(wurl.c_str());
                                return S_OK;
                            })
                            .Get());
                    return S_OK;
                })
                .Get());

        if (SUCCEEDED(hrEnv)) {
            // Pump messages (drives all the async callbacks) until done or the deadline.
            const auto deadline =
                std::chrono::steady_clock::now() + std::chrono::milliseconds(timeoutMs);
            MSG msg;
            while (!done) {
                const auto now = std::chrono::steady_clock::now();
                if (now >= deadline) break;
                const DWORD waitMs = static_cast<DWORD>(
                    std::chrono::duration_cast<std::chrono::milliseconds>(deadline - now).count());
                const DWORD wr =
                    MsgWaitForMultipleObjectsEx(0, nullptr, waitMs, QS_ALLINPUT, MWMO_INPUTAVAILABLE);
                if (wr == WAIT_TIMEOUT) break;
                while (PeekMessageW(&msg, nullptr, 0, 0, PM_REMOVE)) {
                    if (msg.message == WM_QUIT) { done = true; break; }
                    TranslateMessage(&msg);
                    DispatchMessageW(&msg);
                }
            }
        }

        // Ordered teardown: close + release the COM objects (controller/webview ComPtrs
        // destruct at scope exit, BEFORE coGuard's CoUninitialize) and destroy the host.
        if (controller) controller->Close();
        if (host) DestroyWindow(host);
    });
    worker.join();
    return result;
}

}  // namespace stimlab::agent

#else  // !STIMLAB_HAVE_WEBVIEW2 - stubs so web_fetch transparently uses the curl path.

namespace stimlab::agent {

bool webViewRenderAvailable() { return false; }

std::optional<std::string> webFetchRenderedHtml(const std::string&, int) { return std::nullopt; }

}  // namespace stimlab::agent

#endif  // STIMLAB_HAVE_WEBVIEW2
