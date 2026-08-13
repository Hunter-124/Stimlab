// render/Dx11Device.h - minimal DirectX 11 device + swapchain wrapper.
// Owns the device/context/swapchain/backbuffer RTV and handles resize + present.
#pragma once

struct ID3D11Device;
struct ID3D11DeviceContext;
struct IDXGISwapChain;
struct ID3D11RenderTargetView;
struct HWND__;
using HWND = HWND__*;

namespace biocad {

class Dx11Device {
public:
    Dx11Device() = default;
    ~Dx11Device();
    Dx11Device(const Dx11Device&) = delete;
    Dx11Device& operator=(const Dx11Device&) = delete;

    bool init(HWND hwnd);
    void shutdown();

    void resize(unsigned width, unsigned height);  // 0,0 = pull from swapchain
    void beginFrame(const float clearColor[4]);
    void present(bool vsync);

    [[nodiscard]] ID3D11Device* device() const { return device_; }
    [[nodiscard]] ID3D11DeviceContext* context() const { return context_; }
    [[nodiscard]] bool valid() const { return device_ != nullptr; }

private:
    void createRenderTarget();
    void releaseRenderTarget();

    ID3D11Device* device_ = nullptr;
    ID3D11DeviceContext* context_ = nullptr;
    IDXGISwapChain* swapChain_ = nullptr;
    ID3D11RenderTargetView* renderTarget_ = nullptr;
};

}  // namespace biocad
