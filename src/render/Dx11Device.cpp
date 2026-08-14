#include "render/Dx11Device.h"

#ifndef WIN32_LEAN_AND_MEAN
#  define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <d3d11.h>
#include <wincodec.h>

#include <cstring>
#include <vector>

#include <spdlog/spdlog.h>

#pragma comment(lib, "d3d11.lib")
#pragma comment(lib, "dxgi.lib")
#pragma comment(lib, "windowscodecs.lib")

namespace biocad {

Dx11Device::~Dx11Device() { shutdown(); }

bool Dx11Device::init(HWND hwnd) {
    DXGI_SWAP_CHAIN_DESC sd{};
    sd.BufferCount = 2;
    sd.BufferDesc.Width = 0;
    sd.BufferDesc.Height = 0;
    sd.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    sd.BufferDesc.RefreshRate.Numerator = 60;
    sd.BufferDesc.RefreshRate.Denominator = 1;
    sd.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
    sd.OutputWindow = hwnd;
    sd.SampleDesc.Count = 1;
    sd.SampleDesc.Quality = 0;
    sd.Windowed = TRUE;
    sd.SwapEffect = DXGI_SWAP_EFFECT_DISCARD;

    UINT flags = 0;
#if defined(_DEBUG)
    flags |= D3D11_CREATE_DEVICE_DEBUG;
#endif
    const D3D_FEATURE_LEVEL levels[] = {D3D_FEATURE_LEVEL_11_0, D3D_FEATURE_LEVEL_10_0};
    D3D_FEATURE_LEVEL obtained{};

    HRESULT hr = D3D11CreateDeviceAndSwapChain(
        nullptr, D3D_DRIVER_TYPE_HARDWARE, nullptr, flags, levels,
        static_cast<UINT>(sizeof(levels) / sizeof(levels[0])), D3D11_SDK_VERSION,
        &sd, &swapChain_, &device_, &obtained, &context_);

    if (hr == DXGI_ERROR_UNSUPPORTED) {
        // Fall back to the WARP software rasterizer (e.g. headless/CI).
        hr = D3D11CreateDeviceAndSwapChain(
            nullptr, D3D_DRIVER_TYPE_WARP, nullptr, flags, levels,
            static_cast<UINT>(sizeof(levels) / sizeof(levels[0])), D3D11_SDK_VERSION,
            &sd, &swapChain_, &device_, &obtained, &context_);
    }
    if (FAILED(hr)) return false;

    createRenderTarget();
    return true;
}

void Dx11Device::createRenderTarget() {
    if (!swapChain_) return;
    ID3D11Texture2D* backBuffer = nullptr;
    if (SUCCEEDED(swapChain_->GetBuffer(0, IID_PPV_ARGS(&backBuffer))) && backBuffer) {
        device_->CreateRenderTargetView(backBuffer, nullptr, &renderTarget_);
        backBuffer->Release();
    }
}

void Dx11Device::releaseRenderTarget() {
    if (renderTarget_) {
        renderTarget_->Release();
        renderTarget_ = nullptr;
    }
}

void Dx11Device::resize(unsigned width, unsigned height) {
    if (!swapChain_) return;
    releaseRenderTarget();
    swapChain_->ResizeBuffers(0, width, height, DXGI_FORMAT_UNKNOWN, 0);
    createRenderTarget();
}

void Dx11Device::beginFrame(const float clearColor[4]) {
    if (!context_ || !renderTarget_) return;
    context_->OMSetRenderTargets(1, &renderTarget_, nullptr);
    context_->ClearRenderTargetView(renderTarget_, clearColor);
}

void Dx11Device::present(bool vsync) {
    if (swapChain_) swapChain_->Present(vsync ? 1 : 0, 0);
}

namespace {

// Minimal scope guard for the COM interfaces this file touches; the render layer
// deliberately has no WRL/ComPtr dependency.
template <typename T>
struct Released {
    T* p = nullptr;
    ~Released() { if (p) p->Release(); }
    T** put() { return &p; }
    T* operator->() const { return p; }
    explicit operator bool() const { return p != nullptr; }
};

bool fail(const char* what, HRESULT hr) {
    spdlog::error("captureBackBufferPng: {} failed (hr=0x{:08X})", what,
                  static_cast<unsigned>(hr));
    return false;
}

}  // namespace

bool Dx11Device::captureBackBufferPng(const std::filesystem::path& out) {
    if (!swapChain_ || !device_ || !context_) return fail("device", E_POINTER);

    Released<ID3D11Texture2D> back;
    HRESULT hr = swapChain_->GetBuffer(0, IID_PPV_ARGS(back.put()));
    if (FAILED(hr) || !back) return fail("GetBuffer", hr);

    D3D11_TEXTURE2D_DESC desc{};
    back->GetDesc(&desc);
    const UINT width = desc.Width;
    const UINT height = desc.Height;

    // Identical description apart from the CPU-readable staging usage; a staging
    // texture cannot be bound, so BindFlags and MiscFlags must be cleared.
    desc.Usage = D3D11_USAGE_STAGING;
    desc.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
    desc.BindFlags = 0;
    desc.MiscFlags = 0;

    Released<ID3D11Texture2D> staging;
    hr = device_->CreateTexture2D(&desc, nullptr, staging.put());
    if (FAILED(hr) || !staging) return fail("CreateTexture2D(staging)", hr);

    context_->CopyResource(staging.p, back.p);

    D3D11_MAPPED_SUBRESOURCE mapped{};
    hr = context_->Map(staging.p, 0, D3D11_MAP_READ, 0, &mapped);
    if (FAILED(hr)) return fail("Map", hr);

    // The swap chain is DXGI_FORMAT_R8G8B8A8_UNORM, so the WIC pixel format is
    // 32bppRGBA - naming it BGRA here would silently swap red and blue.
    const WICPixelFormatGUID pixelFormat = GUID_WICPixelFormat32bppRGBA;
    const UINT stride = width * 4;
    std::vector<BYTE> pixels(static_cast<std::size_t>(stride) * height);
    const auto* src = static_cast<const BYTE*>(mapped.pData);
    for (UINT y = 0; y < height; ++y) {
        std::memcpy(pixels.data() + static_cast<std::size_t>(y) * stride,
                    src + static_cast<std::size_t>(y) * mapped.RowPitch, stride);
    }
    context_->Unmap(staging.p, 0);

    std::error_code ec;
    if (out.has_parent_path()) std::filesystem::create_directories(out.parent_path(), ec);

    Released<IWICImagingFactory> factory;
    hr = CoCreateInstance(CLSID_WICImagingFactory, nullptr, CLSCTX_INPROC_SERVER,
                          IID_PPV_ARGS(factory.put()));
    if (FAILED(hr) || !factory) return fail("CoCreateInstance(WICImagingFactory)", hr);

    Released<IWICStream> stream;
    hr = factory->CreateStream(stream.put());
    if (FAILED(hr)) return fail("CreateStream", hr);
    hr = stream->InitializeFromFilename(out.wstring().c_str(), GENERIC_WRITE);
    if (FAILED(hr)) return fail("InitializeFromFilename", hr);

    Released<IWICBitmapEncoder> encoder;
    hr = factory->CreateEncoder(GUID_ContainerFormatPng, nullptr, encoder.put());
    if (FAILED(hr)) return fail("CreateEncoder", hr);
    hr = encoder->Initialize(stream.p, WICBitmapEncoderNoCache);
    if (FAILED(hr)) return fail("encoder Initialize", hr);

    Released<IWICBitmapFrameEncode> frame;
    IPropertyBag2* props = nullptr;
    hr = encoder->CreateNewFrame(frame.put(), &props);
    if (props) props->Release();
    if (FAILED(hr)) return fail("CreateNewFrame", hr);
    hr = frame->Initialize(nullptr);
    if (FAILED(hr)) return fail("frame Initialize", hr);
    hr = frame->SetSize(width, height);
    if (FAILED(hr)) return fail("SetSize", hr);

    // Wrap the pixels as a WIC bitmap, then let the encoder negotiate its preferred
    // format: the built-in PNG codec does not accept 32bppRGBA as-is, and the right
    // response to negotiation is a colour-correct CONVERSION, not a refusal.
    Released<IWICBitmap> bitmap;
    hr = factory->CreateBitmapFromMemory(width, height, pixelFormat, stride,
                                         static_cast<UINT>(pixels.size()), pixels.data(),
                                         bitmap.put());
    if (FAILED(hr)) return fail("CreateBitmapFromMemory", hr);

    WICPixelFormatGUID negotiated = pixelFormat;
    hr = frame->SetPixelFormat(&negotiated);
    if (FAILED(hr)) return fail("SetPixelFormat", hr);

    IWICBitmapSource* source = bitmap.p;
    Released<IWICFormatConverter> converter;
    if (!IsEqualGUID(negotiated, pixelFormat)) {
        hr = factory->CreateFormatConverter(converter.put());
        if (FAILED(hr)) return fail("CreateFormatConverter", hr);
        hr = converter->Initialize(bitmap.p, negotiated, WICBitmapDitherTypeNone,
                                   nullptr, 0.0, WICBitmapPaletteTypeCustom);
        if (FAILED(hr)) return fail("converter Initialize", hr);
        source = converter.p;
    }

    hr = frame->WriteSource(source, nullptr);
    if (FAILED(hr)) return fail("WriteSource", hr);
    hr = frame->Commit();
    if (FAILED(hr)) return fail("frame Commit", hr);
    hr = encoder->Commit();
    if (FAILED(hr)) return fail("encoder Commit", hr);

    spdlog::info("captured {}x{} backbuffer to {}", width, height, out.string());
    return true;
}

void Dx11Device::shutdown() {
    releaseRenderTarget();
    if (swapChain_) { swapChain_->Release(); swapChain_ = nullptr; }
    if (context_)   { context_->Release();   context_ = nullptr; }
    if (device_)    { device_->Release();    device_ = nullptr; }
}

}  // namespace biocad
