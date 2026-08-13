// render/MolViewport.cpp - DX11 ball-and-stick / spacefill molecular renderer.
//
// Approach: instanced low-poly meshes (a UV sphere per atom, a unit cylinder per
// bond) drawn with per-instance transform + color into an off-screen RGBA target.
// This is the "acceptable alternative" path from the WP-2 brief; it is robust and
// trivially holds 60fps for the library's <=60-atom molecules. Lighting is a
// simple Lambert + ambient + rim in the pixel shader. The cylinder is split-colored
// by interpolating between the two endpoint colors across its length so a C-N bond
// shows grey on the carbon half and blue on the nitrogen half.
#include "render/MolViewport.h"

#ifndef WIN32_LEAN_AND_MEAN
#  define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <d3d11.h>
#include <d3dcompiler.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#pragma comment(lib, "d3d11.lib")
#pragma comment(lib, "d3dcompiler.lib")

namespace biocad::render {
namespace {

// --------------------------------------------------------------------------
// Tiny math (column-vector convention, row-major storage; we build matrices
// to be used as mul(float4(pos,1), M) in HLSL i.e. row-vector * matrix).
// --------------------------------------------------------------------------
struct V3 { float x, y, z; };

inline V3 sub(const V3& a, const V3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline V3 cross(const V3& a, const V3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline float dot(const V3& a, const V3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline float length(const V3& a) { return std::sqrt(dot(a, a)); }
inline V3 normalize(const V3& a) {
    const float l = length(a);
    return l > 1e-8f ? V3{a.x / l, a.y / l, a.z / l} : V3{0, 0, 1};
}

// Row-major 4x4. m[row][col]. HLSL receives it as a row_major matrix and we do
// mul(rowVector, M).
struct M4 {
    float m[4][4];
};

M4 identity() {
    M4 r{};
    for (int i = 0; i < 4; ++i) r.m[i][i] = 1.0f;
    return r;
}

// out = a * b (row-vector convention: v' = v*a*b applies a first).
M4 mul(const M4& a, const M4& b) {
    M4 r{};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            float s = 0;
            for (int k = 0; k < 4; ++k) s += a.m[i][k] * b.m[k][j];
            r.m[i][j] = s;
        }
    return r;
}

// Right-handed look-at producing a view matrix for row-vector math.
M4 lookAtRH(const V3& eye, const V3& at, const V3& up) {
    const V3 zaxis = normalize(sub(eye, at));     // camera looks down -z
    const V3 xaxis = normalize(cross(up, zaxis));
    const V3 yaxis = cross(zaxis, xaxis);
    M4 r = identity();
    r.m[0][0] = xaxis.x; r.m[1][0] = xaxis.y; r.m[2][0] = xaxis.z;
    r.m[0][1] = yaxis.x; r.m[1][1] = yaxis.y; r.m[2][1] = yaxis.z;
    r.m[0][2] = zaxis.x; r.m[1][2] = zaxis.y; r.m[2][2] = zaxis.z;
    r.m[3][0] = -dot(xaxis, eye);
    r.m[3][1] = -dot(yaxis, eye);
    r.m[3][2] = -dot(zaxis, eye);
    return r;
}

// Right-handed perspective, depth 0..1 (D3D), for row-vector math.
M4 perspectiveRH(float fovY, float aspect, float zn, float zf) {
    const float yScale = 1.0f / std::tan(fovY * 0.5f);
    const float xScale = yScale / aspect;
    M4 r{};
    r.m[0][0] = xScale;
    r.m[1][1] = yScale;
    r.m[2][2] = zf / (zn - zf);
    r.m[2][3] = -1.0f;
    r.m[3][2] = zn * zf / (zn - zf);
    return r;
}

// --------------------------------------------------------------------------
// Mesh generation.
// --------------------------------------------------------------------------
struct MeshVertex {
    float px, py, pz;   // position
    float nx, ny, nz;   // normal
    float t;            // axial parameter 0..1 (used by the cylinder for split color)
};

void makeUvSphere(int stacks, int slices, std::vector<MeshVertex>& verts,
                  std::vector<unsigned>& idx) {
    verts.clear();
    idx.clear();
    for (int i = 0; i <= stacks; ++i) {
        const float v = static_cast<float>(i) / stacks;
        const float phi = v * 3.14159265358979f;       // 0..pi
        const float sp = std::sin(phi), cp = std::cos(phi);
        for (int j = 0; j <= slices; ++j) {
            const float u = static_cast<float>(j) / slices;
            const float theta = u * 6.28318530717959f;  // 0..2pi
            const float st = std::sin(theta), ct = std::cos(theta);
            const float x = sp * ct, y = cp, z = sp * st;
            verts.push_back({x, y, z, x, y, z, 0.0f});  // unit sphere: pos == normal
        }
    }
    const int ring = slices + 1;
    for (int i = 0; i < stacks; ++i) {
        for (int j = 0; j < slices; ++j) {
            const unsigned a = static_cast<unsigned>(i * ring + j);
            const unsigned b = static_cast<unsigned>((i + 1) * ring + j);
            idx.push_back(a); idx.push_back(b); idx.push_back(a + 1);
            idx.push_back(a + 1); idx.push_back(b); idx.push_back(b + 1);
        }
    }
}

// Unit cylinder along +Y, radius 1, from y=0 (t=0) to y=1 (t=1), with caps.
void makeCylinder(int slices, std::vector<MeshVertex>& verts, std::vector<unsigned>& idx) {
    verts.clear();
    idx.clear();
    // Side wall: two rings (bottom t=0, top t=1).
    for (int ring = 0; ring < 2; ++ring) {
        const float y = static_cast<float>(ring);
        for (int j = 0; j <= slices; ++j) {
            const float u = static_cast<float>(j) / slices;
            const float theta = u * 6.28318530717959f;
            const float ct = std::cos(theta), st = std::sin(theta);
            verts.push_back({ct, y, st, ct, 0.0f, st, y});
        }
    }
    const int ring = slices + 1;
    for (int j = 0; j < slices; ++j) {
        const unsigned a = static_cast<unsigned>(j);
        const unsigned b = static_cast<unsigned>(ring + j);
        idx.push_back(a); idx.push_back(b); idx.push_back(a + 1);
        idx.push_back(a + 1); idx.push_back(b); idx.push_back(b + 1);
    }
    // Caps (flat). Bottom center, top center.
    const unsigned baseBottom = static_cast<unsigned>(verts.size());
    verts.push_back({0, 0, 0, 0, -1, 0, 0});
    const unsigned centerBottom = baseBottom;
    const unsigned bottomStart = static_cast<unsigned>(verts.size());
    for (int j = 0; j <= slices; ++j) {
        const float u = static_cast<float>(j) / slices;
        const float theta = u * 6.28318530717959f;
        verts.push_back({std::cos(theta), 0.0f, std::sin(theta), 0, -1, 0, 0});
    }
    for (int j = 0; j < slices; ++j) {
        idx.push_back(centerBottom);
        idx.push_back(bottomStart + static_cast<unsigned>(j) + 1);
        idx.push_back(bottomStart + static_cast<unsigned>(j));
    }
    const unsigned centerTop = static_cast<unsigned>(verts.size());
    verts.push_back({0, 1, 0, 0, 1, 0, 1});
    const unsigned topStart = static_cast<unsigned>(verts.size());
    for (int j = 0; j <= slices; ++j) {
        const float u = static_cast<float>(j) / slices;
        const float theta = u * 6.28318530717959f;
        verts.push_back({std::cos(theta), 1.0f, std::sin(theta), 0, 1, 0, 1});
    }
    for (int j = 0; j < slices; ++j) {
        idx.push_back(centerTop);
        idx.push_back(topStart + static_cast<unsigned>(j));
        idx.push_back(topStart + static_cast<unsigned>(j) + 1);
    }
}

// --------------------------------------------------------------------------
// Per-instance GPU layout. Mirrors AtomInst / a transformed BondInst.
// --------------------------------------------------------------------------
struct GpuAtomInst {
    float cx, cy, cz, radius;     // center + radius
    float r, g, b, a;             // color
};

// For bonds we precompute an orientation basis on the CPU so the vertex shader
// is a cheap basis transform. The cylinder model's +Y maps to `dir`, and X/Z to
// `right`/`fwd`; length scales Y, radius scales X/Z.
struct GpuBondInst {
    float ox, oy, oz, length;     // origin (atom A) + length
    float dx, dy, dz, radius;     // axis direction (unit) + radius
    float rx, ry, rz, _pad0;      // 'right' basis vector
    float fx, fy, fz, _pad1;      // 'forward' basis vector
    float ca_r, ca_g, ca_b, ca_a; // color at A end
    float cb_r, cb_g, cb_b, cb_a; // color at B end
};

struct CBuffer {
    M4 viewProj;
    float lightDir[4];   // xyz + pad
    float camPos[4];
};

// --------------------------------------------------------------------------
// HLSL (compiled at runtime via D3DCompile). One shared pixel shader style.
// --------------------------------------------------------------------------
const char* kAtomShader = R"hlsl(
cbuffer CB : register(b0) {
    row_major float4x4 gViewProj;
    float4 gLightDir;
    float4 gCamPos;
};
struct VSIn {
    float3 pos : POSITION;
    float3 nrm : NORMAL;
    float  t   : TEXCOORD0;
    float4 cr  : INST_CR;   // center.xyz + radius
    float4 col : INST_COL;  // rgba
};
struct VSOut {
    float4 pos : SV_POSITION;
    float3 nrm : NORMAL;
    float3 wpos: TEXCOORD0;
    float4 col : COLOR0;
};
VSOut VSMain(VSIn i) {
    VSOut o;
    float3 world = i.cr.xyz + i.pos * i.cr.w;
    o.pos  = mul(float4(world, 1.0), gViewProj);
    o.nrm  = i.nrm;            // unit sphere: model normal == world normal
    o.wpos = world;
    o.col  = i.col;
    return o;
}
float4 PSMain(VSOut i) : SV_Target {
    float3 N = normalize(i.nrm);
    float3 L = normalize(-gLightDir.xyz);
    float3 Vd = normalize(gCamPos.xyz - i.wpos);
    float diff = saturate(dot(N, L));
    float3 H = normalize(L + Vd);
    float spec = pow(saturate(dot(N, H)), 32.0) * 0.35;
    float rim = pow(1.0 - saturate(dot(N, Vd)), 3.0) * 0.20;
    float amb = 0.30;
    float3 c = i.col.rgb * (amb + diff * 0.85) + spec.xxx + rim.xxx;
    return float4(saturate(c), 1.0);
}
)hlsl";

const char* kBondShader = R"hlsl(
cbuffer CB : register(b0) {
    row_major float4x4 gViewProj;
    float4 gLightDir;
    float4 gCamPos;
};
struct VSIn {
    float3 pos : POSITION;
    float3 nrm : NORMAL;
    float  t   : TEXCOORD0;
    float4 orL : INST_ORL;   // origin.xyz + length
    float4 dirR: INST_DIRR;  // dir.xyz + radius
    float4 right: INST_RIGHT;
    float4 fwd : INST_FWD;
    float4 cA  : INST_CA;
    float4 cB  : INST_CB;
};
struct VSOut {
    float4 pos : SV_POSITION;
    float3 nrm : NORMAL;
    float3 wpos: TEXCOORD0;
    float  t   : TEXCOORD1;
    float4 cA  : COLOR0;
    float4 cB  : COLOR1;
};
VSOut VSMain(VSIn i) {
    VSOut o;
    // Model: X/Z scaled by radius, Y scaled by length, then rotate into basis.
    float3 local = float3(i.pos.x * i.dirR.w, i.pos.y * i.orL.w, i.pos.z * i.dirR.w);
    float3 world = i.orL.xyz
                 + i.right.xyz * local.x
                 + i.dirR.xyz  * local.y
                 + i.fwd.xyz   * local.z;
    float3 nrm = i.nrm.x * i.right.xyz + i.nrm.y * i.dirR.xyz + i.nrm.z * i.fwd.xyz;
    o.pos  = mul(float4(world, 1.0), gViewProj);
    o.nrm  = nrm;
    o.wpos = world;
    o.t    = i.t;
    o.cA   = i.cA;
    o.cB   = i.cB;
    return o;
}
float4 PSMain(VSOut i) : SV_Target {
    float3 N = normalize(i.nrm);
    float3 L = normalize(-gLightDir.xyz);
    float3 Vd = normalize(gCamPos.xyz - i.wpos);
    float diff = saturate(dot(N, L));
    float3 H = normalize(L + Vd);
    float spec = pow(saturate(dot(N, H)), 32.0) * 0.30;
    float amb = 0.30;
    float4 base = (i.t < 0.5) ? i.cA : i.cB;  // hard split at the bond midpoint
    float3 c = base.rgb * (amb + diff * 0.85) + spec.xxx;
    return float4(saturate(c), 1.0);
}
)hlsl";

ID3DBlob* compile(const char* src, const char* entry, const char* target) {
    ID3DBlob* blob = nullptr;
    ID3DBlob* err = nullptr;
    UINT flags = D3DCOMPILE_ENABLE_STRICTNESS | D3DCOMPILE_OPTIMIZATION_LEVEL3;
    HRESULT hr = D3DCompile(src, std::strlen(src), nullptr, nullptr, nullptr, entry, target,
                            flags, 0, &blob, &err);
    if (err) {
        OutputDebugStringA(static_cast<const char*>(err->GetBufferPointer()));
        err->Release();
    }
    if (FAILED(hr)) {
        if (blob) blob->Release();
        return nullptr;
    }
    return blob;
}

template <class T>
void safeRelease(T*& p) {
    if (p) { p->Release(); p = nullptr; }
}

inline void unpackColor(std::uint32_t rgba, float out[4]) {
    out[0] = static_cast<float>(rgba & 0xFF) / 255.0f;           // R
    out[1] = static_cast<float>((rgba >> 8) & 0xFF) / 255.0f;    // G
    out[2] = static_cast<float>((rgba >> 16) & 0xFF) / 255.0f;   // B
    out[3] = static_cast<float>((rgba >> 24) & 0xFF) / 255.0f;   // A
}

}  // namespace

// ---------------------------------------------------------------------------
// CPU scene building.
// ---------------------------------------------------------------------------

// Pack 0xRRGGBB (CPK literal) into 0xAABBGGRR with full alpha.
static std::uint32_t packRGB(std::uint32_t rrggbb) {
    const std::uint32_t r = (rrggbb >> 16) & 0xFF;
    const std::uint32_t g = (rrggbb >> 8) & 0xFF;
    const std::uint32_t b = rrggbb & 0xFF;
    return 0xFF000000u | (b << 16) | (g << 8) | r;
}

std::uint32_t cpkColor(int z) {
    switch (z) {
        case 1:  return packRGB(0xFFFFFF);  // H  white
        case 6:  return packRGB(0x909090);  // C  grey
        case 7:  return packRGB(0x3050F8);  // N  blue
        case 8:  return packRGB(0xFF0D0D);  // O  red
        case 9:  return packRGB(0x90E050);  // F  green
        case 15: return packRGB(0xFF8000);  // P  orange
        case 16: return packRGB(0xFFFF30);  // S  yellow
        case 17: return packRGB(0x1FF01F);  // Cl green
        case 35: return packRGB(0xA62929);  // Br dark red
        case 53: return packRGB(0x940094);  // I  violet
        default: return packRGB(0xFFB5B5);  // default pink
    }
}

MolScene buildMolScene(const chem::Conformer& conf, bool spacefill, bool showH) {
    MolScene scene;
    if (conf.empty()) {
        scene.center = {0, 0, 0};
        scene.radius = 1.0f;
        return scene;
    }

    const int n = conf.size();
    // Map original atom index -> emitted (kept) status, so bonds can be filtered.
    std::vector<bool> kept(static_cast<size_t>(n), false);

    // Atoms.
    for (int i = 0; i < n; ++i) {
        const int z = conf.z[i];
        if (!showH && z == 1) continue;
        kept[static_cast<size_t>(i)] = true;
        AtomInst a;
        a.x = static_cast<float>(conf.pos[i].x);
        a.y = static_cast<float>(conf.pos[i].y);
        a.z = static_cast<float>(conf.pos[i].z);
        if (spacefill) {
            a.r = static_cast<float>(chem::vdwRadius(z));
        } else {
            // Ball-and-stick: shrink to a fraction of the covalent radius; keep
            // hydrogens visibly smaller than heavy atoms.
            const double cov = chem::covalentRadius(z);
            a.r = static_cast<float>(cov * (z == 1 ? 0.42 : 0.42));
            // Floor so tiny atoms remain pickable.
            a.r = std::max(a.r, 0.25f);
        }
        a.rgba = cpkColor(z);
        scene.atoms.push_back(a);
    }

    // Bonds (only in ball-and-stick; spacefill hides them).
    if (!spacefill) {
        for (const auto& bp : conf.bonds) {
            const int ia = bp.first, ib = bp.second;
            if (ia < 0 || ib < 0 || ia >= n || ib >= n) continue;
            if (!kept[static_cast<size_t>(ia)] || !kept[static_cast<size_t>(ib)]) continue;
            BondInst b;
            b.ax = static_cast<float>(conf.pos[ia].x);
            b.ay = static_cast<float>(conf.pos[ia].y);
            b.az = static_cast<float>(conf.pos[ia].z);
            b.bx = static_cast<float>(conf.pos[ib].x);
            b.by = static_cast<float>(conf.pos[ib].y);
            b.bz = static_cast<float>(conf.pos[ib].z);
            b.rgbaA = cpkColor(conf.z[ia]);
            b.rgbaB = cpkColor(conf.z[ib]);
            b.radius = 0.14f;
            scene.bonds.push_back(b);
        }
    }

    // Bounding sphere over the KEPT atoms (so hiding H reframes correctly).
    double cx = 0, cy = 0, cz = 0;
    int cnt = 0;
    for (int i = 0; i < n; ++i) {
        if (!kept[static_cast<size_t>(i)]) continue;
        cx += conf.pos[i].x; cy += conf.pos[i].y; cz += conf.pos[i].z;
        ++cnt;
    }
    if (cnt == 0) { scene.center = {0, 0, 0}; scene.radius = 1.0f; return scene; }
    cx /= cnt; cy /= cnt; cz /= cnt;
    scene.center = {cx, cy, cz};
    double maxR = 0.0;
    for (size_t i = 0; i < scene.atoms.size(); ++i) {
        const double dx = scene.atoms[i].x - cx;
        const double dy = scene.atoms[i].y - cy;
        const double dz = scene.atoms[i].z - cz;
        const double d = std::sqrt(dx * dx + dy * dy + dz * dz) + scene.atoms[i].r;
        maxR = std::max(maxR, d);
    }
    scene.radius = std::max(static_cast<float>(maxR), 1.0f);
    return scene;
}

// ---------------------------------------------------------------------------
// MolViewport.
// ---------------------------------------------------------------------------

MolViewport::MolViewport(ID3D11Device* device, ID3D11DeviceContext* context)
    : dev_(device), ctx_(context) {
    ready_ = (dev_ && ctx_) ? createDeviceResources() : false;
}

MolViewport::~MolViewport() { releaseAll(); }

bool MolViewport::createDeviceResources() {
    // --- Geometry ---
    std::vector<MeshVertex> sv, cv;
    std::vector<unsigned> si, ci;
    makeUvSphere(16, 24, sv, si);
    makeCylinder(20, cv, ci);
    sphereIndexCount_ = static_cast<unsigned>(si.size());
    cylIndexCount_ = static_cast<unsigned>(ci.size());

    auto makeBuffer = [&](const void* data, UINT bytes, UINT bind, ID3D11Buffer** out,
                          bool dynamic) -> bool {
        D3D11_BUFFER_DESC bd{};
        bd.ByteWidth = bytes;
        bd.BindFlags = bind;
        if (dynamic) {
            bd.Usage = D3D11_USAGE_DYNAMIC;
            bd.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
        } else {
            bd.Usage = D3D11_USAGE_IMMUTABLE;
        }
        D3D11_SUBRESOURCE_DATA srd{};
        srd.pSysMem = data;
        return SUCCEEDED(dev_->CreateBuffer(&bd, data ? &srd : nullptr, out));
    };

    if (!makeBuffer(sv.data(), static_cast<UINT>(sv.size() * sizeof(MeshVertex)),
                    D3D11_BIND_VERTEX_BUFFER, &sphereVB_, false)) return false;
    if (!makeBuffer(si.data(), static_cast<UINT>(si.size() * sizeof(unsigned)),
                    D3D11_BIND_INDEX_BUFFER, &sphereIB_, false)) return false;
    if (!makeBuffer(cv.data(), static_cast<UINT>(cv.size() * sizeof(MeshVertex)),
                    D3D11_BIND_VERTEX_BUFFER, &cylVB_, false)) return false;
    if (!makeBuffer(ci.data(), static_cast<UINT>(ci.size() * sizeof(unsigned)),
                    D3D11_BIND_INDEX_BUFFER, &cylIB_, false)) return false;

    // --- Constant buffer ---
    {
        D3D11_BUFFER_DESC bd{};
        bd.ByteWidth = sizeof(CBuffer);
        bd.Usage = D3D11_USAGE_DYNAMIC;
        bd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
        bd.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
        if (FAILED(dev_->CreateBuffer(&bd, nullptr, &cbuf_))) return false;
    }

    // --- Shaders ---
    ID3DBlob* avs = compile(kAtomShader, "VSMain", "vs_5_0");
    ID3DBlob* aps = compile(kAtomShader, "PSMain", "ps_5_0");
    ID3DBlob* bvs = compile(kBondShader, "VSMain", "vs_5_0");
    ID3DBlob* bps = compile(kBondShader, "PSMain", "ps_5_0");
    bool ok = avs && aps && bvs && bps;
    if (ok) {
        ok = SUCCEEDED(dev_->CreateVertexShader(avs->GetBufferPointer(), avs->GetBufferSize(),
                                                nullptr, &atomVS_)) &&
             SUCCEEDED(dev_->CreatePixelShader(aps->GetBufferPointer(), aps->GetBufferSize(),
                                               nullptr, &atomPS_)) &&
             SUCCEEDED(dev_->CreateVertexShader(bvs->GetBufferPointer(), bvs->GetBufferSize(),
                                                nullptr, &bondVS_)) &&
             SUCCEEDED(dev_->CreatePixelShader(bps->GetBufferPointer(), bps->GetBufferSize(),
                                               nullptr, &bondPS_));
    }

    if (ok) {
        const D3D11_INPUT_ELEMENT_DESC atomElems[] = {
            {"POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0},
            {"NORMAL", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D11_INPUT_PER_VERTEX_DATA, 0},
            {"TEXCOORD", 0, DXGI_FORMAT_R32_FLOAT, 0, 24, D3D11_INPUT_PER_VERTEX_DATA, 0},
            {"INST_CR", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 1, 0, D3D11_INPUT_PER_INSTANCE_DATA, 1},
            {"INST_COL", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 1, 16, D3D11_INPUT_PER_INSTANCE_DATA, 1},
        };
        ok = SUCCEEDED(dev_->CreateInputLayout(atomElems, 5, avs->GetBufferPointer(),
                                               avs->GetBufferSize(), &atomLayout_));
    }
    if (ok) {
        const D3D11_INPUT_ELEMENT_DESC bondElems[] = {
            {"POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0},
            {"NORMAL", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D11_INPUT_PER_VERTEX_DATA, 0},
            {"TEXCOORD", 0, DXGI_FORMAT_R32_FLOAT, 0, 24, D3D11_INPUT_PER_VERTEX_DATA, 0},
            {"INST_ORL", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 1, 0, D3D11_INPUT_PER_INSTANCE_DATA, 1},
            {"INST_DIRR", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 1, 16, D3D11_INPUT_PER_INSTANCE_DATA, 1},
            {"INST_RIGHT", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 1, 32, D3D11_INPUT_PER_INSTANCE_DATA, 1},
            {"INST_FWD", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 1, 48, D3D11_INPUT_PER_INSTANCE_DATA, 1},
            {"INST_CA", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 1, 64, D3D11_INPUT_PER_INSTANCE_DATA, 1},
            {"INST_CB", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 1, 80, D3D11_INPUT_PER_INSTANCE_DATA, 1},
        };
        ok = SUCCEEDED(dev_->CreateInputLayout(bondElems, 9, bvs->GetBufferPointer(),
                                               bvs->GetBufferSize(), &bondLayout_));
    }

    safeRelease(avs); safeRelease(aps); safeRelease(bvs); safeRelease(bps);
    if (!ok) return false;

    // --- Render states ---
    {
        D3D11_RASTERIZER_DESC rd{};
        rd.FillMode = D3D11_FILL_SOLID;
        rd.CullMode = D3D11_CULL_BACK;
        rd.FrontCounterClockwise = FALSE;
        rd.DepthClipEnable = TRUE;
        if (FAILED(dev_->CreateRasterizerState(&rd, &raster_))) return false;
    }
    {
        D3D11_DEPTH_STENCIL_DESC dd{};
        dd.DepthEnable = TRUE;
        dd.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ALL;
        dd.DepthFunc = D3D11_COMPARISON_LESS;
        if (FAILED(dev_->CreateDepthStencilState(&dd, &depthState_))) return false;
    }
    return true;
}

bool MolViewport::ensureTarget(int width, int height) {
    if (width <= 0 || height <= 0) return false;
    if (colorTex_ && width == rtW_ && height == rtH_) return true;
    releaseTarget();

    D3D11_TEXTURE2D_DESC td{};
    td.Width = static_cast<UINT>(width);
    td.Height = static_cast<UINT>(height);
    td.MipLevels = 1;
    td.ArraySize = 1;
    td.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    td.SampleDesc.Count = 1;
    td.Usage = D3D11_USAGE_DEFAULT;
    td.BindFlags = D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE;
    if (FAILED(dev_->CreateTexture2D(&td, nullptr, &colorTex_))) return false;
    if (FAILED(dev_->CreateRenderTargetView(colorTex_, nullptr, &colorRtv_))) return false;
    if (FAILED(dev_->CreateShaderResourceView(colorTex_, nullptr, &colorSrv_))) return false;

    D3D11_TEXTURE2D_DESC dd = td;
    dd.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;
    dd.BindFlags = D3D11_BIND_DEPTH_STENCIL;
    if (FAILED(dev_->CreateTexture2D(&dd, nullptr, &depthTex_))) return false;
    if (FAILED(dev_->CreateDepthStencilView(depthTex_, nullptr, &depthDsv_))) return false;

    rtW_ = width;
    rtH_ = height;
    return true;
}

bool MolViewport::ensureAtomInstCapacity(unsigned count) {
    if (count == 0) return true;
    if (atomInstVB_ && count <= atomInstCap_) return true;
    safeRelease(atomInstVB_);
    const unsigned cap = std::max<unsigned>(count, atomInstCap_ ? atomInstCap_ * 2 : 64);
    D3D11_BUFFER_DESC bd{};
    bd.ByteWidth = cap * sizeof(GpuAtomInst);
    bd.Usage = D3D11_USAGE_DYNAMIC;
    bd.BindFlags = D3D11_BIND_VERTEX_BUFFER;
    bd.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
    if (FAILED(dev_->CreateBuffer(&bd, nullptr, &atomInstVB_))) { atomInstCap_ = 0; return false; }
    atomInstCap_ = cap;
    return true;
}

bool MolViewport::ensureBondInstCapacity(unsigned count) {
    if (count == 0) return true;
    if (bondInstVB_ && count <= bondInstCap_) return true;
    safeRelease(bondInstVB_);
    const unsigned cap = std::max<unsigned>(count, bondInstCap_ ? bondInstCap_ * 2 : 64);
    D3D11_BUFFER_DESC bd{};
    bd.ByteWidth = cap * sizeof(GpuBondInst);
    bd.Usage = D3D11_USAGE_DYNAMIC;
    bd.BindFlags = D3D11_BIND_VERTEX_BUFFER;
    bd.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
    if (FAILED(dev_->CreateBuffer(&bd, nullptr, &bondInstVB_))) { bondInstCap_ = 0; return false; }
    bondInstCap_ = cap;
    return true;
}

void MolViewport::resetCamera(const MolScene& scene) {
    target_ = scene.center;
    sceneRadius_ = std::max(scene.radius, 0.5f);
    // Distance to fit the bounding sphere in a ~45deg vertical FOV with margin.
    fitDist_ = sceneRadius_ / std::tan(0.5f * 0.785398f) * 1.25f;
    dist_ = fitDist_;
    yaw_ = 0.6f;
    pitch_ = 0.30f;
}

void MolViewport::orbit(float dxPixels, float dyPixels) {
    constexpr float kSpeed = 0.010f;
    yaw_ -= dxPixels * kSpeed;
    pitch_ += dyPixels * kSpeed;
    const float lim = 1.55334f;  // ~89deg, avoid gimbal at the poles
    pitch_ = std::clamp(pitch_, -lim, lim);
}

void MolViewport::zoom(float wheelDelta) {
    const float factor = std::pow(0.88f, wheelDelta);
    dist_ = std::clamp(dist_ * factor, sceneRadius_ * 0.35f, sceneRadius_ * 12.0f);
}

void MolViewport::pan(float dxPixels, float dyPixels) {
    // Pan in the camera plane, scaled by distance so it feels consistent.
    const float s = dist_ * 0.0018f;
    const float cy = std::cos(yaw_), sy = std::sin(yaw_);
    const float cp = std::cos(pitch_), sp = std::sin(pitch_);
    // Camera basis (matches the eye computation in draw()).
    const V3 fwd{-(cp * sy), -sp, -(cp * cy)};
    const V3 worldUp{0, 1, 0};
    V3 right = normalize(cross(fwd, worldUp));
    V3 up = cross(right, fwd);
    target_.x += (-dxPixels * right.x + dyPixels * up.x) * s;
    target_.y += (-dxPixels * right.y + dyPixels * up.y) * s;
    target_.z += (-dxPixels * right.z + dyPixels * up.z) * s;
}

void* MolViewport::draw(const MolScene& scene, int width, int height) {
    if (!ready_) return nullptr;
    if (!ensureTarget(width, height)) return nullptr;

    // --- Camera ---
    const float cy = std::cos(yaw_), sy = std::sin(yaw_);
    const float cp = std::cos(pitch_), sp = std::sin(pitch_);
    const V3 tgt{static_cast<float>(target_.x), static_cast<float>(target_.y),
                 static_cast<float>(target_.z)};
    const V3 eye{tgt.x + dist_ * cp * sy, tgt.y + dist_ * sp, tgt.z + dist_ * cp * cy};
    const float aspect = static_cast<float>(width) / static_cast<float>(height);
    const float zn = std::max(dist_ - sceneRadius_ * 2.5f, 0.05f);
    const float zf = dist_ + sceneRadius_ * 4.0f + 1.0f;
    const M4 view = lookAtRH(eye, tgt, V3{0, 1, 0});
    const M4 proj = perspectiveRH(0.785398f, aspect, zn, zf);
    const M4 viewProj = mul(view, proj);

    // --- Update constant buffer ---
    {
        D3D11_MAPPED_SUBRESOURCE ms{};
        if (SUCCEEDED(ctx_->Map(cbuf_, 0, D3D11_MAP_WRITE_DISCARD, 0, &ms))) {
            CBuffer cb{};
            cb.viewProj = viewProj;
            // Light from over the camera's shoulder.
            const V3 ld = normalize(V3{eye.x - tgt.x + 1.0f, eye.y - tgt.y + 2.0f,
                                       eye.z - tgt.z + 0.5f});
            cb.lightDir[0] = -ld.x; cb.lightDir[1] = -ld.y; cb.lightDir[2] = -ld.z;
            cb.camPos[0] = eye.x; cb.camPos[1] = eye.y; cb.camPos[2] = eye.z;
            std::memcpy(ms.pData, &cb, sizeof(cb));
            ctx_->Unmap(cbuf_, 0);
        }
    }

    // --- Upload instance data ---
    const unsigned atomCount = static_cast<unsigned>(scene.atoms.size());
    const unsigned bondCount = static_cast<unsigned>(scene.bonds.size());
    if (!ensureAtomInstCapacity(atomCount)) return nullptr;
    if (!ensureBondInstCapacity(bondCount)) return nullptr;

    if (atomCount && atomInstVB_) {
        D3D11_MAPPED_SUBRESOURCE ms{};
        if (SUCCEEDED(ctx_->Map(atomInstVB_, 0, D3D11_MAP_WRITE_DISCARD, 0, &ms))) {
            auto* dst = static_cast<GpuAtomInst*>(ms.pData);
            for (unsigned i = 0; i < atomCount; ++i) {
                const AtomInst& a = scene.atoms[i];
                GpuAtomInst g;
                g.cx = a.x; g.cy = a.y; g.cz = a.z; g.radius = a.r;
                unpackColor(a.rgba, &g.r);
                dst[i] = g;
            }
            ctx_->Unmap(atomInstVB_, 0);
        }
    }
    if (bondCount && bondInstVB_) {
        D3D11_MAPPED_SUBRESOURCE ms{};
        if (SUCCEEDED(ctx_->Map(bondInstVB_, 0, D3D11_MAP_WRITE_DISCARD, 0, &ms))) {
            auto* dst = static_cast<GpuBondInst*>(ms.pData);
            for (unsigned i = 0; i < bondCount; ++i) {
                const BondInst& b = scene.bonds[i];
                const V3 A{b.ax, b.ay, b.az};
                const V3 B{b.bx, b.by, b.bz};
                const V3 axis = sub(B, A);
                const float len = length(axis);
                const V3 dir = normalize(axis);
                // Build an orthonormal basis around dir.
                V3 ref = (std::fabs(dir.y) < 0.99f) ? V3{0, 1, 0} : V3{1, 0, 0};
                V3 right = normalize(cross(ref, dir));
                V3 fwd = cross(dir, right);
                GpuBondInst g{};
                g.ox = A.x; g.oy = A.y; g.oz = A.z; g.length = len;
                g.dx = dir.x; g.dy = dir.y; g.dz = dir.z; g.radius = b.radius;
                g.rx = right.x; g.ry = right.y; g.rz = right.z;
                g.fx = fwd.x; g.fy = fwd.y; g.fz = fwd.z;
                unpackColor(b.rgbaA, &g.ca_r);
                unpackColor(b.rgbaB, &g.cb_r);
                dst[i] = g;
            }
            ctx_->Unmap(bondInstVB_, 0);
        }
    }

    // --- Render to off-screen target ---
    ID3D11RenderTargetView* prevRtv = nullptr;
    ID3D11DepthStencilView* prevDsv = nullptr;
    ctx_->OMGetRenderTargets(1, &prevRtv, &prevDsv);  // save (we restore below)

    const float clear[4] = {0.055f, 0.07f, 0.10f, 1.0f};
    ctx_->OMSetRenderTargets(1, &colorRtv_, depthDsv_);
    ctx_->ClearRenderTargetView(colorRtv_, clear);
    ctx_->ClearDepthStencilView(depthDsv_, D3D11_CLEAR_DEPTH | D3D11_CLEAR_STENCIL, 1.0f, 0);

    D3D11_VIEWPORT vp{};
    vp.Width = static_cast<float>(width);
    vp.Height = static_cast<float>(height);
    vp.MaxDepth = 1.0f;
    ctx_->RSSetViewports(1, &vp);
    ctx_->RSSetState(raster_);
    ctx_->OMSetDepthStencilState(depthState_, 0);
    ctx_->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
    ctx_->VSSetConstantBuffers(0, 1, &cbuf_);
    ctx_->PSSetConstantBuffers(0, 1, &cbuf_);

    const UINT meshStride = sizeof(MeshVertex);
    const UINT meshOffset = 0;

    // Atoms.
    if (atomCount && atomInstVB_) {
        const UINT instStride = sizeof(GpuAtomInst);
        const UINT instOffset = 0;
        ID3D11Buffer* vbs[2] = {sphereVB_, atomInstVB_};
        const UINT strides[2] = {meshStride, instStride};
        const UINT offsets[2] = {meshOffset, instOffset};
        ctx_->IASetInputLayout(atomLayout_);
        ctx_->IASetVertexBuffers(0, 2, vbs, strides, offsets);
        ctx_->IASetIndexBuffer(sphereIB_, DXGI_FORMAT_R32_UINT, 0);
        ctx_->VSSetShader(atomVS_, nullptr, 0);
        ctx_->PSSetShader(atomPS_, nullptr, 0);
        ctx_->DrawIndexedInstanced(sphereIndexCount_, atomCount, 0, 0, 0);
    }

    // Bonds.
    if (bondCount && bondInstVB_) {
        const UINT instStride = sizeof(GpuBondInst);
        const UINT instOffset = 0;
        ID3D11Buffer* vbs[2] = {cylVB_, bondInstVB_};
        const UINT strides[2] = {meshStride, instStride};
        const UINT offsets[2] = {meshOffset, instOffset};
        ctx_->IASetInputLayout(bondLayout_);
        ctx_->IASetVertexBuffers(0, 2, vbs, strides, offsets);
        ctx_->IASetIndexBuffer(cylIB_, DXGI_FORMAT_R32_UINT, 0);
        ctx_->VSSetShader(bondVS_, nullptr, 0);
        ctx_->PSSetShader(bondPS_, nullptr, 0);
        ctx_->DrawIndexedInstanced(cylIndexCount_, bondCount, 0, 0, 0);
    }

    // Restore whatever render target was bound before (the swapchain backbuffer).
    ctx_->OMSetRenderTargets(1, &prevRtv, prevDsv);
    safeRelease(prevRtv);
    safeRelease(prevDsv);

    return static_cast<void*>(colorSrv_);
}

void MolViewport::releaseTarget() {
    safeRelease(colorSrv_);
    safeRelease(colorRtv_);
    safeRelease(colorTex_);
    safeRelease(depthDsv_);
    safeRelease(depthTex_);
    rtW_ = rtH_ = 0;
}

void MolViewport::releaseAll() {
    releaseTarget();
    safeRelease(sphereVB_); safeRelease(sphereIB_);
    safeRelease(cylVB_); safeRelease(cylIB_);
    safeRelease(atomInstVB_); safeRelease(bondInstVB_);
    safeRelease(cbuf_);
    safeRelease(atomLayout_); safeRelease(bondLayout_);
    safeRelease(atomVS_); safeRelease(atomPS_);
    safeRelease(bondVS_); safeRelease(bondPS_);
    safeRelease(raster_); safeRelease(depthState_);
    atomInstCap_ = bondInstCap_ = 0;
    ready_ = false;
}

}  // namespace biocad::render
