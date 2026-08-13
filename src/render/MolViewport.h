// render/MolViewport.h - DirectX 11 in-app 3D molecular viewer (WP-2, Phase D).
//
// Renders a chem::Conformer as ball-and-stick (atoms = CPK-colored spheres with
// element radii; bonds = split-colored cylinders) into an OFF-SCREEN render
// target that is handed to ImGui via ImGui::Image(). Uses GPU instancing: one
// low-poly UV-sphere mesh drawn per atom and one unit cylinder mesh drawn per
// bond, each with a per-instance transform + color. With the ~60-atom molecules
// in the BioCAD library this trivially holds 60fps.
//
// The CPU scene-building (buildMolScene / cpkColor) is deliberately split out as
// pure functions so it is unit-testable without a D3D device. MolViewport::draw()
// consumes a MolScene and is the only part that touches the GPU.
//
// SAFETY SCOPE: this is a viewer. It depicts WHAT A COMPOUND IS (its 3D shape);
// it produces no synthesis / route / manufacturability content of any kind.
#pragma once

#include <cstdint>
#include <vector>

#include "chem/Embed3D.h"  // chem::Conformer, chem::Vec3, covalent/vdw radii

struct ID3D11Device;
struct ID3D11DeviceContext;
struct ID3D11Texture2D;
struct ID3D11RenderTargetView;
struct ID3D11ShaderResourceView;
struct ID3D11DepthStencilView;
struct ID3D11Buffer;
struct ID3D11VertexShader;
struct ID3D11PixelShader;
struct ID3D11InputLayout;
struct ID3D11RasterizerState;
struct ID3D11DepthStencilState;

namespace biocad::render {

// ----------------------------------------------------------------------------
// CPU-side scene representation (GPU-agnostic, unit-testable).
// ----------------------------------------------------------------------------

// One atom to draw: a sphere centered at (x,y,z) with radius r and an RGBA color
// (0xAABBGGRR byte order so it can be memcpy'd straight into a float4 in-shader).
struct AtomInst {
    float x = 0, y = 0, z = 0;
    float r = 0;
    std::uint32_t rgba = 0xFFFFFFFFu;
};

// One bond to draw: a cylinder from atom A (ax,ay,az) to atom B (bx,by,bz),
// split-colored at the midpoint (rgbaA on A's half, rgbaB on B's half).
struct BondInst {
    float ax = 0, ay = 0, az = 0;
    float bx = 0, by = 0, bz = 0;
    std::uint32_t rgbaA = 0xFFFFFFFFu;
    std::uint32_t rgbaB = 0xFFFFFFFFu;
    float radius = 0.15f;
};

// One vertex of a generic indexed triangle mesh: world-space position, world-space normal and
// a packed colour in the SAME 0xAABBGGRR order as AtomInst::rgba, so the mesh path reuses the
// existing colour packing instead of introducing a second convention.
struct MeshVertexRgba {
    float px = 0, py = 0, pz = 0;
    float nx = 0, ny = 0, nz = 0;
    std::uint32_t rgba = 0xFFFFFFFFu;
};

// A triangle list. This is the path the protein cartoon (bio::CartoonMesh) uses: unlike the
// sphere and cylinder instancing above there is no repeated unit mesh to instance, because
// every ribbon vertex is different. Coordinates are already world-space (Angstrom), so the
// vertex shader only applies the view-projection.
struct IndexedMesh {
    std::vector<MeshVertexRgba> vertices;
    std::vector<std::uint32_t>  indices;

    [[nodiscard]] bool empty() const { return indices.empty() || vertices.empty(); }
};

// A whole molecule turned into draw instances plus its bounding sphere (used by
// the camera to auto-fit). `center`/`radius` are in the same Angstrom space.
struct MolScene {
    std::vector<AtomInst> atoms;
    std::vector<BondInst> bonds;
    IndexedMesh           mesh;   // optional ribbon / surface geometry, drawn after the sticks
    chem::Vec3 center;     // geometric center of the heavy+H atoms
    float       radius = 1.0f;  // bounding radius (>= 1, never zero)
};

// Recompute center/radius over everything the scene contains, including the mesh. buildMolScene
// does this for a conformer; a caller that assembles a scene itself (the protein cartoon) needs
// it too, and duplicating the arithmetic in the panel is how the two drift apart.
void computeSceneBounds(MolScene& scene);

// CPK element color as 0xAABBGGRR (alpha = 0xFF). Falls back to pink for any
// element outside the explicit table (matches PyMOL/Jmol conventions).
std::uint32_t cpkColor(int z);

// Build the draw scene from a conformer.
//  - spacefill = true  -> atoms use van-der-Waals radii and NO bonds are emitted.
//  - spacefill = false -> ball-and-stick: atoms use a fraction of the covalent
//                         radius and bonds become cylinders.
//  - showH = false     -> hydrogens (and any bond touching one) are omitted.
// Returns an empty scene (no atoms/bonds, unit bounding sphere) for an empty
// conformer so callers can treat "nothing to draw" uniformly.
MolScene buildMolScene(const chem::Conformer& conf, bool spacefill, bool showH);

// ----------------------------------------------------------------------------
// GPU renderer.
// ----------------------------------------------------------------------------

class MolViewport {
public:
    MolViewport(ID3D11Device* device, ID3D11DeviceContext* context);
    ~MolViewport();
    MolViewport(const MolViewport&) = delete;
    MolViewport& operator=(const MolViewport&) = delete;

    [[nodiscard]] bool valid() const { return ready_; }

    // Frame the given scene so it fits the view (resets orbit to a default pose).
    // Call when the displayed molecule changes; draw() will otherwise preserve
    // the user's orbit/zoom across frames.
    void resetCamera(const MolScene& scene);

    // Render `scene` at the requested pixel size into the off-screen target and
    // return its shader-resource-view as an ImGui texture id (cast to size_t ->
    // ImTextureID). Returns 0 on failure or zero size. The RT is lazily (re)created
    // when the size changes.
    [[nodiscard]] void* draw(const MolScene& scene, int width, int height);

    // Camera controls (call from the panel when the ImGui image is hovered).
    void orbit(float dxPixels, float dyPixels);  // left-drag
    void zoom(float wheelDelta);                  // mouse wheel (1 notch = +-1)
    void pan(float dxPixels, float dyPixels);     // middle/right-drag

private:
    bool createDeviceResources();
    bool ensureTarget(int width, int height);
    void releaseTarget();
    void releaseAll();

    ID3D11Device*        dev_ = nullptr;
    ID3D11DeviceContext* ctx_ = nullptr;
    bool ready_ = false;

    // Off-screen target (color + depth) + its SRV for ImGui.
    ID3D11Texture2D*           colorTex_ = nullptr;
    ID3D11RenderTargetView*    colorRtv_ = nullptr;
    ID3D11ShaderResourceView*  colorSrv_ = nullptr;
    ID3D11Texture2D*           depthTex_ = nullptr;
    ID3D11DepthStencilView*    depthDsv_ = nullptr;
    int rtW_ = 0, rtH_ = 0;

    // Mesh geometry (unit sphere + unit cylinder), shared by all instances.
    ID3D11Buffer* sphereVB_ = nullptr;
    ID3D11Buffer* sphereIB_ = nullptr;
    unsigned      sphereIndexCount_ = 0;
    ID3D11Buffer* cylVB_ = nullptr;
    ID3D11Buffer* cylIB_ = nullptr;
    unsigned      cylIndexCount_ = 0;

    // Per-instance vertex buffers (resized as needed).
    ID3D11Buffer* atomInstVB_ = nullptr;
    unsigned      atomInstCap_ = 0;
    ID3D11Buffer* bondInstVB_ = nullptr;
    unsigned      bondInstCap_ = 0;

    // Generic indexed-mesh path (the protein cartoon). Both buffers are DYNAMIC and grow by
    // doubling, the same policy as the instance buffers above: a ribbon is rebuilt whenever the
    // structure or the secondary-structure assignment changes, not once at load.
    ID3D11Buffer* meshVB_ = nullptr;
    unsigned      meshVBCap_ = 0;   // in vertices
    ID3D11Buffer* meshIB_ = nullptr;
    unsigned      meshIBCap_ = 0;   // in indices

    // Constant buffer (view-proj + light dir).
    ID3D11Buffer* cbuf_ = nullptr;

    // Shaders + state.
    ID3D11InputLayout*       atomLayout_ = nullptr;
    ID3D11InputLayout*       bondLayout_ = nullptr;
    ID3D11VertexShader*      atomVS_ = nullptr;
    ID3D11PixelShader*       atomPS_ = nullptr;
    ID3D11VertexShader*      bondVS_ = nullptr;
    ID3D11PixelShader*       bondPS_ = nullptr;
    ID3D11RasterizerState*   raster_ = nullptr;
    ID3D11DepthStencilState* depthState_ = nullptr;
    ID3D11InputLayout*       meshLayout_ = nullptr;
    ID3D11VertexShader*      meshVS_ = nullptr;
    ID3D11PixelShader*       meshPS_ = nullptr;
    // The ribbon is drawn with culling OFF and two-sided shading. A cartoon cross-section is a
    // thin closed tube whose triangle winding depends on the sign of the local frame, which the
    // geometry stage does not fix; culling it would drop faces on half the ribbon. A flat ribbon
    // is also legitimately viewed from behind.
    ID3D11RasterizerState*   rasterNoCull_ = nullptr;

    // Orbit camera state (spherical around a target point).
    float yaw_ = 0.6f;       // radians
    float pitch_ = 0.35f;    // radians
    float dist_ = 8.0f;      // camera distance (Angstrom)
    float fitDist_ = 8.0f;   // distance computed by resetCamera (for clip planes)
    chem::Vec3 target_{};    // look-at point
    float sceneRadius_ = 4.0f;

    bool ensureAtomInstCapacity(unsigned count);
    bool ensureBondInstCapacity(unsigned count);
    bool ensureMeshCapacity(unsigned vertexCount, unsigned indexCount);
};

}  // namespace biocad::render
