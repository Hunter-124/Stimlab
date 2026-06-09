// ui/AppShell.h - top-level UI: host dockspace, left navigator, central content
// router, right assistant. Holds the shared UiState and the assistant->UI
// highlight mechanism (the agent can pulse any panel to answer "how do I...").
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "contracts/Services.h"
#include "data/Domain.h"

struct ID3D11Device;
struct ID3D11DeviceContext;

namespace stimlab {

namespace render { class MolViewport; }


struct PanelInfo {
    std::string id;     // stable key used by the content router + highlight
    std::string label;  // display label
    std::string help;   // one-line description shown by the assistant
};

struct UiState {
    std::string selectedMolecule = "amphetamine";
    std::string activePanel = "Dashboard";
    std::string dockTarget;                 // chosen docking target
    std::string highlight;                  // panel id currently pulsing ("" = none)
    double      highlightStart = -1000.0;   // ImGui time when highlight began
    bool        showAssistant = true;
    bool        showAbout = false;
    bool        quitRequested = false;
    std::vector<std::string> assistantLog;
};

class AppShell {
public:
    explicit AppShell(Services services);
    ~AppShell();  // defined in the .cpp (MolViewport is incomplete here)

    void draw();  // call once per frame between NewFrame() and Render()

    UiState&  state()    { return state_; }
    Services& services() { return svc_; }
    [[nodiscard]] Molecule currentMolecule() const;

    // 3D viewer wiring (WP-2). WinMain hands us the live DX11 device so panels can
    // render the off-screen molecular viewport via viewer().
    void setRenderDevice(ID3D11Device* device, ID3D11DeviceContext* context);
    // Lazily constructs the viewport when a device has been set; returns nullptr
    // if no device is available (panels must guard).
    [[nodiscard]] render::MolViewport* viewer();

    // Assistant -> UI bridge.
    void requestHighlight(const std::string& panelId, const std::string& explanation);
    void frameHighlightCurrentWindow(const std::string& panelId);  // call inside a panel
    [[nodiscard]] bool isHighlighted(const std::string& panelId) const;

    [[nodiscard]] const std::vector<PanelInfo>& panels() const { return panels_; }

private:
    void drawMainMenuBar();
    void drawNavigator();
    void drawContent();
    void drawAssistant();
    void drawAboutModal();

    Services svc_;
    UiState  state_;
    std::vector<PanelInfo> panels_;

    ID3D11Device*        renderDev_ = nullptr;
    ID3D11DeviceContext* renderCtx_ = nullptr;
    std::unique_ptr<render::MolViewport> viewer_;
};

}  // namespace stimlab
