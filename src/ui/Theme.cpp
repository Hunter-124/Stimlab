#include "ui/Theme.h"

#include "data/Domain.h"

namespace stimlab::theme {

namespace {
ImVec4 rgba(int r, int g, int b, int a = 255) {
    return ImVec4(r / 255.0f, g / 255.0f, b / 255.0f, a / 255.0f);
}
}  // namespace

void apply() {
    ImGuiStyle& s = ImGui::GetStyle();

    s.WindowRounding   = 7.0f;
    s.ChildRounding    = 7.0f;
    s.FrameRounding    = 6.0f;
    s.PopupRounding    = 6.0f;
    s.GrabRounding     = 6.0f;
    s.TabRounding      = 6.0f;
    s.ScrollbarRounding = 9.0f;
    s.WindowBorderSize = 1.0f;
    s.FrameBorderSize  = 0.0f;
    s.WindowPadding    = ImVec2(14, 12);
    s.FramePadding     = ImVec2(11, 6);
    s.ItemSpacing      = ImVec2(10, 8);
    s.ItemInnerSpacing = ImVec2(8, 6);
    s.ScrollbarSize    = 13.0f;
    s.GrabMinSize      = 11.0f;
    s.WindowTitleAlign = ImVec2(0.02f, 0.5f);

    ImVec4* c = s.Colors;
    const ImVec4 bg0 = rgba(17, 21, 28);     // window bg
    const ImVec4 bg1 = rgba(23, 28, 37);     // child/frame
    const ImVec4 bg2 = rgba(31, 38, 50);     // hovered
    const ImVec4 accent = rgba(45, 212, 191);
    const ImVec4 accentDim = rgba(45, 212, 191, 140);
    const ImVec4 text = rgba(226, 232, 240);
    const ImVec4 textDim = rgba(148, 163, 184);

    c[ImGuiCol_Text]                = text;
    c[ImGuiCol_TextDisabled]        = textDim;
    c[ImGuiCol_WindowBg]            = bg0;
    c[ImGuiCol_ChildBg]             = rgba(20, 25, 33);
    c[ImGuiCol_PopupBg]             = rgba(20, 25, 33);
    c[ImGuiCol_Border]              = rgba(44, 52, 66);
    c[ImGuiCol_BorderShadow]        = rgba(0, 0, 0, 0);
    c[ImGuiCol_FrameBg]             = bg1;
    c[ImGuiCol_FrameBgHovered]      = bg2;
    c[ImGuiCol_FrameBgActive]       = rgba(38, 47, 61);
    c[ImGuiCol_TitleBg]             = rgba(14, 18, 24);
    c[ImGuiCol_TitleBgActive]       = rgba(18, 23, 31);
    c[ImGuiCol_TitleBgCollapsed]    = rgba(14, 18, 24);
    c[ImGuiCol_MenuBarBg]           = rgba(14, 18, 24);
    c[ImGuiCol_ScrollbarBg]         = rgba(14, 18, 24);
    c[ImGuiCol_ScrollbarGrab]       = rgba(44, 52, 66);
    c[ImGuiCol_ScrollbarGrabHovered]= rgba(60, 70, 88);
    c[ImGuiCol_ScrollbarGrabActive] = accentDim;
    c[ImGuiCol_CheckMark]           = accent;
    c[ImGuiCol_SliderGrab]          = accentDim;
    c[ImGuiCol_SliderGrabActive]    = accent;
    c[ImGuiCol_Button]              = rgba(34, 42, 55);
    c[ImGuiCol_ButtonHovered]       = rgba(45, 56, 73);
    c[ImGuiCol_ButtonActive]        = accentDim;
    c[ImGuiCol_Header]              = rgba(34, 42, 55);
    c[ImGuiCol_HeaderHovered]       = rgba(45, 56, 73);
    c[ImGuiCol_HeaderActive]        = accentDim;
    c[ImGuiCol_Separator]           = rgba(44, 52, 66);
    c[ImGuiCol_SeparatorHovered]    = accentDim;
    c[ImGuiCol_SeparatorActive]     = accent;
    c[ImGuiCol_Tab]                 = rgba(23, 28, 37);
    c[ImGuiCol_TabHovered]          = accentDim;
    c[ImGuiCol_TabActive]           = rgba(38, 47, 61);
    c[ImGuiCol_TabUnfocused]        = rgba(20, 25, 33);
    c[ImGuiCol_TabUnfocusedActive]  = rgba(31, 38, 50);
    c[ImGuiCol_PlotLines]           = accent;
    c[ImGuiCol_PlotHistogram]       = accent;
    c[ImGuiCol_DockingPreview]      = accentDim;
    c[ImGuiCol_TextSelectedBg]      = accentDim;
    c[ImGuiCol_NavHighlight]        = accent;
}

ImVec4 verdictColor(int verdict) {
    switch (static_cast<Verdict>(verdict)) {
        case Verdict::Good:   return rgba(74, 222, 128);
        case Verdict::Warn:   return rgba(250, 204, 21);
        case Verdict::Danger: return rgba(248, 113, 113);
        case Verdict::Info:
        default:              return rgba(96, 165, 250);
    }
}

}  // namespace stimlab::theme
