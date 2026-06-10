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

    // Rounded, airy, modern. A touch more padding/spacing improves readability.
    s.WindowRounding    = 10.0f;
    s.ChildRounding     = 9.0f;
    s.FrameRounding     = 7.0f;
    s.PopupRounding     = 8.0f;
    s.GrabRounding      = 7.0f;
    s.TabRounding       = 8.0f;
    s.ScrollbarRounding = 10.0f;
    s.WindowBorderSize  = 1.0f;
    s.ChildBorderSize   = 1.0f;
    s.FrameBorderSize   = 0.0f;
    s.WindowPadding     = ImVec2(16, 14);
    s.FramePadding      = ImVec2(12, 7);
    s.CellPadding       = ImVec2(8, 6);
    s.ItemSpacing       = ImVec2(10, 9);
    s.ItemInnerSpacing  = ImVec2(8, 6);
    s.ScrollbarSize     = 13.0f;
    s.GrabMinSize       = 12.0f;
    s.WindowTitleAlign  = ImVec2(0.02f, 0.5f);
    s.SeparatorTextBorderSize = 2.0f;
    s.SeparatorTextPadding    = ImVec2(18, 4);

    ImVec4* c = s.Colors;
    const ImVec4 bg0 = rgba(15, 19, 26);     // window bg (a hair deeper for contrast)
    const ImVec4 bg1 = rgba(23, 28, 37);     // child/frame
    const ImVec4 bg2 = rgba(33, 41, 54);     // hovered
    const ImVec4 accent = rgba(45, 212, 191);
    const ImVec4 accentDim = rgba(45, 212, 191, 150);
    const ImVec4 accentSoft = rgba(45, 212, 191, 46);
    const ImVec4 text = rgba(228, 234, 242);
    const ImVec4 textDim = rgba(143, 158, 178);

    c[ImGuiCol_Text]                = text;
    c[ImGuiCol_TextDisabled]        = textDim;
    c[ImGuiCol_WindowBg]            = bg0;
    c[ImGuiCol_ChildBg]             = rgba(20, 25, 33);
    c[ImGuiCol_PopupBg]             = rgba(20, 25, 33);
    c[ImGuiCol_Border]              = rgba(46, 55, 70);
    c[ImGuiCol_BorderShadow]        = rgba(0, 0, 0, 0);
    c[ImGuiCol_FrameBg]             = bg1;
    c[ImGuiCol_FrameBgHovered]      = bg2;
    c[ImGuiCol_FrameBgActive]       = rgba(40, 50, 65);
    c[ImGuiCol_TitleBg]             = rgba(13, 17, 23);
    c[ImGuiCol_TitleBgActive]       = rgba(17, 22, 30);
    c[ImGuiCol_TitleBgCollapsed]    = rgba(13, 17, 23);
    c[ImGuiCol_MenuBarBg]           = rgba(13, 17, 23);
    c[ImGuiCol_ScrollbarBg]         = rgba(13, 17, 23, 0);
    c[ImGuiCol_ScrollbarGrab]       = rgba(46, 55, 70);
    c[ImGuiCol_ScrollbarGrabHovered]= rgba(62, 74, 92);
    c[ImGuiCol_ScrollbarGrabActive] = accentDim;
    c[ImGuiCol_CheckMark]           = accent;
    c[ImGuiCol_SliderGrab]          = accentDim;
    c[ImGuiCol_SliderGrabActive]    = accent;
    c[ImGuiCol_Button]              = rgba(35, 44, 58);
    c[ImGuiCol_ButtonHovered]       = rgba(47, 59, 77);
    c[ImGuiCol_ButtonActive]        = accentDim;
    c[ImGuiCol_Header]              = accentSoft;
    c[ImGuiCol_HeaderHovered]       = rgba(47, 59, 77);
    c[ImGuiCol_HeaderActive]        = accentDim;
    c[ImGuiCol_Separator]           = rgba(40, 48, 61);
    c[ImGuiCol_SeparatorHovered]    = accentDim;
    c[ImGuiCol_SeparatorActive]     = accent;
    c[ImGuiCol_Tab]                 = rgba(22, 27, 36);
    c[ImGuiCol_TabHovered]          = accentDim;
    c[ImGuiCol_TabSelected]         = rgba(40, 50, 65);
    c[ImGuiCol_TabSelectedOverline] = accent;
    c[ImGuiCol_TabDimmed]           = rgba(19, 24, 32);
    c[ImGuiCol_TabDimmedSelected]   = rgba(30, 37, 49);
    c[ImGuiCol_PlotLines]           = accent;
    c[ImGuiCol_PlotHistogram]       = accent;
    c[ImGuiCol_TableHeaderBg]       = rgba(26, 32, 42);
    c[ImGuiCol_TableBorderStrong]   = rgba(46, 55, 70);
    c[ImGuiCol_TableBorderLight]    = rgba(34, 41, 53);
    c[ImGuiCol_TableRowBg]          = rgba(0, 0, 0, 0);
    c[ImGuiCol_TableRowBgAlt]       = rgba(255, 255, 255, 6);
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

void sectionHeader(const char* label) {
    ImGui::Dummy(ImVec2(0, 1));
    ImGui::PushStyleColor(ImGuiCol_Text, ImGui::ColorConvertU32ToFloat4(kAccent));
    ImGui::TextUnformatted(label);
    ImGui::PopStyleColor();
    ImGui::PushStyleColor(ImGuiCol_Separator, kBorder);
    ImGui::Separator();
    ImGui::PopStyleColor();
    ImGui::Dummy(ImVec2(0, 2));
}

void pill(const char* text, ImU32 bg, ImU32 fg) {
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const ImVec2 pad(9, 3);
    const ImVec2 sz = ImGui::CalcTextSize(text);
    const ImVec2 p0 = ImGui::GetCursorScreenPos();
    const ImVec2 p1(p0.x + sz.x + pad.x * 2, p0.y + sz.y + pad.y * 2);
    const float rounding = (p1.y - p0.y) * 0.5f;
    dl->AddRectFilled(p0, p1, bg, rounding);
    dl->AddText(ImVec2(p0.x + pad.x, p0.y + pad.y), fg, text);
    ImGui::Dummy(ImVec2(sz.x + pad.x * 2, sz.y + pad.y * 2));
}

}  // namespace stimlab::theme
