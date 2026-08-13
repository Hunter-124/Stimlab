#include "ui/Theme.h"

#include "data/Domain.h"

namespace biocad::theme {

namespace {
ImVec4 rgba(int r, int g, int b, int a = 255) {
    return ImVec4(r / 255.0f, g / 255.0f, b / 255.0f, a / 255.0f);
}

// Pick semibold face if loaded, otherwise fall back to body.
ImFont* semiFace() {
    const Fonts& f = fonts();
    return (f.semi != nullptr) ? f.semi : f.body;
}
}  // namespace

// ---- Typography ------------------------------------------------------------

Fonts& fonts() {
    static Fonts f;
    return f;
}

void pushTitle()       { ImGui::PushFont(semiFace(), 28.0f); }
void pushH2()          { ImGui::PushFont(semiFace(), 21.0f); }
void pushValue()       { ImGui::PushFont(semiFace(), 25.0f); }
void pushSmallStrong() { ImGui::PushFont(semiFace(), 13.0f); }
void popFont()         { ImGui::PopFont(); }

// ---- Theme application -----------------------------------------------------

void apply() {
    ImGuiStyle& s = ImGui::GetStyle();

    // Rounded, airy, modern.
    s.WindowRounding    = 10.0f;
    s.ChildRounding     = 10.0f;
    s.FrameRounding     = 8.0f;
    s.PopupRounding     = 8.0f;
    s.GrabRounding      = 7.0f;
    s.TabRounding       = 8.0f;
    s.ScrollbarRounding = 10.0f;
    s.WindowBorderSize  = 1.0f;
    s.ChildBorderSize   = 1.0f;
    s.FrameBorderSize   = 0.0f;
    s.WindowPadding     = ImVec2(16, 12);
    s.FramePadding      = ImVec2(11,  7);
    s.CellPadding       = ImVec2(8,   6);
    s.ItemSpacing       = ImVec2(10,  8);
    s.ItemInnerSpacing  = ImVec2(8,   6);
    s.ScrollbarSize     = 12.0f;
    s.GrabMinSize       = 12.0f;
    s.WindowTitleAlign  = ImVec2(0.02f, 0.5f);
    s.SeparatorTextBorderSize = 2.0f;
    s.SeparatorTextPadding    = ImVec2(18, 4);

    ImVec4* c = s.Colors;

    // Local helpers
    const ImVec4 bgSunken      = rgba(9,   11,  16);
    const ImVec4 bg            = rgba(13,  16,  23);
    const ImVec4 surface       = rgba(20,  24,  33);
    const ImVec4 surfaceHi     = rgba(28,  34,  45);
    const ImVec4 surfaceActive = rgba(36,  44,  58);
    const ImVec4 border        = rgba(38,  45,  58);
    const ImVec4 borderStrong  = rgba(52,  62,  79);
    const ImVec4 text          = rgba(197, 205, 217);
    const ImVec4 textDim       = rgba(138, 149, 167);
    const ImVec4 primary       = rgba(110, 98,  246);
    const ImVec4 primaryBright = rgba(139, 130, 255);
    const ImVec4 primarySoft   = rgba(110, 98,  246, 38);
    const ImVec4 transparent   = rgba(0,   0,   0,   0);

    c[ImGuiCol_Text]                = text;
    c[ImGuiCol_TextDisabled]        = textDim;
    c[ImGuiCol_WindowBg]            = bg;
    c[ImGuiCol_ChildBg]             = surface;
    c[ImGuiCol_PopupBg]             = surface;
    c[ImGuiCol_Border]              = border;
    c[ImGuiCol_BorderShadow]        = transparent;
    c[ImGuiCol_FrameBg]             = surface;
    c[ImGuiCol_FrameBgHovered]      = surfaceHi;
    c[ImGuiCol_FrameBgActive]       = surfaceActive;
    c[ImGuiCol_TitleBg]             = bgSunken;
    c[ImGuiCol_TitleBgActive]       = bgSunken;
    c[ImGuiCol_TitleBgCollapsed]    = bgSunken;
    c[ImGuiCol_MenuBarBg]           = bgSunken;
    c[ImGuiCol_ScrollbarBg]         = rgba(9, 11, 16, 0);   // transparent over kBgSunken
    c[ImGuiCol_ScrollbarGrab]       = borderStrong;
    c[ImGuiCol_ScrollbarGrabHovered]= primary;
    c[ImGuiCol_ScrollbarGrabActive] = primaryBright;
    c[ImGuiCol_CheckMark]           = primaryBright;
    c[ImGuiCol_SliderGrab]          = primary;
    c[ImGuiCol_SliderGrabActive]    = primaryBright;
    c[ImGuiCol_Button]              = surfaceHi;
    c[ImGuiCol_ButtonHovered]       = surfaceActive;
    c[ImGuiCol_ButtonActive]        = primary;
    c[ImGuiCol_Header]              = primarySoft;
    c[ImGuiCol_HeaderHovered]       = surfaceHi;
    c[ImGuiCol_HeaderActive]        = primary;
    c[ImGuiCol_Separator]           = border;
    c[ImGuiCol_SeparatorHovered]    = primary;
    c[ImGuiCol_SeparatorActive]     = primaryBright;
    c[ImGuiCol_Tab]                 = surface;
    c[ImGuiCol_TabHovered]          = surfaceHi;
    c[ImGuiCol_TabSelected]         = surfaceActive;
    c[ImGuiCol_TabSelectedOverline] = primaryBright;
    c[ImGuiCol_TabDimmed]           = surface;
    c[ImGuiCol_TabDimmedSelected]   = surfaceHi;
    c[ImGuiCol_PlotLines]           = primaryBright;
    c[ImGuiCol_PlotHistogram]       = primaryBright;
    c[ImGuiCol_TableHeaderBg]       = surfaceHi;
    c[ImGuiCol_TableBorderStrong]   = border;
    c[ImGuiCol_TableBorderLight]    = surface;
    c[ImGuiCol_TableRowBg]          = transparent;
    c[ImGuiCol_TableRowBgAlt]       = rgba(255, 255, 255, 6);
    c[ImGuiCol_DockingPreview]      = primary;
    c[ImGuiCol_TextSelectedBg]      = primarySoft;
    c[ImGuiCol_NavHighlight]        = primary;
}

ImVec4 verdictColor(int verdict) {
    switch (static_cast<Verdict>(verdict)) {
        case Verdict::Good:   return rgba(52,  211, 153);
        case Verdict::Warn:   return rgba(251, 191, 36);
        case Verdict::Danger: return rgba(248, 113, 113);
        case Verdict::Info:
        default:              return rgba(96,  165, 250);
    }
}

// ---- Widgets ---------------------------------------------------------------

void sectionHeader(const char* label) {
    ImGui::Dummy(ImVec2(0, 2));
    pushSmallStrong();
    ImGui::PushStyleColor(ImGuiCol_Text, ImGui::ColorConvertU32ToFloat4(kTextDim));
    ImGui::TextUnformatted(label);
    ImGui::PopStyleColor();
    popFont();
    ImGui::PushStyleColor(ImGuiCol_Separator, ImGui::ColorConvertU32ToFloat4(kBorder));
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

void badge(const char* text, ImU32 fg, ImU32 bg) {
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

void verdictBadge(int verdict, const char* text) {
    const ImVec4 col = verdictColor(verdict);
    const ImU32 fg = ImGui::ColorConvertFloat4ToU32(col);
    // Build a faint background at ~20% alpha using the same RGB
    const ImU32 bg = IM_COL32(
        static_cast<int>(col.x * 255.0f),
        static_cast<int>(col.y * 255.0f),
        static_cast<int>(col.z * 255.0f),
        38);
    badge(text, fg, bg);
}

bool beginCard(const char* id, ImVec2 size, bool border) {
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(16, 13));
    bool visible = ImGui::BeginChild(id, size,
        border ? ImGuiChildFlags_Borders : ImGuiChildFlags_None,
        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
    ImGui::PopStyleVar();  // pop immediately — padding is captured at Begin
    return visible;
}

void endCard() {
    ImGui::EndChild();
}

void metricCard(const char* title, const char* value, const char* sub,
                ImU32 valueColor, float width, float height) {
    beginCard(title, ImVec2(width, height));

    // Title in small-caps dim style
    pushSmallStrong();
    ImGui::PushStyleColor(ImGuiCol_Text, ImGui::ColorConvertU32ToFloat4(kTextDim));
    ImGui::TextUnformatted(title);
    ImGui::PopStyleColor();
    popFont();

    ImGui::Spacing();

    // Big value in valueColor
    pushValue();
    ImGui::PushStyleColor(ImGuiCol_Text, ImGui::ColorConvertU32ToFloat4(valueColor));
    ImGui::TextUnformatted(value);
    ImGui::PopStyleColor();
    popFont();

    ImGui::Spacing();

    // Sub text in dim
    pushSmallStrong();
    ImGui::PushStyleColor(ImGuiCol_Text, ImGui::ColorConvertU32ToFloat4(kTextDim));
    ImGui::TextUnformatted(sub);
    ImGui::PopStyleColor();
    popFont();

    endCard();
}

}  // namespace biocad::theme
