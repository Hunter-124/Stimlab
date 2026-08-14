#include "ui/Theme.h"

#include "data/Domain.h"

namespace biocad::theme {

namespace {
ImVec4 rgba(int r, int g, int b, int a = 255) {
    return ImVec4(r / 255.0f, g / 255.0f, b / 255.0f, a / 255.0f);
}

ImVec4 u32ToVec4(ImU32 u) { return ImGui::ColorConvertU32ToFloat4(u); }

// Pick semibold face if loaded, otherwise fall back to body.
ImFont* semiFace() {
    const Fonts& f = fonts();
    return (f.semi != nullptr) ? f.semi : f.body;
}

ImFont* monoFace() {
    const Fonts& f = fonts();
    return (f.mono != nullptr) ? f.mono : f.body;
}

// Small-caps style label: 13px semibold, dim.
void smallCaps(const char* label, ImU32 color) {
    pushSmallStrong();
    ImGui::PushStyleColor(ImGuiCol_Text, u32ToVec4(color));
    ImGui::TextUnformatted(label);
    ImGui::PopStyleColor();
    popFont();
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
void pushMono()        { ImGui::PushFont(monoFace(), 15.0f); }
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

    const ImVec4 bgSunken      = u32ToVec4(kBgSunken);
    const ImVec4 bg            = u32ToVec4(kBg);
    const ImVec4 surface       = u32ToVec4(kSurface);
    const ImVec4 surfaceHi     = u32ToVec4(kSurfaceHi);
    const ImVec4 surfaceActive = u32ToVec4(kSurfaceActive);
    const ImVec4 border        = u32ToVec4(kBorder);
    const ImVec4 borderStrong  = u32ToVec4(kBorderStrong);
    const ImVec4 text          = u32ToVec4(kText);
    const ImVec4 textDim       = u32ToVec4(kTextDim);
    const ImVec4 primary       = u32ToVec4(kPrimary);
    const ImVec4 primaryBright = u32ToVec4(kPrimaryBright);
    const ImVec4 primarySoft   = u32ToVec4(kPrimarySoft);
    const ImVec4 transparent   = rgba(0, 0, 0, 0);

    c[ImGuiCol_Text]                = text;
    c[ImGuiCol_TextDisabled]        = textDim;
    c[ImGuiCol_WindowBg]            = bg;
    c[ImGuiCol_ChildBg]             = surface;
    c[ImGuiCol_PopupBg]             = surfaceHi;
    c[ImGuiCol_Border]              = border;
    c[ImGuiCol_BorderShadow]        = transparent;
    c[ImGuiCol_FrameBg]             = surface;
    c[ImGuiCol_FrameBgHovered]      = surfaceHi;
    c[ImGuiCol_FrameBgActive]       = surfaceActive;
    c[ImGuiCol_TitleBg]             = bgSunken;
    c[ImGuiCol_TitleBgActive]       = bgSunken;
    c[ImGuiCol_TitleBgCollapsed]    = bgSunken;
    c[ImGuiCol_MenuBarBg]           = bgSunken;
    c[ImGuiCol_ScrollbarBg]         = transparent;
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
    c[ImGuiCol_PlotHistogram]       = primary;
    c[ImGuiCol_TableHeaderBg]       = surfaceActive;
    c[ImGuiCol_TableBorderStrong]   = border;
    c[ImGuiCol_TableBorderLight]    = surfaceHi;
    c[ImGuiCol_TableRowBg]          = transparent;
    c[ImGuiCol_TableRowBgAlt]       = rgba(255, 255, 255, 5);
    c[ImGuiCol_DockingPreview]      = primary;
    c[ImGuiCol_TextSelectedBg]      = primarySoft;
    c[ImGuiCol_NavHighlight]        = primary;
}

ImVec4 verdictColor(int verdict) {
    switch (static_cast<Verdict>(verdict)) {
        case Verdict::Good:   return u32ToVec4(kGood);
        case Verdict::Warn:   return u32ToVec4(kWarn);
        case Verdict::Danger: return u32ToVec4(kDanger);
        case Verdict::Info:
        default:              return u32ToVec4(kInfo);
    }
}

// Provenance colours: green measured / blue predicted / purple model /
// amber heuristic / grey not-computed. Deliberately distinct from the verdict
// scale: provenance says how a number was obtained, not whether it is good news.
ImVec4 provenanceColor(Provenance p) {
    switch (p) {
        case Provenance::Measured:    return u32ToVec4(kGood);
        case Provenance::Predicted:   return u32ToVec4(kInfo);
        case Provenance::Model:       return u32ToVec4(kAccent2);
        case Provenance::Heuristic:   return u32ToVec4(kWarn);
        case Provenance::NotComputed: return u32ToVec4(kTextDim);
    }
    return u32ToVec4(kTextDim);
}

// ---- Widgets ---------------------------------------------------------------

void sectionHeader(const char* label) {
    ImGui::Dummy(ImVec2(0, 2));
    // 3px accent tick ahead of the small-caps label: sections scan as sections.
    const ImVec2 p0 = ImGui::GetCursorScreenPos();
    const float h = ImGui::GetTextLineHeight();
    ImGui::GetWindowDrawList()->AddRectFilled(
        p0, ImVec2(p0.x + 3.0f, p0.y + h), kPrimary, 1.5f);
    ImGui::SetCursorPosX(ImGui::GetCursorPosX() + 10.0f);
    smallCaps(label, kTextDim);
    ImGui::PushStyleColor(ImGuiCol_Separator, u32ToVec4(kBorder));
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
    const ImU32 bg = IM_COL32(
        static_cast<int>(col.x * 255.0f),
        static_cast<int>(col.y * 255.0f),
        static_cast<int>(col.z * 255.0f),
        38);
    badge(text, fg, bg);
}

bool beginCard(const char* id, ImVec2 size, bool border) {
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(16, 13));
    ImGuiChildFlags cf = border ? ImGuiChildFlags_Borders : ImGuiChildFlags_None;
    // size.y == 0 means "as tall as the content": without AutoResizeY a zero-height
    // child eats ALL remaining window space and everything drawn after the card
    // lands off-screen.
    if (size.y == 0.0f) cf |= ImGuiChildFlags_AutoResizeY;
    bool visible = ImGui::BeginChild(id, size, cf,
        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
    ImGui::PopStyleVar();  // pop immediately — padding is captured at Begin
    return visible;
}

void endCard() {
    ImGui::EndChild();
}

bool beginTitledCard(const char* id, const char* title, ImVec2 size,
                     const char* rightNote) {
    const bool visible = beginCard(id, size);
    if (!visible) return false;
    if (rightNote != nullptr) {
        const float noteW = ImGui::CalcTextSize(rightNote).x;
        const float x = ImGui::GetCursorPosX() + ImGui::GetContentRegionAvail().x - noteW;
        smallCaps(title, kTextDim);
        ImGui::SameLine();
        ImGui::SetCursorPosX(x);
        ImGui::PushStyleColor(ImGuiCol_Text, u32ToVec4(kTextFaint));
        ImGui::TextUnformatted(rightNote);
        ImGui::PopStyleColor();
    } else {
        smallCaps(title, kTextDim);
    }
    ImGui::Dummy(ImVec2(0, 4));
    return true;
}

void metricCard(const char* title, const char* value, const char* sub,
                ImU32 valueColor, float width, float height) {
    const bool visible = beginCard(title, ImVec2(width, height));
    if (!visible) {
        endCard();
        return;
    }

    smallCaps(title, kTextDim);

    ImGui::Spacing();

    pushValue();
    ImGui::PushStyleColor(ImGuiCol_Text, u32ToVec4(valueColor));
    ImGui::TextUnformatted(value);
    ImGui::PopStyleColor();
    popFont();

    ImGui::Spacing();

    smallCaps(sub, kTextDim);

    endCard();
}

void kvRow(const char* label, const char* value, float labelWidth) {
    ImGui::PushStyleColor(ImGuiCol_Text, u32ToVec4(kTextDim));
    ImGui::TextUnformatted(label);
    ImGui::PopStyleColor();
    ImGui::SameLine(labelWidth);
    ImGui::PushStyleColor(ImGuiCol_Text, u32ToVec4(kTextHi));
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextUnformatted(value);
    ImGui::PopTextWrapPos();
    ImGui::PopStyleColor();
}

void statusDot(ImU32 color, const char* tooltip) {
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const ImVec2 p = ImGui::GetCursorScreenPos();
    const float r = 4.0f;
    const ImVec2 c(p.x + r + 1.0f, p.y + ImGui::GetTextLineHeight() * 0.5f);
    dl->AddCircleFilled(c, r, color);
    ImGui::Dummy(ImVec2(r * 2.0f + 2.0f, ImGui::GetTextLineHeight()));
    if (tooltip != nullptr && ImGui::IsItemHovered())
        ImGui::SetTooltip("%s", tooltip);
}

void bubble(const char* text, ImU32 bg, ImU32 fg, bool alignRight,
            float maxWidthFraction) {
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const float avail = ImGui::GetContentRegionAvail().x;
    const float wrapW = avail * maxWidthFraction;
    const ImVec2 pad(11, 8);
    // Measure the wrapped text first so the rect fits it exactly.
    const ImVec2 sz = ImGui::CalcTextSize(text, nullptr, false, wrapW - pad.x * 2.0f);
    const ImVec2 box(sz.x + pad.x * 2.0f, sz.y + pad.y * 2.0f);

    float x = ImGui::GetCursorScreenPos().x;
    if (alignRight) x += avail - box.x;
    const ImVec2 p0(x, ImGui::GetCursorScreenPos().y);
    const ImVec2 p1(p0.x + box.x, p0.y + box.y);
    dl->AddRectFilled(p0, p1, bg, 8.0f);
    dl->AddText(nullptr, 0.0f, ImVec2(p0.x + pad.x, p0.y + pad.y), fg, text,
                nullptr, wrapW - pad.x * 2.0f);
    ImGui::SetCursorScreenPos(ImVec2(ImGui::GetCursorScreenPos().x, p1.y + 4.0f));
    ImGui::Dummy(ImVec2(box.x, 0));
}

}  // namespace biocad::theme
