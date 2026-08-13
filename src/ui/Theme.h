// ui/Theme.h - one clean, modern dark theme + shared accent colors + helpers.
#pragma once

#include <imgui.h>

namespace biocad::theme {

// ---- Palette ---------------------------------------------------------------
// Background / surface tokens
inline constexpr ImU32 kBgSunken      = IM_COL32(9,   11,  16,  255);
inline constexpr ImU32 kBg            = IM_COL32(13,  16,  23,  255);
inline constexpr ImU32 kSurface       = IM_COL32(20,  24,  33,  255);
inline constexpr ImU32 kSurfaceHi     = IM_COL32(28,  34,  45,  255);
inline constexpr ImU32 kSurfaceActive = IM_COL32(36,  44,  58,  255);  // NEW
inline constexpr ImU32 kBorder        = IM_COL32(38,  45,  58,  255);
inline constexpr ImU32 kBorderStrong  = IM_COL32(52,  62,  79,  255);  // NEW

// Text tokens
inline constexpr ImU32 kTextHi        = IM_COL32(233, 237, 244, 255);
inline constexpr ImU32 kText          = IM_COL32(197, 205, 217, 255);  // NEW
inline constexpr ImU32 kTextDim       = IM_COL32(138, 149, 167, 255);
inline constexpr ImU32 kTextFaint     = IM_COL32(95,  105, 122, 255);  // NEW

// Primary / accent tokens
inline constexpr ImU32 kPrimary       = IM_COL32(110, 98,  246, 255);  // NEW, indigo
inline constexpr ImU32 kPrimaryBright = IM_COL32(139, 130, 255, 255);  // NEW
inline constexpr ImU32 kPrimarySoft   = IM_COL32(110, 98,  246, 38);   // NEW
inline constexpr ImU32 kPrimaryFaint  = IM_COL32(110, 98,  246, 20);   // NEW
inline constexpr ImU32 kAccent2       = IM_COL32(192, 132, 252, 255);  // NEW, violet secondary

// Back-compat tokens — repointed from teal to indigo
inline constexpr ImU32 kAccent        = IM_COL32(110, 98,  246, 255);  // was teal, now indigo
inline constexpr ImU32 kAccentSoft    = IM_COL32(110, 98,  246, 70);   // was teal

// Semantic / verdict tokens
inline constexpr ImU32 kGood          = IM_COL32(52,  211, 153, 255);
inline constexpr ImU32 kInfo          = IM_COL32(96,  165, 250, 255);
inline constexpr ImU32 kWarn          = IM_COL32(251, 191, 36,  255);
inline constexpr ImU32 kDanger        = IM_COL32(248, 113, 113, 255);
inline constexpr ImU32 kHighlight     = IM_COL32(245, 159, 64,  255);

// ---- Typography ------------------------------------------------------------
struct Fonts { ImFont* body = nullptr; ImFont* semi = nullptr; };
Fonts& fonts();          // returns static instance; WinMain assigns into it after loading

void pushTitle();        // 28 px semibold (or body fallback)
void pushH2();           // 21 px semibold
void pushValue();        // 25 px semibold
void pushSmallStrong();  // 13 px semibold
void popFont();          // ImGui::PopFont()

// ---- Theme application -----------------------------------------------------
void apply();

ImVec4 verdictColor(int verdict);  // 0=Info 1=Good 2=Warn 3=Danger

// ---- Shared widgets --------------------------------------------------------

// Accent-tinted small-caps section heading followed by a thin rule.
void sectionHeader(const char* label);

// Rounded "pill" / chip drawn at the current cursor.
void pill(const char* text, ImU32 bg = kAccentSoft, ImU32 fg = kTextHi);

// General rounded pill badge.
void badge(const char* text, ImU32 fg = kTextHi, ImU32 bg = kPrimarySoft);

// Maps verdict (0=Info..3=Danger) to colored badge.
void verdictBadge(int verdict, const char* text);

// Card container: BeginChild that never grows a scrollbar. Pair with endCard().
bool beginCard(const char* id, ImVec2 size, bool border = true);
void endCard();

// Fixed-height metric tile.
void metricCard(const char* title, const char* value, const char* sub,
                ImU32 valueColor = kTextHi, float width = 0.0f, float height = 112.0f);

}  // namespace biocad::theme
