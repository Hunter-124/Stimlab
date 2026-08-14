// ui/Theme.h - one clean, modern dark theme + shared accent colors + helpers.
#pragma once

#include <imgui.h>

#include "data/Domain.h"

namespace biocad::theme {

// ---- Palette ---------------------------------------------------------------
// Background / surface tokens. Three distinct elevation levels so a card reads
// as a card against the window, and a hovered control reads as interactive.
inline constexpr ImU32 kBgSunken      = IM_COL32(7,   9,   13,  255);
inline constexpr ImU32 kBg            = IM_COL32(11,  14,  20,  255);
inline constexpr ImU32 kSurface       = IM_COL32(17,  21,  30,  255);
inline constexpr ImU32 kSurfaceHi     = IM_COL32(25,  30,  42,  255);
inline constexpr ImU32 kSurfaceActive = IM_COL32(34,  41,  56,  255);
inline constexpr ImU32 kBorder        = IM_COL32(34,  41,  56,  255);
inline constexpr ImU32 kBorderStrong  = IM_COL32(52,  62,  83,  255);

// Text tokens
inline constexpr ImU32 kTextHi        = IM_COL32(235, 239, 247, 255);
inline constexpr ImU32 kText          = IM_COL32(198, 206, 219, 255);
inline constexpr ImU32 kTextDim       = IM_COL32(141, 152, 171, 255);
inline constexpr ImU32 kTextFaint     = IM_COL32(94,  104, 122, 255);

// Primary / accent tokens
inline constexpr ImU32 kPrimary       = IM_COL32(124, 111, 250, 255);  // indigo
inline constexpr ImU32 kPrimaryBright = IM_COL32(150, 141, 255, 255);
inline constexpr ImU32 kPrimarySoft   = IM_COL32(124, 111, 250, 44);
inline constexpr ImU32 kPrimaryFaint  = IM_COL32(124, 111, 250, 22);
inline constexpr ImU32 kOnPrimary     = IM_COL32(246, 246, 255, 255);
inline constexpr ImU32 kAccent2       = IM_COL32(192, 132, 252, 255);  // violet

// Back-compat tokens
inline constexpr ImU32 kAccent        = kPrimary;
inline constexpr ImU32 kAccentSoft    = IM_COL32(124, 111, 250, 78);

// Semantic / verdict tokens
inline constexpr ImU32 kGood          = IM_COL32(52,  211, 153, 255);
inline constexpr ImU32 kInfo          = IM_COL32(96,  165, 250, 255);
inline constexpr ImU32 kWarn          = IM_COL32(251, 191, 36,  255);
inline constexpr ImU32 kDanger        = IM_COL32(248, 113, 113, 255);
inline constexpr ImU32 kHighlight     = IM_COL32(245, 159, 64,  255);

// ---- Icons (Nerd Font private-use glyphs merged into the UI fonts) ---------
// Kept to the Font-Awesome-era range (E000-F2E8) that every Nerd Font ships.
namespace icon {
inline constexpr const char* kHome      = "\xEF\x80\x95";  // F015
inline constexpr const char* kSearch    = "\xEF\x80\x82";  // F002
inline constexpr const char* kFlask     = "\xEF\x83\x83";  // F0C3
inline constexpr const char* kCube      = "\xEF\x86\xB2";  // F1B2
inline constexpr const char* kCubes     = "\xEF\x86\xB3";  // F1B3
inline constexpr const char* kBook      = "\xEF\x80\xAD";  // F02D
inline constexpr const char* kPencil    = "\xEF\x81\x80";  // F040
inline constexpr const char* kChevronR  = "\xEF\x81\x94";  // F054
inline constexpr const char* kChevronD  = "\xEF\x81\xB8";  // F078
inline constexpr const char* kClose     = "\xEF\x80\x8D";  // F00D
inline constexpr const char* kSend      = "\xEF\x87\x98";  // F1D8
inline constexpr const char* kPlay      = "\xEF\x81\x8B";  // F04B
inline constexpr const char* kCog       = "\xEF\x80\x93";  // F013
inline constexpr const char* kCogs      = "\xEF\x82\x85";  // F085
inline constexpr const char* kHistory   = "\xEF\x87\x9A";  // F1DA
inline constexpr const char* kDatabase  = "\xEF\x87\x80";  // F1C0
inline constexpr const char* kSitemap   = "\xEF\x83\xA8";  // F0E8
inline constexpr const char* kAnchor    = "\xEF\x84\xBD";  // F13D
inline constexpr const char* kRoad      = "\xEF\x80\x98";  // F018
inline constexpr const char* kFilter    = "\xEF\x82\xB0";  // F0B0
inline constexpr const char* kCrosshair = "\xEF\x81\x9B";  // F05B
inline constexpr const char* kExchange  = "\xEF\x83\xAC";  // F0EC
inline constexpr const char* kTint      = "\xEF\x81\x83";  // F043
inline constexpr const char* kThermo    = "\xEF\x8B\x89";  // F2C9
inline constexpr const char* kShield    = "\xEF\x84\xB2";  // F132
inline constexpr const char* kWarning   = "\xEF\x81\xB1";  // F071
inline constexpr const char* kUsers     = "\xEF\x83\x80";  // F0C0
inline constexpr const char* kBolt      = "\xEF\x83\xA7";  // F0E7
inline constexpr const char* kCycle     = "\xEF\x80\xA1";  // F021
inline constexpr const char* kColumns   = "\xEF\x83\x9B";  // F0DB
inline constexpr const char* kBranch    = "\xEF\x84\xA6";  // F126
inline constexpr const char* kAsterisk  = "\xEF\x81\xA9";  // F069
inline constexpr const char* kChart     = "\xEF\x88\x81";  // F201
inline constexpr const char* kMagic     = "\xEF\x83\x90";  // F0D0
inline constexpr const char* kBulb      = "\xEF\x83\xAB";  // F0EB
inline constexpr const char* kAlignLeft = "\xEF\x80\xB6";  // F036
inline constexpr const char* kConnect   = "\xEF\x88\x8E";  // F20E
inline constexpr const char* kListOl    = "\xEF\x83\x8B";  // F0CB
inline constexpr const char* kBraille   = "\xEF\x8A\xA1";  // F2A1
inline constexpr const char* kBalance   = "\xEF\x89\x8E";  // F24E
inline constexpr const char* kGavel     = "\xEF\x83\xA3";  // F0E3
inline constexpr const char* kTasks     = "\xEF\x82\xAE";  // F0AE
inline constexpr const char* kCopy      = "\xEF\x83\x85";  // F0C5
inline constexpr const char* kHeart     = "\xEF\x88\x9E";  // F21E
inline constexpr const char* kInfo      = "\xEF\x81\x9A";  // F05A
inline constexpr const char* kRobot     = "\xEF\x95\x8B";  // F54B (falls back to kBolt if absent)
inline constexpr const char* kExport    = "\xEF\x82\x8E";  // F08E
inline constexpr const char* kDot       = "\xEF\x84\x91";  // F111
}  // namespace icon

// ---- Typography ------------------------------------------------------------
struct Fonts { ImFont* body = nullptr; ImFont* semi = nullptr; ImFont* mono = nullptr; };
Fonts& fonts();          // returns static instance; WinMain assigns into it after loading

void pushTitle();        // 28 px semibold (or body fallback)
void pushH2();           // 21 px semibold
void pushValue();        // 25 px semibold
void pushSmallStrong();  // 13 px semibold
void pushMono();         // 15 px monospace (SMILES / sequences / raw data)
void popFont();          // ImGui::PopFont()

// ---- Theme application -----------------------------------------------------
void apply();

ImVec4 verdictColor(int verdict);  // 0=Info 1=Good 2=Warn 3=Danger

// Colour for a number's provenance tier (green/blue/purple/amber/grey).
ImVec4 provenanceColor(Provenance p);

// ---- Shared widgets --------------------------------------------------------

// Accent-tick + small-caps section heading followed by a thin rule.
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

// Titled card: small-caps header + optional dim right-aligned note, then content.
// Pairs with endCard().
bool beginTitledCard(const char* id, const char* title, ImVec2 size,
                     const char* rightNote = nullptr);

// Fixed-height metric tile.
void metricCard(const char* title, const char* value, const char* sub,
                ImU32 valueColor = kTextHi, float width = 0.0f, float height = 112.0f);

// Label/value row: dim label column, high-contrast value (wrapped).
void kvRow(const char* label, const char* value, float labelWidth = 120.0f);

// Filled status dot + optional tooltip on hover.
void statusDot(ImU32 color, const char* tooltip = nullptr);

// Chat bubble: wrapped text inside a rounded filled rect. alignRight puts the
// bubble on the right edge of the content region (user messages).
void bubble(const char* text, ImU32 bg, ImU32 fg, bool alignRight = false,
            float maxWidthFraction = 0.86f);

}  // namespace biocad::theme
