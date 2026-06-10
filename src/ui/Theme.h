// ui/Theme.h - one clean, modern dark theme + shared accent colors + helpers.
#pragma once

#include <imgui.h>

namespace stimlab::theme {

// Shared palette (also used by panels for charts / highlights).
inline constexpr ImU32 kAccent     = IM_COL32(45, 212, 191, 255);   // teal
inline constexpr ImU32 kAccentSoft = IM_COL32(45, 212, 191, 70);
inline constexpr ImU32 kGood       = IM_COL32(74, 222, 128, 255);   // green
inline constexpr ImU32 kInfo       = IM_COL32(96, 165, 250, 255);   // blue
inline constexpr ImU32 kWarn       = IM_COL32(250, 204, 21, 255);   // amber
inline constexpr ImU32 kDanger     = IM_COL32(248, 113, 113, 255);  // red
inline constexpr ImU32 kHighlight  = IM_COL32(255, 191, 64, 255);   // pulsing highlight

// Surface / text tokens (for cards, panels, depictions drawn via ImDrawList).
inline constexpr ImU32 kBg        = IM_COL32(17, 21, 28, 255);
inline constexpr ImU32 kSurface   = IM_COL32(23, 28, 37, 255);
inline constexpr ImU32 kSurfaceHi = IM_COL32(31, 38, 50, 255);
inline constexpr ImU32 kBorder    = IM_COL32(44, 52, 66, 255);
inline constexpr ImU32 kTextHi    = IM_COL32(226, 232, 240, 255);
inline constexpr ImU32 kTextDim   = IM_COL32(148, 163, 184, 255);

void apply();

ImVec4 verdictColor(int verdict);  // 0=Info 1=Good 2=Warn 3=Danger

// ---- Small shared widgets (consistent look across the navigator + panels) ----

// Accent-tinted small-caps section heading followed by a thin rule. Use to break a
// panel into clearly labelled blocks instead of a wall of TextDisabled().
void sectionHeader(const char* label);

// Rounded "pill" / chip drawn at the current cursor (e.g. the active compound,
// a legal-status tag). Advances the layout by its own size. `bg`/`fg` default to
// the accent palette.
void pill(const char* text, ImU32 bg = kAccentSoft, ImU32 fg = kTextHi);

}  // namespace stimlab::theme
