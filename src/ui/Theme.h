// ui/Theme.h - one clean, modern dark theme + shared accent colors.
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

void apply();

ImVec4 verdictColor(int verdict);  // 0=Info 1=Good 2=Warn 3=Danger

}  // namespace stimlab::theme
