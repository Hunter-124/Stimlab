# ---- stimlab_ui : ImGui app shell + feature panels ----
add_library(stimlab_ui STATIC
  Theme.cpp
  AppShell.cpp
  Panels.cpp
)

target_include_directories(stimlab_ui PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(stimlab_ui PUBLIC
  stimlab_contracts
  stimlab_data
  stimlab_core
  stimlab_render
  stimlab_chem
  imgui::imgui
  implot::implot
)

target_compile_features(stimlab_ui PUBLIC cxx_std_20)
