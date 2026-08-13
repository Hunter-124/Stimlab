# ---- biocad_ui : ImGui app shell + feature panels ----
add_library(biocad_ui STATIC
  Theme.cpp
  AppShell.cpp
  Panels.cpp
)

target_include_directories(biocad_ui PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(biocad_ui PUBLIC
  biocad_contracts
  biocad_bio
  biocad_data
  biocad_core
  biocad_render
  biocad_chem
  biocad_modules
  biocad_workflow
  biocad_agent
  biocad_pkpd
  imgui::imgui
  implot::implot
)

target_compile_features(biocad_ui PUBLIC cxx_std_20)
