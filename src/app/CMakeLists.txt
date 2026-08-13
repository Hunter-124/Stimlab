# ---- BioCAD : the single GUI executable (WIN32 subsystem) ----
add_executable(BioCAD WIN32
  WinMain.cpp
)

# Compile the Windows resource script (app icon) into the exe.
if(WIN32)
  enable_language(RC)
  target_sources(BioCAD PRIVATE BioCAD.rc)
endif()

target_link_libraries(BioCAD PRIVATE
  biocad_ui
  biocad_render
  biocad_modules
  biocad_core
  biocad_data
  imgui::imgui
  implot::implot
)

set_target_properties(BioCAD PROPERTIES OUTPUT_NAME "BioCAD")
