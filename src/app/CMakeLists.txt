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

# The compound/target catalog is data, so it has to travel with the binary: copy
# assets/packs next to the exe after every build. packs::builtinPackDir() looks
# here first, then falls back to the source tree for a dev run from elsewhere.
add_custom_command(TARGET BioCAD POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          "${CMAKE_SOURCE_DIR}/assets/packs"
          "$<TARGET_FILE_DIR:BioCAD>/assets/packs"
  COMMENT "Copying data packs beside BioCAD.exe")
