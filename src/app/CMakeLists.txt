# ---- StimLab : the single GUI executable (WIN32 subsystem) ----
add_executable(StimLab WIN32
  WinMain.cpp
)

# Compile the Windows resource script (app icon) into the exe.
if(WIN32)
  enable_language(RC)
  target_sources(StimLab PRIVATE StimLab.rc)
endif()

target_link_libraries(StimLab PRIVATE
  stimlab_ui
  stimlab_render
  stimlab_modules
  stimlab_core
  stimlab_data
  imgui::imgui
  implot::implot
)

set_target_properties(StimLab PROPERTIES OUTPUT_NAME "StimLab")
