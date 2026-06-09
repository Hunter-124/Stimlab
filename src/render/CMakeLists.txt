# ---- stimlab_render : DirectX 11 device/swapchain + 3D molecular viewport ----
add_library(stimlab_render STATIC
  Dx11Device.cpp
  MolViewport.cpp
)

target_include_directories(stimlab_render PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

# MolViewport needs the chem conformer/radii (stimlab_chem) and ImGui texture
# types (ImTextureID via ImGui::Image in the panels). d3dcompiler compiles the
# viewport's HLSL at runtime.
target_link_libraries(stimlab_render PUBLIC
  stimlab_chem
  imgui::imgui
)

if(WIN32)
  target_link_libraries(stimlab_render PUBLIC d3d11 dxgi d3dcompiler)
endif()

target_compile_features(stimlab_render PUBLIC cxx_std_20)
