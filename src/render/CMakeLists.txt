# ---- biocad_render : DirectX 11 device/swapchain + 3D molecular viewport ----
add_library(biocad_render STATIC
  Dx11Device.cpp
  MolViewport.cpp
)

target_include_directories(biocad_render PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

# MolViewport needs the chem conformer/radii (biocad_chem) and ImGui texture
# types (ImTextureID via ImGui::Image in the panels). d3dcompiler compiles the
# viewport's HLSL at runtime.
target_link_libraries(biocad_render PUBLIC
  biocad_chem
  imgui::imgui
  spdlog::spdlog
)

if(WIN32)
  # windowscodecs = WIC, the PNG encoder behind Dx11Device::captureBackBufferPng.
  target_link_libraries(biocad_render PUBLIC d3d11 dxgi d3dcompiler windowscodecs)
endif()

target_compile_features(biocad_render PUBLIC cxx_std_20)
