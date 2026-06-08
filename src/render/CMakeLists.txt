# ---- stimlab_render : DirectX 11 device/swapchain ----
add_library(stimlab_render STATIC
  Dx11Device.cpp
)

target_include_directories(stimlab_render PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

if(WIN32)
  target_link_libraries(stimlab_render PUBLIC d3d11 dxgi)
endif()

target_compile_features(stimlab_render PUBLIC cxx_std_20)
