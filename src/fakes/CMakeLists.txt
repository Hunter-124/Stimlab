# ---- stimlab_fakes : thick fakes implementing the module contracts ----
add_library(stimlab_fakes STATIC
  FakeBackend.cpp
)

target_include_directories(stimlab_fakes PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(stimlab_fakes PUBLIC
  stimlab_contracts
  stimlab_data
)

target_compile_features(stimlab_fakes PUBLIC cxx_std_20)
