# ---- biocad_fakes : thick fakes implementing the module contracts ----
add_library(biocad_fakes STATIC
  FakeBackend.cpp
)

target_include_directories(biocad_fakes PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(biocad_fakes PUBLIC
  biocad_contracts
  biocad_bio
  biocad_data
  biocad_packs
  biocad_pkpd
)

target_compile_features(biocad_fakes PUBLIC cxx_std_20)
