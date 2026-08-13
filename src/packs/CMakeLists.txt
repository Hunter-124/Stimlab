# ---- biocad_packs : versioned JSON data packs (compounds + receptor targets) ----
add_library(biocad_packs STATIC
  Pack.cpp
)

target_include_directories(biocad_packs PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(biocad_packs PUBLIC
  biocad_contracts
  biocad_data
  biocad_core
  nlohmann_json::nlohmann_json
)

target_compile_features(biocad_packs PUBLIC cxx_std_20)
