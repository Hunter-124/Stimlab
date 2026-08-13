# ---- biocad_storage : SQLite-backed persistence (run history) ----
add_library(biocad_storage STATIC
  RunStore.cpp
)

target_include_directories(biocad_storage PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(biocad_storage PUBLIC
  biocad_contracts
  biocad_data
  biocad_core
  SQLite::SQLite3
)

target_compile_features(biocad_storage PUBLIC cxx_std_20)
