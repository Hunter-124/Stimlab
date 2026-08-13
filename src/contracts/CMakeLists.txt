# ---- biocad_contracts : header-only frozen interfaces ----
add_library(biocad_contracts INTERFACE)

target_include_directories(biocad_contracts INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(biocad_contracts INTERFACE
  biocad_data
  biocad_bio
  biocad_chem
)
