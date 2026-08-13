# ---- biocad_pkpd : pharmacodynamic fits and PK/PD exposure simulation ----
add_library(biocad_pkpd STATIC
  PkEngine.cpp
  Pharmacodynamics.cpp
  Fits.cpp
)

target_include_directories(biocad_pkpd PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(biocad_pkpd PUBLIC
  biocad_numeric
  biocad_contracts
  biocad_data
  biocad_core
)

target_compile_features(biocad_pkpd PUBLIC cxx_std_20)
