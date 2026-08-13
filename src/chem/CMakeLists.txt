# ---- biocad_chem : in-house cheminformatics core ----
# (no chemistry deps; Eigen3 is a header-only linear-algebra helper for Embed3D)
add_library(biocad_chem STATIC
  Periodic.cpp
  Smiles.cpp
  Descriptors.cpp
  Analysis.cpp
  Embed3D.cpp
  Rings.cpp
  Aromaticity.cpp
  Canonical.cpp
  Smarts.cpp
)

target_include_directories(biocad_chem PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(biocad_chem PUBLIC
  Eigen3::Eigen
)

target_compile_features(biocad_chem PUBLIC cxx_std_20)
