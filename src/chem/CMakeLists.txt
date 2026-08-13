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
  Alerts.cpp
  Crippen.cpp
)

target_include_directories(biocad_chem PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

# biocad_core supplies the asset-root resolver (core/Assets.h) that the SMARTS
# rule packs are read through; nlohmann_json parses them. This must never become
# a dependency on biocad_packs: biocad_contracts already depends on biocad_chem,
# so that direction would be a cycle.
target_link_libraries(biocad_chem PUBLIC
  Eigen3::Eigen
  biocad_core
  nlohmann_json::nlohmann_json
)

target_compile_features(biocad_chem PUBLIC cxx_std_20)
