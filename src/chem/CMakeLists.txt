# ---- stimlab_chem : in-house cheminformatics core ----
# (no chemistry deps; Eigen3 is a header-only linear-algebra helper for Embed3D)
add_library(stimlab_chem STATIC
  Periodic.cpp
  Smiles.cpp
  Descriptors.cpp
  Analysis.cpp
  Embed3D.cpp
)

target_include_directories(stimlab_chem PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(stimlab_chem PUBLIC
  Eigen3::Eigen
)

target_compile_features(stimlab_chem PUBLIC cxx_std_20)
