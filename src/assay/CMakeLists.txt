# ---- biocad_assay : real experimental data - import, QC, fitting, design ----
add_library(biocad_assay STATIC
  Dataset.cpp
  Qc.cpp
  Fits.cpp
  Biophysics.cpp
  Design.cpp
)

target_include_directories(biocad_assay PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(biocad_assay PUBLIC
  biocad_data
  biocad_core
  biocad_numeric
  nlohmann_json::nlohmann_json
)

target_compile_features(biocad_assay PUBLIC cxx_std_20)
