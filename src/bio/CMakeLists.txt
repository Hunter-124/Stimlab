add_library(biocad_bio STATIC
  Structure.cpp
  PdbReader.cpp
  CifReader.cpp
  Fasta.cpp
  Score.cpp
  Align.cpp
  Superpose.cpp
)

target_include_directories(biocad_bio PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_link_libraries(biocad_bio PUBLIC
  biocad_data
  biocad_core
  Eigen3::Eigen
  nlohmann_json::nlohmann_json
)

target_compile_features(biocad_bio PUBLIC cxx_std_20)
