# ---- biocad_data : JSON-serializable domain types (the lingua franca) ----
add_library(biocad_data STATIC
  Domain.cpp
)

target_include_directories(biocad_data PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(biocad_data PUBLIC
  nlohmann_json::nlohmann_json
)

target_compile_features(biocad_data PUBLIC cxx_std_20)
