# ---- stimlab_data : JSON-serializable domain types (the lingua franca) ----
add_library(stimlab_data STATIC
  Domain.cpp
)

target_include_directories(stimlab_data PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(stimlab_data PUBLIC
  nlohmann_json::nlohmann_json
)

target_compile_features(stimlab_data PUBLIC cxx_std_20)
