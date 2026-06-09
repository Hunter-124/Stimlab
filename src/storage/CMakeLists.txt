# ---- stimlab_storage : SQLite-backed persistence (run history) ----
add_library(stimlab_storage STATIC
  RunStore.cpp
)

target_include_directories(stimlab_storage PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(stimlab_storage PUBLIC
  stimlab_contracts
  stimlab_data
  stimlab_core
  SQLite::SQLite3
)

target_compile_features(stimlab_storage PUBLIC cxx_std_20)
