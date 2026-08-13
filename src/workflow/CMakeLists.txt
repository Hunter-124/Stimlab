# ---- biocad_workflow : DAG engine + job/thread pool (WP-D) ----
add_library(biocad_workflow STATIC
  JobSystem.cpp
  Dag.cpp
)

target_include_directories(biocad_workflow PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(biocad_workflow PUBLIC
  biocad_core
)

target_compile_features(biocad_workflow PUBLIC cxx_std_20)
