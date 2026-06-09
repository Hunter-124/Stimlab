# ---- stimlab_workflow : DAG engine + job/thread pool (WP-D) ----
add_library(stimlab_workflow STATIC
  JobSystem.cpp
  Dag.cpp
)

target_include_directories(stimlab_workflow PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(stimlab_workflow PUBLIC
  stimlab_core
)

target_compile_features(stimlab_workflow PUBLIC cxx_std_20)
