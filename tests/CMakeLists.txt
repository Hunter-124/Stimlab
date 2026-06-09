# ---- Catch2 unit tests ----
find_package(Catch2 CONFIG REQUIRED)

add_executable(stimlab_tests
  test_core.cpp
  test_fakes.cpp
  test_chem.cpp
  test_modules.cpp
  test_embed3d.cpp
  test_molviewport.cpp
  test_docking.cpp
  test_storage.cpp
  test_agent.cpp
  test_workflow.cpp
)

target_link_libraries(stimlab_tests PRIVATE
  Catch2::Catch2WithMain
  stimlab_core
  stimlab_data
  stimlab_fakes
  stimlab_chem
  stimlab_modules
  stimlab_render
  stimlab_storage
  stimlab_agent
  stimlab_workflow
)

target_compile_features(stimlab_tests PRIVATE cxx_std_20)

include(Catch)
catch_discover_tests(stimlab_tests)
