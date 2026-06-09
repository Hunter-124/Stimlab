# ---- stimlab_modules : REAL analysis backend (chem-engine backed) ----
add_library(stimlab_modules STATIC
  RealBackend.cpp
  docking/Presets.cpp
  docking/PdbqtWriter.cpp
  docking/VinaParser.cpp
  docking/EngineLocator.cpp
  docking/Backends.cpp
)

target_include_directories(stimlab_modules PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(stimlab_modules PUBLIC
  stimlab_contracts
  stimlab_data
  stimlab_chem
  stimlab_core
)

target_compile_features(stimlab_modules PUBLIC cxx_std_20)
