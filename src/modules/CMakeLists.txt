# ---- stimlab_modules : REAL analysis backend (chem-engine backed) ----
add_library(stimlab_modules STATIC
  RealBackend.cpp
  Pipelines.cpp
  docking/Presets.cpp
  docking/PdbqtWriter.cpp
  docking/VinaParser.cpp
  docking/EngineLocator.cpp
  docking/ReceptorPrep.cpp
  docking/Provisioning.cpp
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
  stimlab_storage
  stimlab_workflow
)

target_compile_features(stimlab_modules PUBLIC cxx_std_20)

# First-party CUDA GPU docking backend (CudaBackend.cpp host + CudaScore.cu kernel),
# compiled only in the CUDA build. The kernel is pure POD so nvcc never sees C++20.
if(STIMLAB_ENABLE_CUDA)
  target_sources(stimlab_modules PRIVATE
    docking/CudaBackend.cpp
    docking/CudaScore.cu)
  # Link the CUDA runtime that MATCHES the C runtime: the static cudart for the fully
  # static /MT build (so the GPU exe is self-contained and needs no cudart64_*.dll on
  # PATH), the import-lib + cudart64_*.dll otherwise. Keyed off the vcpkg triplet, which
  # encodes static vs dynamic; CUDA::cudart_static is itself built against the static CRT,
  # so it pairs with the /MT host the windows-cuda-static preset selects.
  if(VCPKG_TARGET_TRIPLET MATCHES "static")
    target_link_libraries(stimlab_modules PUBLIC CUDA::cudart_static)
  else()
    target_link_libraries(stimlab_modules PUBLIC CUDA::cudart)
  endif()
endif()
