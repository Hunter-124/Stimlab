# ---- biocad_modules : REAL analysis backend (chem-engine backed) ----
add_library(biocad_modules STATIC
  RealBackend.cpp
  BioModules.cpp
  Pipelines.cpp
  Metabolites.cpp
  AlertsModule.cpp
  IonizationModule.cpp
  NucleicModule.cpp
  AssayModule.cpp
  SimulationModule.cpp
  BiologicsModule.cpp
  PopulationPkModule.cpp
  MechanismModule.cpp
  VariantModule.cpp
  docking/Presets.cpp
  docking/PdbqtWriter.cpp
  docking/VinaParser.cpp
  docking/EngineLocator.cpp
  docking/ReceptorPrep.cpp
  docking/Provisioning.cpp
  docking/Backends.cpp
)

target_include_directories(biocad_modules PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(biocad_modules PUBLIC
  biocad_contracts
  biocad_bio
  biocad_assay
  biocad_data
  biocad_packs
  biocad_chem
  biocad_core
  biocad_pkpd
  biocad_sim
  biocad_storage
  biocad_workflow
  # The mechanism module queries public REST APIs (ChEMBL, Reactome) through the ONE
  # HTTP client the assistant's web tools already own (agent::apiGet), so there is a
  # single rate limiter, cache and User-Agent rather than two.
  biocad_agent
)

target_compile_features(biocad_modules PUBLIC cxx_std_20)

# First-party CUDA GPU docking backend (CudaBackend.cpp host + CudaScore.cu kernel),
# compiled only in the CUDA build. The kernel is pure POD so nvcc never sees C++20.
if(BIOCAD_ENABLE_CUDA)
  target_sources(biocad_modules PRIVATE
    docking/CudaBackend.cpp
    docking/CudaScore.cu)
  # Link the CUDA runtime that MATCHES the C runtime: the static cudart for the fully
  # static /MT build (so the GPU exe is self-contained and needs no cudart64_*.dll on
  # PATH), the import-lib + cudart64_*.dll otherwise. Keyed off the vcpkg triplet, which
  # encodes static vs dynamic; CUDA::cudart_static is itself built against the static CRT,
  # so it pairs with the /MT host the windows-cuda-static preset selects.
  if(VCPKG_TARGET_TRIPLET MATCHES "static")
    target_link_libraries(biocad_modules PUBLIC CUDA::cudart_static)
  else()
    target_link_libraries(biocad_modules PUBLIC CUDA::cudart)
  endif()
endif()
