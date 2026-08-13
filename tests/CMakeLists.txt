# ---- Catch2 unit tests ----
find_package(Catch2 CONFIG REQUIRED)

add_executable(biocad_tests
  test_core.cpp
  test_backend.cpp
  test_chem.cpp
  test_chem_groups.cpp
  test_chem_rings.cpp
  test_chem_formula.cpp
  test_chem_canonical.cpp
  test_chem_smarts.cpp
  test_chem_alerts.cpp
  test_chem_crippen.cpp
  test_modules.cpp
  test_embed3d.cpp
  test_molviewport.cpp
  test_docking.cpp
  test_storage.cpp
  test_agent.cpp
  test_workflow.cpp
  test_pipelines.cpp
  test_manifest.cpp
  test_packs.cpp
  test_provenance.cpp
  test_pkpd_fits.cpp
  test_numeric.cpp
  test_pkpd_engine.cpp
  test_bio_score.cpp
  test_bio_align.cpp
  test_bio_io.cpp
  test_metabolites.cpp
  test_chem_speciation.cpp
  test_chem_solubility.cpp
  test_ionization_module.cpp
)

target_link_libraries(biocad_tests PRIVATE
  Catch2::Catch2WithMain
  biocad_core
  biocad_data
  biocad_numeric
  biocad_bio
  biocad_pkpd
  biocad_packs
  biocad_chem
  biocad_modules
  biocad_render
  biocad_storage
  biocad_agent
  biocad_workflow
  biocad_ui
)

target_compile_features(biocad_tests PRIVATE cxx_std_20)

# The BLOSUM62 matrix pack is data on disk; tests read the in-tree copy.
target_compile_definitions(biocad_tests PRIVATE
  BIOCAD_ASSETS_DIR="${CMAKE_SOURCE_DIR}/assets"
  BIOCAD_TEST_FIXTURES="${CMAKE_CURRENT_SOURCE_DIR}/fixtures"
)

include(Catch)
catch_discover_tests(biocad_tests)
