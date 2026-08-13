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
  test_sim_random.cpp
  test_sim_nca.cpp
  test_sim_ddi.cpp
  test_bio_score.cpp
  test_bio_align.cpp
  test_bio_conservation.cpp
  test_bio_rotamer.cpp
  test_variant_module.cpp
  test_bio_connectivity.cpp
  test_bio_cartoon.cpp
  test_bio_io.cpp
  test_bio_nuc.cpp
  test_bio_oligo.cpp
  test_bio_codon.cpp
  test_metabolites.cpp
  test_chem_speciation.cpp
  test_chem_solubility.cpp
  test_ionization_module.cpp
  test_nucleic_module.cpp
  test_assay_fits.cpp
  test_assay_dataset.cpp
  test_assay_qc.cpp
  test_assay_module.cpp
  test_assay_biophysics.cpp
  test_network.cpp
  test_kinetics.cpp
  test_sbml.cpp
  test_enrichment.cpp
  test_biologics_numbering.cpp
  test_biologics_mass.cpp
  test_biologics_interface.cpp
  test_mechanism.cpp
)

target_link_libraries(biocad_tests PRIVATE
  Catch2::Catch2WithMain
  biocad_core
  biocad_data
  biocad_numeric
  biocad_bio
  biocad_pkpd
  biocad_sim
  biocad_assay
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
