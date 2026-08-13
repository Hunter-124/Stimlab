# ---- biocad_sim : population PK, uncertainty, drug interactions (Phase 13)
#      and reaction-network / systems-biology simulation (Phase 14) ----
#
# One target, two phases, because both are simulation over the SAME seeded random
# stream: a Gillespie trajectory and a population band that disagree about what
# "seed 42" means are not reproducible in the sense either phase claims.
add_library(biocad_sim STATIC
  Random.cpp
  Population.cpp
  Nca.cpp
  Ddi.cpp
  # Phase 14: reaction networks, solvers, chemical kinetics, SBML, flux, control
  # analysis and enrichment. Flux.cpp is compiled unconditionally - the LP and the
  # mass/charge balance gate are ordinary code - and only the MODULE that exposes it
  # is gated on BIOCAD_ENABLE_FBA, so the solver cannot rot behind an off switch.
  Network.cpp
  Solvers.cpp
  Kinetics.cpp
  Sbml.cpp
  Flux.cpp
  Control.cpp
  Enrichment.cpp
)

target_include_directories(biocad_sim PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(biocad_sim PUBLIC
  biocad_numeric
  biocad_pkpd
  biocad_data
  biocad_core
  biocad_chem      # chem::parseFormula, for the flux mass/charge balance
  Eigen3::Eigen
  pugixml::pugixml  # SBML. NOT libSBML: its LGPL relinking obligation is incompatible
                    # with shipping a single fully static Windows binary.
  nlohmann_json::nlohmann_json
)

target_compile_features(biocad_sim PUBLIC cxx_std_20)
