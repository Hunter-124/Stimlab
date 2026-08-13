# ---- biocad_numeric : the one least-squares fitter and the one ODE integrator ----
add_library(biocad_numeric STATIC
  Optimize.cpp
  Ode.cpp
)

target_include_directories(biocad_numeric PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(biocad_numeric PUBLIC
  Eigen3::Eigen
)

target_compile_features(biocad_numeric PUBLIC cxx_std_20)
