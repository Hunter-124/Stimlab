# ---- stimlab_contracts : header-only frozen interfaces ----
add_library(stimlab_contracts INTERFACE)

target_include_directories(stimlab_contracts INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(stimlab_contracts INTERFACE
  stimlab_data
  stimlab_chem
)
