# ---- stimlab_core : foundation (paths, logging, errors, hashing, events) ----
add_library(stimlab_core STATIC
  AppPaths.cpp
  Log.cpp
  Config.cpp
  Secrets.cpp
  Manifest.cpp
)

target_include_directories(stimlab_core PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/..   # allow #include "core/Xxx.h"
)

target_link_libraries(stimlab_core PUBLIC
  nlohmann_json::nlohmann_json
  spdlog::spdlog
  fmt::fmt
)

if(WIN32)
  # DbgHelp -> minidumps; Crypt32 -> DPAPI secret store; Shell32 -> known folders
  target_link_libraries(stimlab_core PUBLIC DbgHelp Crypt32 Shell32 Ole32)
endif()

target_compile_features(stimlab_core PUBLIC cxx_std_20)
