cmake_minimum_required(VERSION 3.27)

project(BioCAD
  VERSION 0.2.0
  DESCRIPTION "BioCAD - native workstation for molecular, protein, and pharmacological analysis"
  LANGUAGES CXX)

# ---- Global options -------------------------------------------------------
option(BIOCAD_ENABLE_CUDA    "Enable CUDA-accelerated engines (needs nvcc)" OFF)
option(BIOCAD_ENABLE_SCIENCE "Link the heavy science stack (RDKit/curl)"    OFF)
option(BIOCAD_BUILD_TESTS    "Build the Catch2 test suite"                  ON)
option(BIOCAD_ENABLE_WEBVIEW2 "Optional WebView2 headless JS-render web_fetch (Windows; needs the Evergreen WebView2 Runtime at runtime)" OFF)

# ---- Language / compiler standard ----------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

if(MSVC)
  # Scope the C++ warning/conformance flags to CXX so nvcc (CUDA) is not handed MSVC-only
  # flags like /permissive- that it does not understand.
  add_compile_options(
    "$<$<COMPILE_LANGUAGE:CXX>:/W4>"
    "$<$<COMPILE_LANGUAGE:CXX>:/permissive->"
    "$<$<COMPILE_LANGUAGE:CXX>:/Zc:__cplusplus>"
    "$<$<COMPILE_LANGUAGE:CXX>:/Zc:preprocessor>"
    "$<$<COMPILE_LANGUAGE:CXX>:/utf-8>"
    "$<$<COMPILE_LANGUAGE:CXX>:/EHsc>"
    "$<$<COMPILE_LANGUAGE:CXX>:/MP>")
  add_compile_definitions(NOMINMAX WIN32_LEAN_AND_MEAN UNICODE _UNICODE
                          _CRT_SECURE_NO_WARNINGS
                          _SILENCE_ALL_CXX20_DEPRECATION_WARNINGS
                          _SILENCE_ALL_CXX23_DEPRECATION_WARNINGS)
endif()

# ---- Feature defines passed to code --------------------------------------
if(BIOCAD_ENABLE_SCIENCE)
  add_compile_definitions(BIOCAD_HAVE_SCIENCE=1)
endif()
if(BIOCAD_ENABLE_WEBVIEW2)
  # The rendered fetch augments the curl web tools, so it requires the science feature.
  if(NOT BIOCAD_ENABLE_SCIENCE)
    message(FATAL_ERROR "BIOCAD_ENABLE_WEBVIEW2 requires BIOCAD_ENABLE_SCIENCE.")
  endif()
  add_compile_definitions(BIOCAD_HAVE_WEBVIEW2=1)
endif()
if(BIOCAD_ENABLE_CUDA)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 86)  # RTX 3080 Ti = sm_86 (CUDA 13 dropped pre-sm_75)
  endif()
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  add_compile_definitions(BIOCAD_HAVE_CUDA=1)
endif()

# ---- Third-party dependencies (vcpkg) ------------------------------------
find_package(nlohmann_json CONFIG REQUIRED)
find_package(spdlog        CONFIG REQUIRED)
find_package(fmt           CONFIG REQUIRED)
find_package(SQLite3       REQUIRED)
find_package(Eigen3        CONFIG REQUIRED)
find_package(imgui         CONFIG REQUIRED)
find_package(implot        CONFIG REQUIRED)
find_package(yaml-cpp      CONFIG REQUIRED)

# yaml-cpp target name differs across versions; normalize to yaml-cpp::yaml-cpp
if(NOT TARGET yaml-cpp::yaml-cpp AND TARGET yaml-cpp)
  add_library(yaml-cpp::yaml-cpp ALIAS yaml-cpp)
endif()

# ---- Generated version header (single source: project(BioCAD VERSION ...)) --
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/src/core/Version.h.in"
  "${CMAKE_BINARY_DIR}/generated/core/Version.h"
  @ONLY)

# The data packs must sit beside every binary that reads them (the app and the test
# runner share ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}). Copy at configure time so a test
# build that never builds the exe still finds a catalog; src/app re-copies after each
# build so edits during development take effect.
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/assets/packs"
     DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/assets")

# ---- Project modules ------------------------------------------------------
add_subdirectory(src/core)
add_subdirectory(src/data)
add_subdirectory(src/numeric)
add_subdirectory(src/chem)
add_subdirectory(src/contracts)
add_subdirectory(src/bio)
add_subdirectory(src/pkpd)
add_subdirectory(src/packs)
add_subdirectory(src/storage)
add_subdirectory(src/workflow)
add_subdirectory(src/agent)
add_subdirectory(src/modules)
add_subdirectory(src/render)
add_subdirectory(src/ui)
add_subdirectory(src/app)

if(BIOCAD_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

message(STATUS "BioCAD ${PROJECT_VERSION}  |  CUDA=${BIOCAD_ENABLE_CUDA}  Science=${BIOCAD_ENABLE_SCIENCE}  Tests=${BIOCAD_BUILD_TESTS}")
