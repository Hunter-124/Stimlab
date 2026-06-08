cmake_minimum_required(VERSION 3.27)

project(StimLab
  VERSION 0.1.0
  DESCRIPTION "Native CNS-stimulant computational drug-discovery suite"
  LANGUAGES CXX)

# ---- Global options -------------------------------------------------------
option(STIMLAB_ENABLE_CUDA    "Enable CUDA-accelerated engines (needs nvcc)" OFF)
option(STIMLAB_ENABLE_SCIENCE "Link the heavy science stack (RDKit/curl)"    OFF)
option(STIMLAB_BUILD_TESTS    "Build the Catch2 test suite"                  ON)

# ---- Language / compiler standard ----------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

if(MSVC)
  add_compile_options(/W4 /permissive- /Zc:__cplusplus /Zc:preprocessor /utf-8 /EHsc /MP)
  add_compile_definitions(NOMINMAX WIN32_LEAN_AND_MEAN UNICODE _UNICODE
                          _CRT_SECURE_NO_WARNINGS
                          _SILENCE_ALL_CXX20_DEPRECATION_WARNINGS
                          _SILENCE_ALL_CXX23_DEPRECATION_WARNINGS)
endif()

# ---- Feature defines passed to code --------------------------------------
if(STIMLAB_ENABLE_SCIENCE)
  add_compile_definitions(STIMLAB_HAVE_SCIENCE=1)
endif()
if(STIMLAB_ENABLE_CUDA)
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  add_compile_definitions(STIMLAB_HAVE_CUDA=1)
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

# ---- Project modules ------------------------------------------------------
add_subdirectory(src/core)
add_subdirectory(src/data)
add_subdirectory(src/chem)
add_subdirectory(src/contracts)
add_subdirectory(src/fakes)
add_subdirectory(src/modules)
add_subdirectory(src/render)
add_subdirectory(src/ui)
add_subdirectory(src/app)

if(STIMLAB_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

message(STATUS "StimLab ${PROJECT_VERSION}  |  CUDA=${STIMLAB_ENABLE_CUDA}  Science=${STIMLAB_ENABLE_SCIENCE}  Tests=${STIMLAB_BUILD_TESTS}")
