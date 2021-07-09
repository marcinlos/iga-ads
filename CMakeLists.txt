cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(ADS
  VERSION 0.1.0
  DESCRIPTION "Isogeometric Finite Element Method framework"
  LANGUAGES CXX
)

# We use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules/")

# Default to Release build
if (NOT CMAKE_BUILD_TYPE AND NOT GENERATOR_IS_MULTI_CONFIG)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()

# Options
option(ADS_USE_MUMPS "Use MUMPS solver" OFF)
option(ADS_USE_GALOIS "Use Galois framework" OFF)
option(ADS_BUILD_PROBLEMS "Build example problems" ON)
option(ADS_BUILD_TOOLS "Build supporting applications" ON)
option(ADS_BUILD_TESTING "Build tests (even as a subproject)" OFF)
set(ADS_USE_SANITIZERS "" CACHE STRING "List of sanitizers to use (semicolon separated)")

include(CTest)
include(GNUInstallDirs)

# Imaginary libraries to propagate settings
add_library(ads-options-public INTERFACE)
add_library(ads-options-private INTERFACE)
target_compile_features(ads-options-public INTERFACE cxx_std_17)

# Compiler warning options
include(cmake/Warnings.cmake)

# Sanitizers
include(cmake/Sanitizers.cmake)

# Version information
include(cmake/Version.cmake)

# Target library definition
add_library(ADS)
add_library(ADS::ADS ALIAS ADS)
set_target_properties(ADS PROPERTIES OUTPUT_NAME ads)

target_sources(ADS
  PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}/ads/version.cpp"
    src/ads/form_matrix.cpp
    src/ads/bspline/bspline.cpp
    src/ads/quad/gauss_data.cpp
    src/ads/basis_data.cpp
    src/ads/simulation/dimension.cpp
    src/ads/simulation/simulation_base.cpp
    src/ads/simulation/simulation_1d.cpp
    src/ads/simulation/simulation_2d.cpp
    src/ads/simulation/simulation_3d.cpp
)

target_link_libraries(ADS PUBLIC ads-options-public)
target_link_libraries(ADS PRIVATE ads-options-private)

target_include_directories(ADS
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

# --------------------------------------------------------------------
# Libraries
# --------------------------------------------------------------------

add_subdirectory(external)

# Numerical libraries
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_link_libraries(ADS PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

# Boost
find_package(Boost 1.58.0 REQUIRED)
target_link_libraries(ADS PUBLIC Boost::boost)

# fmt
find_package(fmt 7.1 REQUIRED)
target_link_libraries(ADS PUBLIC fmt::fmt)

# MUMPS
if (ADS_USE_MUMPS)
  find_package(MUMPS REQUIRED)
  target_link_libraries(ADS PUBLIC MUMPS::MUMPS)
endif()

# Galois
if (ADS_USE_GALOIS)
  find_package(Galois 6.0 REQUIRED)
  target_link_libraries(ADS PUBLIC Galois::shmem)
  target_sources(ADS PRIVATE src/ads/executor/galois.cpp)
endif()

# Helper for defining programs using ads
include(cmake/AddProgram.cmake)

# --------------------------------------------------------------------
# Other targets
# --------------------------------------------------------------------

# Tests
if ((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR ADS_BUILD_TESTING) AND BUILD_TESTING)
  add_subdirectory(tests)
endif()

# Example problems
if (ADS_BUILD_PROBLEMS)
  add_subdirectory(examples)
endif()

# Tools
if (ADS_BUILD_TOOLS)
  add_subdirectory(tools)
endif()

# --------------------------------------------------------------------
# Installation
# --------------------------------------------------------------------

include(cmake/Installation.cmake)
include(cmake/Packaging.cmake)
