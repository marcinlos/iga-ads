cmake_minimum_required(VERSION 3.23 FATAL_ERROR)
project(ADS
  VERSION 0.1.0
  DESCRIPTION "Isogeometric Finite Element Method framework"
  LANGUAGES CXX
)

# General project-wide settings
include(cmake/Setup.cmake)

# Version information
include(cmake/Version.cmake)

# --------------------------------------------------------------------
# Build options
# --------------------------------------------------------------------

option(ADS_USE_MUMPS "Use MUMPS solver. Values: { ON, OFF }." OFF)
option(ADS_USE_GALOIS "Use Galois framework. Values: { ON, OFF }." OFF)

option(
  ADS_BUILD_PROBLEMS
  "Build example problems. Default: ${PROJECT_IS_TOP_LEVEL}. Values: { ON, OFF }."
  ${PROJECT_IS_TOP_LEVEL}
)

option(
  ADS_BUILD_TOOLS
  "Build supporting applications. Default: ${PROJECT_IS_TOP_LEVEL}. Values: { ON, OFF }."
  ${PROJECT_IS_TOP_LEVEL}
)

option(
  ADS_BUILD_TESTS
  "Enable building tests. Default: ${PROJECT_IS_TOP_LEVEL}. Values: { ON, OFF }."
  ${PROJECT_IS_TOP_LEVEL}
)

option(ADS_ENABLE_COVERAGE "Code coverage" OFF)
set(ADS_USE_SANITIZERS "" CACHE STRING "List of sanitizers to use (semicolon separated)")

# --------------------------------------------------------------------
# Other settings
# --------------------------------------------------------------------

# Compiler warning options
include(cmake/Warnings.cmake)

# Sanitizers
include(cmake/Sanitizers.cmake)

# Code coverage
include(cmake/Coverage.cmake)

# --------------------------------------------------------------------
# Target library definition
# --------------------------------------------------------------------

add_library(ADS)
add_library(ADS::ADS ALIAS ADS)
set_target_properties(ADS PROPERTIES OUTPUT_NAME ads)

add_library(ads-objects OBJECT)
target_link_libraries(ADS PUBLIC ads-objects)

# Apply project settings
target_link_libraries(ads-objects PUBLIC ads-options-public)
target_link_libraries(ads-objects PRIVATE ads-options-private)

target_sources(ads-objects
  PUBLIC
    FILE_SET HEADERS
    BASE_DIRS
      include
      ${CMAKE_CURRENT_BINARY_DIR}/include
)

# --------------------------------------------------------------------
# Dependencies
# --------------------------------------------------------------------

# Numerical libraries
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_link_libraries(ads-objects PUBLIC LAPACK::LAPACK BLAS::BLAS)

find_package(Boost 1.70 CONFIG REQUIRED)
target_link_libraries(ads-objects PUBLIC Boost::boost)

find_package(fmt 7.1 REQUIRED)
target_link_libraries(ads-objects PUBLIC fmt::fmt)

if (ADS_USE_MUMPS)
  find_package(MUMPS REQUIRED)
  target_link_libraries(ads-objects PUBLIC MUMPS::MUMPS)
endif()

if (ADS_USE_GALOIS)
  find_package(Galois 6.0 REQUIRED)
  target_link_libraries(ads-objects PUBLIC Galois::shmem)
endif()

# --------------------------------------------------------------------
# Subdirectories
# --------------------------------------------------------------------

add_subdirectory(src/ads)
add_subdirectory(include/ads)

if(ADS_BUILD_PROBLEMS)
    add_subdirectory(examples)
endif()

if(ADS_BUILD_TOOLS)
    add_subdirectory(tools)
endif()

if(ADS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests/ads)
endif()

# --------------------------------------------------------------------
# Installation & distribution
# --------------------------------------------------------------------

include(cmake/Installation.cmake)
include(cmake/Packaging.cmake)
