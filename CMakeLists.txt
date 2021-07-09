cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(ADS
  VERSION 0.1.0
  DESCRIPTION "Isogeometric Finite Element Method framework"
  LANGUAGES CXX
)

# General project-wide settings
include(cmake/Setup.cmake)

# --------------------------------------------------------------------
# Build options
# --------------------------------------------------------------------

option(ADS_USE_MUMPS "Use MUMPS solver" OFF)
option(ADS_USE_GALOIS "Use Galois framework" OFF)
option(ADS_BUILD_PROBLEMS "Build example problems" ON)
option(ADS_BUILD_TOOLS "Build supporting applications" ON)
option(ADS_BUILD_TESTING "Build tests (even as a subproject)" OFF)
set(ADS_USE_SANITIZERS "" CACHE STRING "List of sanitizers to use (semicolon separated)")

# --------------------------------------------------------------------
# Other settings
# --------------------------------------------------------------------

# Compiler warning options
include(cmake/Warnings.cmake)

# Sanitizers
include(cmake/Sanitizers.cmake)

# Version information
include(cmake/Version.cmake)

# --------------------------------------------------------------------
# Target library definition
# --------------------------------------------------------------------

add_library(ADS)
add_library(ADS::ADS ALIAS ADS)
set_target_properties(ADS PROPERTIES OUTPUT_NAME ads)

# Apply project settings
target_link_libraries(ADS PUBLIC ads-options-public)
target_link_libraries(ADS PRIVATE ads-options-private)

target_include_directories(ADS
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

# --------------------------------------------------------------------
# Dependencies
# --------------------------------------------------------------------

# Bundled with the project
add_subdirectory(external)

# Numerical libraries
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_link_libraries(ADS PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

find_package(Boost 1.58.0 REQUIRED)
target_link_libraries(ADS PUBLIC Boost::boost)

find_package(fmt 7.1 REQUIRED)
target_link_libraries(ADS PUBLIC fmt::fmt)

if (ADS_USE_MUMPS)
  find_package(MUMPS REQUIRED)
  target_link_libraries(ADS PUBLIC MUMPS::MUMPS)
endif()

if (ADS_USE_GALOIS)
  find_package(Galois 6.0 REQUIRED)
  target_link_libraries(ADS PUBLIC Galois::shmem)
  target_sources(ADS PRIVATE src/ads/executor/galois.cpp)
endif()

# --------------------------------------------------------------------
# Subdirectories
# --------------------------------------------------------------------

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(examples)
add_subdirectory(tools)

# --------------------------------------------------------------------
# Installation & distribution
# --------------------------------------------------------------------

include(cmake/Installation.cmake)
include(cmake/Packaging.cmake)
