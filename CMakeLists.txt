cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(ADS
  VERSION 0.1.0
  LANGUAGES CXX)

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

include(CTest)
include(GNUInstallDirs)

# Imaginary library to propagate settings
add_library(project-options INTERFACE)
target_compile_features(project-options INTERFACE cxx_std_17)

# Compiler warning options
include(cmake/Warnings.cmake)

# Sanitizers
include(cmake/Sanitizers.cmake)

# Version information
include(cmake/Version.cmake)

# Target library definition
add_library(ads)
add_library(ads::ads ALIAS ads)

target_sources(ads
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

target_link_libraries(ads PUBLIC project-options)
target_link_libraries(ads PRIVATE project-warnings)

target_include_directories(ads
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

# --------------------------------------------------------------------
# Libraries
# --------------------------------------------------------------------

add_subdirectory(libs)

# Numerical libraries
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_link_libraries(ads PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

# Boost
find_package(Boost 1.58.0 REQUIRED)
target_link_libraries(ads PUBLIC Boost::boost)

# fmt
find_package(fmt REQUIRED)
target_link_libraries(ads PUBLIC fmt::fmt)

# MUMPS
if (ADS_USE_MUMPS)
  find_package(MUMPS REQUIRED)
  target_link_libraries(ads PUBLIC MUMPS::MUMPS)
endif()

# Galois
if (ADS_USE_GALOIS)
  find_package(Galois 6.0 REQUIRED)
  target_link_libraries(ads PUBLIC Galois::shmem)
  target_sources(ads PRIVATE src/ads/executor/galois.cpp)
endif()

# Helper for defining programs using ads
include(cmake/AddProgram.cmake)

# --------------------------------------------------------------------
# Other targets
# --------------------------------------------------------------------

# Tests
if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  add_subdirectory(test)
endif()

# Problems
if (ADS_BUILD_PROBLEMS)
  add_subdirectory(problems)
endif()

# Tools
add_subdirectory(tools)

# --------------------------------------------------------------------
# Installation
# --------------------------------------------------------------------

install(
  TARGETS ads project-options project-warnings
  EXPORT ads-targets
  LIBRARY
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT ads-runtime
  ARCHIVE
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT ads-devel
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
  EXPORT ads-targets
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/ads"
  NAMESPACE ads::
  COMPONENT ads-devel
)

configure_file(cmake/ads-config.cmake.in ads-config.cmake @ONLY)

install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/ads-config.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/ads"
  COMPONENT ads-devel
)

install(
  DIRECTORY cmake/Modules
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/ads"
  COMPONENT ads-devel
)

install(
  DIRECTORY include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT ads-devel
)

