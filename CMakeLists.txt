cmake_minimum_required(VERSION 3.25)
project(
    ADS
    VERSION 0.1.0
    DESCRIPTION "Isogeometric Finite Element Method framework"
    LANGUAGES CXX
)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules/")

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

option(
    ADS_DEV_SANITIZERS
    "Build with sanitizers enabled. Default: OFF. Values: { ON, OFF }."
    OFF
)

option(ADS_ENABLE_COVERAGE "Code coverage" OFF)

# --------------------------------------------------------------------
# Other settings
# --------------------------------------------------------------------

# Code coverage
include(cmake/Coverage.cmake)

# --------------------------------------------------------------------
# Target library definition
# --------------------------------------------------------------------

add_library(ads)
add_library(ads::ads ALIAS ads)
target_compile_features(ads INTERFACE cxx_std_17)

target_sources(
    ads
    PUBLIC
        FILE_SET HEADERS
            BASE_DIRS #
                include
                "${CMAKE_CURRENT_BINARY_DIR}/include"
)

# --------------------------------------------------------------------
# Dependencies
# --------------------------------------------------------------------

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Boost 1.70 CONFIG REQUIRED)
find_package(fmt 7.1 CONFIG REQUIRED)

target_link_libraries(
    ads
    PUBLIC #
        LAPACK::LAPACK
        BLAS::BLAS
        Boost::boost
        fmt::fmt
)

if(ADS_USE_MUMPS)
    find_package(MUMPS REQUIRED)
    target_link_libraries(ads PUBLIC MUMPS::MUMPS)
endif()

if(ADS_USE_GALOIS)
    find_package(Galois 6.0 CONFIG REQUIRED)
    target_link_libraries(ads PUBLIC Galois::shmem)
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

if(PROJECT_IS_TOP_LEVEL)
    include(cmake/Installation.cmake)
    include(cmake/Packaging.cmake)
endif()
