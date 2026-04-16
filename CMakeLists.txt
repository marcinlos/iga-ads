cmake_minimum_required(VERSION 3.25)
project(
    ADS
    VERSION 0.1.0
    DESCRIPTION "Isogeometric Finite Element Method framework"
    LANGUAGES CXX
)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules/")

enable_testing()

include(cmake/version.cmake)
include(cmake/options.cmake)

# Target library definition
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

include(cmake/coverage.cmake)
include(cmake/dependencies.cmake)

# Subdirectories
add_subdirectory(src/ads)
add_subdirectory(include/ads)

if(ADS_BUILD_PROBLEMS)
    add_subdirectory(examples)
endif()

if(ADS_BUILD_TOOLS)
    add_subdirectory(tools)
endif()

if(ADS_BUILD_TESTS)
    add_subdirectory(tests/ads)
endif()

# Installation
if(PROJECT_IS_TOP_LEVEL)
    include(cmake/install.cmake)
    include(cmake/package.cmake)
endif()
