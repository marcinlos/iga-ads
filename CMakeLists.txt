cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(ADS LANGUAGES C CXX)

# We use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/Modules/")

# Default to Release build
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()

# Options
option(USE_MUMPS "Use MUMPS solver")
option(USE_GALOIS "Use Galois framework")
option(SKIP_PROBLEMS "Skip compiling example problems")
option(COMPILE_TESTS "Enable unit tests")

# Imaginary library to propagate settings
add_library(project_options INTERFACE)

# Compiler warning options
include(cmake/Warnings.cmake)

# Sanitizers
# Enable address sanitizer
add_compile_options("$<$<CONFIG:Debug>:-fsanitize=address,undefined;-fno-omit-frame-pointer>")
add_link_options("$<$<CONFIG:Debug>:-fsanitize=address,undefined>")


# Target library definition
add_library(ads
  STATIC
  src/ads/form_matrix.cpp
  src/ads/bspline/bspline.cpp
  src/ads/quad/gauss_data.cpp
  src/ads/basis_data.cpp
  src/ads/simulation/dimension.cpp
  src/ads/simulation/simulation_base.cpp
  src/ads/simulation/simulation_1d.cpp
  src/ads/simulation/simulation_2d.cpp
  src/ads/simulation/simulation_3d.cpp)

target_link_libraries(ads PRIVATE project_options)
target_include_directories(ads PUBLIC ${PROJECT_SOURCE_DIR}/include)

# --------------------------------------------------------------------
# Libraries
# --------------------------------------------------------------------

include(libs/CMakeLists.txt)

# Numerical libraries
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_link_libraries(ads PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

# Boost
find_package(Boost REQUIRED)
target_link_libraries(ads PUBLIC Boost::boost)

# fmt
find_package(fmt REQUIRED)
target_link_libraries(ads PUBLIC fmt::fmt)

# MUMPS
if (USE_MUMPS)
  find_package(MUMPS)
  target_link_libraries(ads PUBLIC MUMPS::MUMPS)
endif()

# Galois
if (USE_GALOIS)
  find_package(Galois 6.0 REQUIRED)
  target_link_libraries(ads PUBLIC Galois::shmem)
  target_sources(ads PRIVATE src/ads/executor/galois.cpp)
endif()

# Helper for defining programs using ads
include(cmake/AddProgram.cmake)

# Tests
if (COMPILE_TESTS)
  include(CTest)
  add_subdirectory(test)
endif()

# Problems
if (NOT SKIP_PROBLEMS)
  add_subdirectory(problems)
endif()

# Tools
add_program(error GALOIS SRC src/tools/error.cpp)

