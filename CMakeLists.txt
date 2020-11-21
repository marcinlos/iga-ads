cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(ADS LANGUAGES C CXX)

# We use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/Modules/")

# Default to Release build
if(NOT CMAKE_BUILD_TYPE)
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
target_include_directories(ads PUBLIC ${PROJECT_SOURCE_DIR}/src)

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


# --------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------
if (COMPILE_TESTS)
  include(CTest)
  add_subdirectory(test)
endif()

# --------------------------------------------------------------------
# Problems
# --------------------------------------------------------------------

function(define_problem name)
  add_executable(${name} ${ARGN})
  target_link_libraries(${name} PRIVATE ads project_options)
endfunction()

if (NOT SKIP_PROBLEMS)

  if (USE_GALOIS)
    define_problem(heat_2d
      src/problems/heat/heat_2d.cpp)

    define_problem(flow
      src/problems/flow/main.cpp)

    define_problem(linear_elasticity
      src/problems/elasticity/main.cpp)

    define_problem(elasticity_victor
      src/problems/elasticity/victor.cpp)

    define_problem(elasticity_pouria
    src/problems/elasticity/elasticity_pouria.cpp)

    define_problem(stokes_split
      src/problems/stokes/main_split.cpp)

    define_problem(validation
      src/problems/validation/main.cpp)

    define_problem(pollution
      src/problems/pollution/main.cpp)

    define_problem(pollution_dpg
      src/problems/pollution/dpg.cpp)

    define_problem(pollution_dpg_3d
      src/problems/pollution/dpg3d.cpp)

    define_problem(victor
      src/problems/victor/main.cpp)

    define_problem(erikkson
      src/problems/erikkson/main.cpp)

    define_problem(erikkson_mumps
      src/problems/erikkson/main_mumps.cpp)

    define_problem(erikkson_nonstationary
      src/problems/erikkson/main_nonstationary.cpp)

    define_problem(stokes
      src/problems/stokes/main.cpp)

    define_problem(stokes_dg
      src/problems/stokes/main_dg.cpp)

    define_problem(stokes_nonstationary
      src/problems/stokes/main_nonstationary.cpp)

    define_problem(stokes_projection
      src/problems/stokes/main_projection.cpp)

    define_problem(demkowicz
      src/problems/demkowicz/main.cpp)

    define_problem(cg
      src/problems/cg/main.cpp
      src/problems/cg/shishkin.cpp)
    target_link_libraries(cg PRIVATE clara::clara)

    define_problem(maxwell
      src/problems/maxwell/main.cpp)
    target_link_libraries(maxwell PRIVATE clara::clara)

    define_problem(maxwell_head
      src/problems/maxwell/main_head.cpp)
    target_link_libraries(maxwell_head PRIVATE clara::clara)

    define_problem(error
      src/tools/error.cpp)

    define_problem(tumor
      src/problems/tumor/vasculature/plot.cpp
      src/problems/tumor/vasculature/vasculature.cpp
      src/problems/tumor/tumor.cpp
      src/problems/tumor.cpp)

    define_problem(tumor_3d
      src/problems/tumor/3d/main.cpp)

  endif()

  define_problem(heat_1d
      src/problems/heat/heat_1d.cpp)

  define_problem(heat_3d
    src/problems/heat/heat_3d.cpp)

  define_problem(implicit
    src/problems/implicit/main.cpp)

  define_problem(coupled
    src/problems/implicit/coupled.cpp)

  define_problem(dg_laplace
    src/problems/dg/laplace.cpp)

endif()
