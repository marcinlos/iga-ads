cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(ADS LANGUAGES C CXX)


# We use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to Release build
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(SYSTEM "${PROJECT_SOURCE_DIR}/lib")


set(SRC src/ads)
set(ads_SRC
  ${SRC}/form_matrix.cpp
  ${SRC}/bspline/bspline.cpp
  ${SRC}/quad/gauss_data.cpp
  ${SRC}/basis_data.cpp
  ${SRC}/simulation/dimension.cpp
  ${SRC}/simulation/simulation_base.cpp
  ${SRC}/simulation/simulation_1d.cpp
  ${SRC}/simulation/simulation_2d.cpp
  ${SRC}/simulation/simulation_3d.cpp)

# --------------------------------------------------------------------
# Libraries
# --------------------------------------------------------------------
set(ads_LIBS)

# Numerical libraries
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)


list(APPEND ads_LIBS ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

# MUMPS
set(USE_MUMPS OFF CACHE BOOL "Use MUMPS solver")

if (USE_MUMPS)
  set(MUMPS_LIBRARIES "dmumps")
  list(APPEND ads_LIBS ${MUMPS_LIBRARIES})
  find_package(MPI)
  include_directories(SYSTEM ${MPI_CXX_INCLUDE_DIRS})
  list(APPEND ads_LIBS ${MPI_CXX_LIBRARIES})
endif()


find_package(Boost REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

# Galois
set(USE_GALOIS OFF CACHE BOOL "Use Galois framework")

if (USE_GALOIS)
  find_package(Galois 6.0 REQUIRED)
  list(APPEND ads_LIBS Galois::shmem)
  list(APPEND ads_SRC ${SRC}/executor/galois.cpp)
endif()

if (CMAKE_COMPILER_IS_GNUCXX)
  set(ignored_warnings "-Wno-narrowing")
  set(warnings "-Wall -Wextra -pedantic -Woverloaded-virtual -Wnon-virtual-dtor ${ignored_warnings}")

  set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} ${optimization} ${warnings} ${std_version}")
endif()

# Target library definition
add_library(ads STATIC ${ads_SRC})
target_link_libraries(ads ${ads_LIBS})

# --------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------

set(COMPILE_TESTS OFF CACHE BOOL "Enable unit tests")

if (COMPILE_TESTS)
  # Unit testing
  enable_testing()

  function(define_test name)
    add_executable(${name}-test ${ARGN})
    set_target_properties(${name}-test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)
    target_link_libraries(${name}-test ads unittest pthread)
    add_test(NAME ${name} COMMAND ${name}-test -v)
  endfunction()

  set(TSRC test)

  define_test(bspline
    ${TSRC}/ads/bspline/bspline_test.cpp
    ${TSRC}/ads/bspline/eval_test.cpp)
  define_test(tensor ${TSRC}/ads/lin/tensor_test.cpp)
  define_test(solver ${TSRC}/ads/lin/banded_solver_test.cpp)
  define_test(dense_solver ${TSRC}/ads/lin/dense_solver_test.cpp)
  define_test(array ${TSRC}/ads/util/multi_array_test.cpp)

endif()

# --------------------------------------------------------------------
# Problems
# --------------------------------------------------------------------

option(SKIP_PROBLEMS "skip compiling example problems" OFF)

function(define_problem name)
  add_executable(${name} ${ARGN})
  target_link_libraries(${name} ads)
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

    define_problem(maxwell
      src/problems/maxwell/main.cpp)

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


endif()
