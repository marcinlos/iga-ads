cmake_minimum_required(VERSION 3.1)
project(ADS C CXX)


# We use C++14
set (CMAKE_CXX_STANDARD 14)

# Default to Release build
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()

include_directories("${PROJECT_SOURCE_DIR}/src")

set(SRC src/ads)
set(ads_SRC
  ${SRC}/mass_matrix.cpp
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

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Galois
set(USE_GALOIS OFF CACHE BOOL "Use Galois framework")

if (USE_GALOIS)
  find_package(Galois REQUIRED)
  include_directories(${Galois_INCLUDE_DIRS})
  set(CMAKE_CXX_COMPILER ${Galois_CXX_COMPILER})
  set(CMAKE_CXX_FLAGS "${Galois_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")

  list(APPEND ads_LIBS ${Galois_LIBRARIES})
  list(APPEND ads_SRC ${SRC}/executor/galois.cpp)
endif()

if (CMAKE_COMPILER_IS_GNUCXX)
  set(ignored_warnings "-Wno-narrowing -Wno-deprecated-declarations")
  set(warnings "-Wall -Wextra -pedantic ${ignored_warnings}")

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

    define_problem(cahn_hilliard
    src/problems/cahn_hilliard/ch_2d.cpp)

    define_problem(validation
      src/problems/validation/main.cpp)
 endif()

  define_problem(heat_1d
      src/problems/heat/heat_1d.cpp)

  define_problem(heat_3d
    src/problems/heat/heat_3d.cpp)

endif()
