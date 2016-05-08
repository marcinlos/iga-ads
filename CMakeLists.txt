cmake_minimum_required(VERSION 2.8)
project(ADS CXX)


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(SRC src/ads)

include_directories("${PROJECT_SOURCE_DIR}/src")

# We use C++14
set(std_version "-std=c++14")
set(optimization "-O2")
set(ignored_warnings "-Wno-narrowing -Wno-deprecated-declarations")
set(warnings "-Wall -Wextra -pedantic ${ignored_warnings}")
set(debug_opts "-g")

set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} ${debug_opts} ${optimization} ${warnings} ${std_version}")

# Packages
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

list(APPEND numerical_libs lapack blas)

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(ads
    ${SRC}/bspline/bspline.cpp
    ${SRC}/quad/gauss_data.cpp
    src/main.cpp)

target_link_libraries(ads ${numerical_libs})

# Unit testing
enable_testing()

set(TSRC test)

add_executable(runUnitTests
    ${SRC}/bspline/bspline.cpp
    ${TSRC}/ads/bspline/bspline_test.cpp
    ${TSRC}/ads/bspline/eval_test.cpp
    ${TSRC}/ads/lin/tensor_test.cpp
    ${TSRC}/ads/lin/banded_solver_test.cpp
    ${TSRC}/ads/util/multi_array_test.cpp)

target_link_libraries(runUnitTests ${numerical_libs} unittest pthread)

add_test(NAME AllTests COMMAND $<TARGET_FILE:runUnitTests>)

