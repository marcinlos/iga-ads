cmake_minimum_required(VERSION 2.8)
project(ADS)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

# We use C++11
add_definitions(-std=c++11)

# Coverage
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage -pg")

# gprof
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall -g -pedantic ")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(SRC src/ads)

include_directories("${PROJECT_BINARY_DIR}/src")

add_executable(ads
    src/main.cpp)
    
target_link_libraries(ads lapack blas)

# Unit testing
# enable_testing()

# set(TSRC test/ads)

# add_executable(runUnitTests
#     ${TSRC}/example.cpp)
    
# target_link_libraries(runUnitTests unittest)
    
