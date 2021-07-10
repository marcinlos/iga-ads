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

# Detect whether ADS is build as top-level project
if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(ADS_IS_TOP_LEVEL TRUE)
else()
  set(ADS_IS_TOP_LEVEL FALSE)
endif()

# Imaginary libraries to propagate settings
add_library(ads-options-public INTERFACE)
add_library(ads-options-private INTERFACE)
target_compile_features(ads-options-public INTERFACE cxx_std_17)

# Ensure BUILD_TESTING is defined
include(CTest)
