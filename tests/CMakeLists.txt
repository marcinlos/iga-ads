# Only build tests if:
# a) we are top-level project and BUILD_TESTING is TRUE
# b) both BUILD_TESTING and ADS_BUILD_TESTING are TRUE
if (NOT (BUILD_TESTING AND (ADS_IS_TOP_LEVEL OR ADS_BUILD_TESTING)))
  return()
endif()

find_package(Catch2 REQUIRED)
include(Catch)

add_executable(ads-suite)
set_target_properties(ads-suite PROPERTIES OUTPUT_NAME suite)

target_sources(ads-suite
  PRIVATE
    main.cpp
    ads/bspline/bspline_test.cpp
    ads/bspline/eval_test.cpp
    ads/basis_data_test.cpp
    ads/util/multi_array_test.cpp
    ads/lin/band_solve_test.cpp
    ads/lin/dense_solve_test.cpp
    ads/lin/tensor_test.cpp
)

target_include_directories(ads-suite PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(ads-suite PRIVATE ads-objects ads-options-private Catch2::Catch2)

catch_discover_tests(ads-suite PROPERTIES LABELS "ADS")
