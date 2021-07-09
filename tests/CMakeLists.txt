
find_package(Catch2 REQUIRED)
include(Catch)

add_executable(ads-suite)
set_target_properties(ads-suite PROPERTIES OUTPUT_NAME suite)

target_sources(ads-suite
  PRIVATE
    main.cpp
    ads/bspline/bspline_test.cpp
    ads/bspline/eval_test.cpp
    ads/util/multi_array_test.cpp
    ads/lin/band_solve_test.cpp
    ads/lin/dense_solve_test.cpp
    ads/lin/tensor_test.cpp
)

target_include_directories(ads-suite PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(ads-suite PRIVATE ADS ads-options-private Catch2::Catch2)

catch_discover_tests(ads-suite PROPERTIES LABELS "ADS")
