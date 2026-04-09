find_package(Catch2 3.8 CONFIG REQUIRED)

add_executable(ads.tests)

target_sources(ads.tests
  PRIVATE
    bspline/bspline_test.cpp
    bspline/eval_test.cpp
    basis_data_test.cpp
    util/multi_array_test.cpp
    lin/band_solve_test.cpp
    lin/dense_solve_test.cpp
    lin/tensor_test.cpp
    solver_test.cpp
)

target_link_libraries(ads.tests PRIVATE ADS::ADS  Catch2::Catch2WithMain)

include(Catch)
catch_discover_tests(
    ads.tests
    EXTRA_ARGS --colour-mode ansi
    DISCOVERY_MODE PRE_TEST
)
