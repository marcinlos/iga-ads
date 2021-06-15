
find_package(Catch2 REQUIRED)
include(Catch)

add_executable(suite
    ads/bspline/bspline_test.cpp
    ads/bspline/eval_test.cpp
    ads/util/multi_array_test.cpp
        ads/lin/band_solve_test.cpp
        ads/lin/dense_solve_test.cpp
    ads/lin/tensor_test.cpp
    main.cpp
    )

target_include_directories(suite
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(suite PRIVATE ads project_options Catch2::Catch2)

catch_discover_tests(suite)
