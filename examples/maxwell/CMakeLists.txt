add_library(
    ads.example.maxwell.common
    OBJECT
    state.cpp
    spaces.cpp
    results.cpp
    maxwell_base.cpp
)

target_link_libraries(
    ads.example.maxwell.common
    PUBLIC #
        ADS::ADS
        ads-options-private
        bfg::lyra
)

add_example(
    maxwell_uniform
    GALOIS
    SRC main_uniform.cpp
    LIBS ads.example.maxwell.common
)

add_example(
    maxwell_head
    MUMPS
    GALOIS
    SRC #
        main_head.cpp
        head_data.cpp
        maxwell_head_problem.cpp
    LIBS ads.example.maxwell.common
)

add_example(
    maxwell_ads
    GALOIS
    SRC #
        main_ads.cpp
        head_data.cpp
        maxwell_head_problem.cpp
    LIBS ads.example.maxwell.common
)
