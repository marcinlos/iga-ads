include(AddExample)

add_example(
    cg
    MUMPS
    GALOIS
    SRC #
        main.cpp
        shishkin.cpp
    LIBS bfg::lyra
)
