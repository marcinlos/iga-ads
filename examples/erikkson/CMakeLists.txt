add_example(erikkson GALOIS SRC main.cpp)

add_example(erikkson_mumps MUMPS GALOIS SRC main_mumps.cpp)

add_example(
    erikkson_inverse
    MUMPS
    GALOIS
    SRC main_inverse.cpp inverse.cpp
    LIBS bfg::lyra
)

add_example(
    inverse_postprocess
    MUMPS
    GALOIS
    SRC inverse_postprocess.cpp inverse.cpp
    LIBS bfg::lyra
)

add_example(erikkson_nonstationary MUMPS GALOIS SRC main_nonstationary.cpp)
