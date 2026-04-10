include(AddExample)

add_example(stokes_split MUMPS GALOIS SRC main_split.cpp)
add_example(stokes MUMPS GALOIS SRC main.cpp)
add_example(stokes_dg MUMPS GALOIS SRC main_dg.cpp)
add_example(stokes_nonstationary MUMPS GALOIS SRC main_nonstationary.cpp)
add_example(stokes_projection MUMPS GALOIS SRC main_projection.cpp)
