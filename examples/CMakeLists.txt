
add_program(heat_1d
  SRC
  heat/heat_1d.cpp)

add_program(heat_2d GALOIS
  SRC
  heat/heat_2d.cpp)

add_program(heat_3d
  SRC
  heat/heat_3d.cpp)

add_program(flow GALOIS
  SRC
  flow/main.cpp)

add_program(linear_elasticity GALOIS
  SRC
  elasticity/main.cpp)

add_program(elasticity_victor GALOIS
  SRC
  elasticity/victor.cpp)

add_program(elasticity_pouria GALOIS
  SRC
  elasticity/elasticity_pouria.cpp)

add_program(stokes_split MUMPS GALOIS
  SRC
  stokes/main_split.cpp)

add_program(validation GALOIS
  SRC
  validation/main.cpp)

add_program(pollution GALOIS
  SRC
  pollution/main.cpp)

add_program(pollution_dpg GALOIS
  SRC
  pollution/dpg.cpp)

add_program(pollution_dpg_3d GALOIS
  SRC
  pollution/dpg3d.cpp)

add_program(victor GALOIS
  SRC
  victor/main.cpp)

add_program(erikkson GALOIS
  SRC
  erikkson/main.cpp)

add_program(erikkson_mumps MUMPS GALOIS
  SRC
  erikkson/main_mumps.cpp)

add_program(erikkson_nonstationary MUMPS GALOIS
  SRC
  erikkson/main_nonstationary.cpp)

add_program(stokes MUMPS GALOIS
  SRC
  stokes/main.cpp)

add_program(stokes_dg MUMPS GALOIS
  SRC
  stokes/main_dg.cpp)

add_program(stokes_nonstationary MUMPS GALOIS
  SRC
  stokes/main_nonstationary.cpp)

add_program(stokes_projection MUMPS GALOIS
  SRC
  stokes/main_projection.cpp)

add_program(demkowicz GALOIS
  SRC
  demkowicz/main.cpp)

add_program(cg MUMPS GALOIS
  SRC
  cg/main.cpp
  cg/shishkin.cpp
  LIBS clara::clara)

add_program(maxwell GALOIS
  SRC
  maxwell/main.cpp
  LIBS clara::clara)

add_program(maxwell_head MUMPS GALOIS
  SRC
  maxwell/main_head.cpp
  LIBS clara::clara)

add_program(tumor GALOIS
  SRC
  tumor/vasculature/plot.cpp
  tumor/vasculature/vasculature.cpp
  tumor/tumor.cpp
  tumor/tumor_main.cpp)

add_program(tumor_3d GALOIS
  SRC
  tumor/3d/main.cpp)

add_program(implicit GALOIS
  SRC
  implicit/main.cpp)

add_program(coupled GALOIS
  SRC
  implicit/coupled.cpp)

add_program(dg_laplace MUMPS GALOIS
  SRC
  dg/laplace.cpp)

