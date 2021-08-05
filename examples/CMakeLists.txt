if (NOT ADS_BUILD_PROBLEMS)
  return()
endif()

include(AddProgram)

# Custom target that depends on all example applications
add_custom_target(ads-examples)

# Helper function prefixing target names with "ads-example-"
function(add_example name)
  set(TARGET_NAME "ads-example-${name}")
  add_program(${TARGET_NAME} ${ARGN})

  if (TARGET ${TARGET_NAME})
    set_target_properties(${TARGET_NAME} PROPERTIES OUTPUT_NAME ${name})
    add_dependencies(ads-examples ${TARGET_NAME})
  endif()
endfunction()

# --------------------------------------------------------------------
# Problem definitions
# --------------------------------------------------------------------

add_example(heat_1d
  SRC
  heat/heat_1d.cpp)

add_example(heat_2d GALOIS
  SRC
  heat/heat_2d.cpp)

add_example(heat_3d
  SRC
  heat/heat_3d.cpp)

add_example(scalability GALOIS
  SRC
  scalability/main.cpp)

add_example(flow GALOIS
  SRC
  flow/main.cpp)

add_example(linear_elasticity GALOIS
  SRC
  elasticity/main.cpp)

add_example(elasticity_victor GALOIS
  SRC
  elasticity/victor.cpp)

add_example(elasticity_pouria GALOIS
  SRC
  elasticity/elasticity_pouria.cpp)

add_example(stokes_split MUMPS GALOIS
  SRC
  stokes/main_split.cpp)

add_example(validation GALOIS
  SRC
  validation/main.cpp)

add_example(pollution GALOIS
  SRC
  pollution/main.cpp)

add_example(pollution_dpg GALOIS
  SRC
  pollution/dpg.cpp)

add_example(pollution_dpg_3d GALOIS
  SRC
  pollution/dpg3d.cpp)

add_example(victor GALOIS
  SRC
  victor/main.cpp)

add_example(erikkson GALOIS
  SRC
  erikkson/main.cpp)

add_example(erikkson_mumps MUMPS GALOIS
  SRC
  erikkson/main_mumps.cpp)

add_example(erikkson_nonstationary MUMPS GALOIS
  SRC
  erikkson/main_nonstationary.cpp)

add_example(multistep GALOIS
  SRC
  multistep/main.cpp
  multistep/scheme.cpp)

add_example(stokes MUMPS GALOIS
  SRC
  stokes/main.cpp)

add_example(stokes_dg MUMPS GALOIS
  SRC
  stokes/main_dg.cpp)

add_example(stokes_nonstationary MUMPS GALOIS
  SRC
  stokes/main_nonstationary.cpp)

add_example(stokes_projection MUMPS GALOIS
  SRC
  stokes/main_projection.cpp)

add_example(demkowicz GALOIS
  SRC
  demkowicz/main.cpp)

add_example(cg MUMPS GALOIS
  SRC
  cg/main.cpp
  cg/shishkin.cpp
  LIBS clara::clara)

add_example(maxwell GALOIS
  SRC
  maxwell/main.cpp
  LIBS clara::clara)

add_example(maxwell_head MUMPS GALOIS
  SRC
  maxwell/main_head.cpp
  maxwell/head_data.cpp
  maxwell/maxwell_head_problem.cpp
  LIBS clara::clara)

add_example(tumor GALOIS
  SRC
  tumor/vasculature/plot.cpp
  tumor/vasculature/vasculature.cpp
  tumor/tumor.cpp
  tumor/tumor_main.cpp)

add_example(tumor_3d GALOIS
  SRC
  tumor/3d/vasculature_parser.cpp
  tumor/3d/main.cpp)

add_example(implicit GALOIS
  SRC
  implicit/main.cpp)

add_example(coupled GALOIS
  SRC
  implicit/coupled.cpp)

add_example(dg_laplace MUMPS GALOIS
  SRC
  dg/laplace.cpp)
