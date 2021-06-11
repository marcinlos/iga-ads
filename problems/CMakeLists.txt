
function(define_problem name)
  add_executable(${name} ${ARGN})
  target_link_libraries(${name} PRIVATE ads project_options)
endfunction()

if (USE_GALOIS)
  define_problem(heat_2d
    heat/heat_2d.cpp)

  define_problem(flow
    flow/main.cpp)

  define_problem(linear_elasticity
    elasticity/main.cpp)

  define_problem(elasticity_victor
    elasticity/victor.cpp)

  define_problem(elasticity_pouria
    elasticity/elasticity_pouria.cpp)

  define_problem(stokes_split
    stokes/main_split.cpp)

  define_problem(validation
    validation/main.cpp)

  define_problem(pollution
    pollution/main.cpp)

  define_problem(pollution_dpg
    pollution/dpg.cpp)

  define_problem(pollution_dpg_3d
    pollution/dpg3d.cpp)

  define_problem(victor
    victor/main.cpp)

  define_problem(erikkson
    erikkson/main.cpp)

  define_problem(erikkson_mumps
    erikkson/main_mumps.cpp)

  define_problem(erikkson_nonstationary
    erikkson/main_nonstationary.cpp)

  define_problem(stokes
    stokes/main.cpp)

  define_problem(stokes_dg
    stokes/main_dg.cpp)

  define_problem(stokes_nonstationary
    stokes/main_nonstationary.cpp)

  define_problem(stokes_projection
    stokes/main_projection.cpp)

  define_problem(demkowicz
    demkowicz/main.cpp)

  define_problem(cg
    cg/main.cpp
    cg/shishkin.cpp)
  target_link_libraries(cg PRIVATE clara::clara)

  define_problem(maxwell
    maxwell/main.cpp)
  target_link_libraries(maxwell PRIVATE clara::clara)

  define_problem(maxwell_head
    maxwell/main_head.cpp)
  target_link_libraries(maxwell_head PRIVATE clara::clara)

  define_problem(tumor
    tumor/vasculature/plot.cpp
    tumor/vasculature/vasculature.cpp
    tumor/tumor.cpp
          tumor/tumor_main.cpp)

  define_problem(tumor_3d
    tumor/3d/main.cpp)

endif()

define_problem(heat_1d
  heat/heat_1d.cpp)

define_problem(heat_3d
  heat/heat_3d.cpp)

define_problem(implicit
  implicit/main.cpp)

define_problem(coupled
  implicit/coupled.cpp)

define_problem(dg_laplace
  dg/laplace.cpp)

