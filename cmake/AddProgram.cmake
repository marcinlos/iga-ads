function(add_program name)
  cmake_parse_arguments(PROGRAM "MUMPS;GALOIS" "LIBS" "SRC" ${ARGN})

  set(OK TRUE)
  if (PROGRAM_GALOIS AND NOT USE_GALOIS)
    set(OK FALSE)
  endif()
  if (PROGRAM_MUMPS AND NOT USE_MUMPS)
    set(OK FALSE)
  endif()

  if (OK)
    add_executable(${name} ${PROGRAM_SRC})
    target_link_libraries(${name} PRIVATE ads project_options ${PROGRAM_LIBS})
  endif()
endfunction()

