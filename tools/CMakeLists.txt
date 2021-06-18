function(add_tool name)
  add_program(${name} ${ARGN})

  if (TARGET ${name})
    install(TARGETS ${name} COMPONENT ads-tools)
  endif()
endfunction()

add_tool(ads-error GALOIS SRC error.cpp)

