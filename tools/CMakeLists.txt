function(add_tool name)
  set(TARGET_NAME "ads-${name}")

  add_program(${TARGET_NAME} ${ARGN})

  if (TARGET ${TARGET_NAME})
    install(TARGETS ${TARGET_NAME} COMPONENT ads-tools)
  endif()
endfunction()

add_tool(error GALOIS SRC error.cpp)

