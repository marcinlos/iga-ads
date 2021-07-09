if (NOT ADS_BUILD_TOOLS)
  return()
endif()

include(AddProgram)

# Helper function prefixing target and executable names with "ads-"
function(add_tool name)
  set(TARGET_NAME "ads-${name}")

  add_program(${TARGET_NAME} ${ARGN})

  if (TARGET ${TARGET_NAME})
    install(TARGETS ${TARGET_NAME} COMPONENT ads-tools)
  endif()
endfunction()

# --------------------------------------------------------------------
# Tool definitions
# --------------------------------------------------------------------

add_tool(error GALOIS SRC error.cpp)
