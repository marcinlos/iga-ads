if (NOT ADS_BUILD_TOOLS)
  return()
endif()

include(AddProgram)

# Custom target that depends on all tools
add_custom_target(ads-tools)

# Helper function prefixing target and executable names with "ads-"
function(add_tool name)
  set(TARGET_NAME "ads-${name}")

  add_program(${TARGET_NAME} ${ARGN})

  if (TARGET ${TARGET_NAME})
    install(TARGETS ${TARGET_NAME} COMPONENT ads-tools)
    add_dependencies(ads-tools ${TARGET_NAME})
  endif()
endfunction()

# --------------------------------------------------------------------
# Tool definitions
# --------------------------------------------------------------------

add_tool(error GALOIS SRC error.cpp)
