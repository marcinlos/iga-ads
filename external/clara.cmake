add_library(clara::clara INTERFACE IMPORTED GLOBAL)

set_target_properties(clara::clara
  PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_LIST_DIR}/clara)
