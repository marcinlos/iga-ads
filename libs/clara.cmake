add_library(clara::clara INTERFACE IMPORTED)

set_target_properties(clara::clara
    PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${LIBS}/clara)
