include_guard(GLOBAL)

# Helper function for defining example applications
function(add_example name)
    set(optional_deps MUMPS GALOIS)
    set(single_value_args "")
    set(multi_value_args LIBS SRC)

    cmake_parse_arguments(
        arg
        "${optional_deps}"
        "${single_value_args}"
        "${multi_value_args}"
        ${ARGN}
    )

    if(arg_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "Bad arguments: ${arg_UNPARSED_ARGUMENTS}")
    endif()

    set(_target_name "ads.example.${name}")
    set(_define_target TRUE)

    foreach(_dep IN LISTS optional_deps)
        if(arg_${_dep} AND NOT ADS_USE_${_dep})
            set(_define_target FALSE)
        endif()
    endforeach()

    if(_define_target)
        add_executable(${_target_name} ${arg_SRC})
        target_link_libraries(
            ${_target_name}
            PRIVATE #
                ADS::ADS
                ads-options-private
                ${arg_LIBS}
        )
        set_target_properties(
            ${_target_name}
            PROPERTIES #
                OUTPUT_NAME ${name}
                RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin"
        )
        add_dependencies(ads.examples ${_target_name})
    endif()
endfunction()
