# Custom target that depends on all tools
add_custom_target(ads.tools)

if(ADS_USE_GALOIS)
    add_executable(ads.tools.error error.cpp)
    target_link_libraries(
        ads.tools.error
        PRIVATE #
            ADS::ADS
            ads-options-private
    )
    set_target_properties(
        ads.tools.error
        PROPERTIES #
            OUTPUT_NAME ads-error
    )
    add_dependencies(ads.tools ads.tools.error)
endif()

# Installation rules for tools
get_target_property(ADS_TOOLS_TARGETS ads.tools MANUALLY_ADDED_DEPENDENCIES)

if(ADS_TOOLS_TARGETS)
    install(TARGETS ${ADS_TOOLS_TARGETS} COMPONENT ads-dev)
endif()
