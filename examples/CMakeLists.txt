# Custom target that depends on all example applications
add_custom_target(ads.examples)

# Used by many examples to parse command arguments
find_package(lyra REQUIRED)

# Problem definitions
add_subdirectory(heat)
add_subdirectory(poisson)
add_subdirectory(validation)
add_subdirectory(scalability)
add_subdirectory(flow)
add_subdirectory(elasticity)
add_subdirectory(stokes)
add_subdirectory(pollution)
add_subdirectory(victor)
add_subdirectory(erikkson)
add_subdirectory(multistep)
add_subdirectory(demkowicz)
add_subdirectory(cg)
add_subdirectory(maxwell)
add_subdirectory(tumor)
add_subdirectory(implicit)
add_subdirectory(dg)
add_subdirectory(co2_sequestration)
add_subdirectory(fire)
