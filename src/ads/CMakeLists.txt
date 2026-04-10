configure_file(version.cpp.in version.cpp @ONLY)

target_sources(
    ads-objects
    PRIVATE
        basis_data.cpp
        form_matrix.cpp
        bspline/bspline.cpp
        executor/galois.cpp
        quad/gauss_data.cpp
        simulation/dimension.cpp
        simulation/simulation_base.cpp
        simulation/simulation_1d.cpp
        simulation/simulation_2d.cpp
        simulation/simulation_3d.cpp
        "${CMAKE_CURRENT_BINARY_DIR}/version.cpp"
)
