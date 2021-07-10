target_sources(ADS
  PRIVATE
    ads/basis_data.cpp
    ads/form_matrix.cpp
    ads/bspline/bspline.cpp
    ads/executor/galois.cpp
    ads/quad/gauss_data.cpp
    ads/simulation/dimension.cpp
    ads/simulation/simulation_base.cpp
    ads/simulation/simulation_1d.cpp
    ads/simulation/simulation_2d.cpp
    ads/simulation/simulation_3d.cpp
)

# Add configured source file with version information
target_sources(ADS PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/ads/version.cpp")
