list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules/")

# Imaginary libraries to propagate settings
add_library(ads-options-public INTERFACE)
add_library(ads-options-private INTERFACE)
target_compile_features(ads-options-public INTERFACE cxx_std_17)
