include(GNUInstallDirs)

# Define an export set for the main library
install(
    TARGETS ADS ads-objects ads-options-public ads-options-private
    EXPORT ads-targets
    LIBRARY COMPONENT ads-rt
    ARCHIVE COMPONENT ads-dev
    FILE_SET HEADERS COMPONENT ads-dev
)

# Install targets from the export set (library and target definitions)
install(
    EXPORT ads-targets
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/ads"
    NAMESPACE ADS::
    COMPONENT ads-dev
)

include(CMakePackageConfigHelpers)

configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/ads-config.cmake.in"
    ads-config.cmake
    @ONLY
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/ads-version.cmake"
    COMPATIBILITY ExactVersion
)

# Copy package configuration files
install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/ads-config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/ads-version.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/ads"
    COMPONENT ads-dev
)

# Copy find modules defined in the project
install(
    DIRECTORY cmake/Modules
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/ads"
    COMPONENT ads-dev
    FILES_MATCHING
    PATTERN "Find*"
)
