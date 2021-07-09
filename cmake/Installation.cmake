
# Define an export set for the main library
install(
  TARGETS ADS ads-options-public ads-options-private
  EXPORT ads-targets
  LIBRARY
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT ads-runtime
  ARCHIVE
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT ads-devel
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install targets from the export set (library and target definitions)
install(
  EXPORT ads-targets
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/ads"
  NAMESPACE ADS::
  COMPONENT ads-devel
)

include(CMakePackageConfigHelpers)

configure_file(cmake/ads-config.cmake.in ads-config.cmake @ONLY)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/ads-config-version.cmake"
  COMPATIBILITY SameMajorVersion
)

# Copy package configuration files
install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/ads-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/ads-config-version.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/ads"
  COMPONENT ads-devel
)

# Copy find modules defined in the project
install(
  DIRECTORY cmake/Modules
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/ads"
  COMPONENT ads-devel
)

# Copy public header files
install(
  DIRECTORY include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT ads-devel
)
