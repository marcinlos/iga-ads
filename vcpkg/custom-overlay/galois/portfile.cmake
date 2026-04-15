# vcpkg_check_linkage(ONLY_STATIC_LIBRARY)

vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO IntelligentSoftwareSystems/Galois
    REF "release-${VERSION}"
    SHA512 634982846bc219f56d21d40f4099bf6f48cfeed95c2b23bc49e5ca396b4eb6efdb9ccd1d1b43d52d33f7d01ba599d84262ad79273ee6513985a5a00b7be44b29
    PATCHES #
        add-missing-includes.patch
        fix-cpuinfo-parsing.patch
        restrict-targets.patch
        build-without-llvm.patch
        skip-tests.patch
        remove-legacy-vars.patch
        find-boost-by-config.patch
        no-need-for-boost-iostreams.patch
)

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
)

vcpkg_cmake_build(TARGET galois_shmem)
vcpkg_cmake_install()

vcpkg_cmake_config_fixup(PACKAGE_NAME "galois" CONFIG_PATH "lib/cmake/Galois")

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")

file(
    INSTALL "${CMAKE_CURRENT_LIST_DIR}/usage"
    DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}"
)
vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE.txt")
