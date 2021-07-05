# Possibly modified below
set(ADS_VERSION_FULL "${ADS_VERSION}")

find_package(Git)

if (Git_FOUND)
  # Read last commit SHA
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
    RESULT_VARIABLE exit_code
    OUTPUT_VARIABLE ADS_GIT_COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )

  if (exit_code)
    # If the command failed, the most likely reason is that
    # the project is not in Git repository (e.g. tarball release)
    message(STATUS "Failed to get git commit SHA: ${exit_code}")
  else()
    # Is HEAD exactly on a version tag?
    execute_process(
      COMMAND ${GIT_EXECUTABLE} describe --exact-match --match "v[0-9]*" HEAD
      RESULT_VARIABLE inexact_version
      OUTPUT_QUIET ERROR_QUIET
    )

    # If not, append GIT commit SHA to full version string
    if (inexact_version)
      set(ADS_VERSION_FULL "${ADS_VERSION}+${ADS_GIT_COMMIT}")
    endif()

  endif()
else()
  message(STATUS "Git not available, no commit SHA")
endif()

configure_file(src/ads/version.cpp.in ads/version.cpp @ONLY)
