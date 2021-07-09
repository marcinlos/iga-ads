# Possibly modified below
set(ADS_VERSION_FULL "${ADS_VERSION}")

find_package(Git)

if (Git_FOUND)
  # Determine if we are in the (correct) Git repository
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --show-toplevel
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    RESULT_VARIABLE exit_code
    OUTPUT_VARIABLE ADS_GIT_WORKTREE_ROOT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )

  if (exit_code)
    message(STATUS "No git repository found")
  elseif(NOT ADS_GIT_WORKTREE_ROOT STREQUAL PROJECT_SOURCE_DIR)
    message(STATUS "ADS git repository unavailable")
  else()
    # Read last commit SHA
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
      WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
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
  endif()
else()
  message(STATUS "Git not available, no commit SHA")
endif()

configure_file(src/ads/version.cpp.in src/ads/version.cpp @ONLY)
