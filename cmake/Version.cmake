find_package(Git)

if (NOT Git_FOUND)
  message("Git not available, no commit SHA")
else()
  # Read last commit SHA
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
    RESULT_VARIABLE result
    OUTPUT_VARIABLE ADS_GIT_COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  if (result)
    # If the command failed, the most likely reason is that
    # the project is not in Git repository (e.g. tarball release)
    message("Failed to get git commit SHA: ${result}")
  endif()
endif()

configure_file(src/ads/version.cpp.in ads/version.cpp @ONLY)
