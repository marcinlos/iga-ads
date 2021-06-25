find_package(Git REQUIRED)

# Read last commit SHA
execute_process(
  COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
  RESULT_VARIABLE result
  OUTPUT_VARIABLE ADS_GIT_COMMIT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if (result)
  message(FATAL_ERROR "Failed to get git commit SHA: ${result}")
endif()

configure_file(src/ads/version.cpp.in ads/version.cpp @ONLY)
