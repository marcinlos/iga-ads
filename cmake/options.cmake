include_guard(GLOBAL)

option(
    ADS_USE_MUMPS
    "Use the MUMPS linear solver. Default: OFF. Values: { ON, OFF }."
    OFF
)

option(
    ADS_USE_GALOIS
    "Use the Galois framework. Default: OFF. Values: { ON, OFF }."
    OFF
)

option(
    ADS_BUILD_PROBLEMS
    "Build example problems. Default: ${PROJECT_IS_TOP_LEVEL}. Values: { ON, OFF }."
    ${PROJECT_IS_TOP_LEVEL}
)

option(
    ADS_BUILD_TOOLS
    "Build supporting tools. Default: ${PROJECT_IS_TOP_LEVEL}. Values: { ON, OFF }."
    ${PROJECT_IS_TOP_LEVEL}
)

option(
    ADS_BUILD_TESTS
    "Build tests. Default: ${PROJECT_IS_TOP_LEVEL}. Values: { ON, OFF }."
    ${PROJECT_IS_TOP_LEVEL}
)

option(
    ADS_DEV_SANITIZERS
    "Enable sanitizers. Default: OFF. Values: { ON, OFF }."
    OFF
)

option(
    ADS_DEV_COVERAGE
    "Gather code coverage data. Default: OFF. Values: { ON, OFF}."
    OFF
)
