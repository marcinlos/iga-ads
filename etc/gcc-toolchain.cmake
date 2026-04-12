include_guard(GLOBAL)

include("${CMAKE_CURRENT_LIST_DIR}/common.cmake")

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

string(
    JOIN " "
    WARNING_FLAGS #
    "${WARNING_FLAGS}"
    -Wduplicated-cond
    -Wduplicated-branches
    -Wuseless-cast
    -Wlogical-op
    -Wnoexcept
)

set(BASE_FLAGS "${WARNING_FLAGS} ${SANITIZER_FLAGS}")

set(DEBUG_FLAGS "-O0 -fno-inline -g3 ${BASE_FLAGS}")

set(CMAKE_C_FLAGS_DEBUG_INIT "${DEBUG_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG_INIT "${DEBUG_FLAGS}")

set(RELEASE_FLAGS "-Ofast -march=native ${BASE_FLAGS}")

set(CMAKE_C_FLAGS_RELEASE_INIT "${RELEASE_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE_INIT "${RELEASE_FLAGS}")
