
set(COMMON_WARNINGS
  -Wall -Wextra -Wpedantic
  -Wnon-virtual-dtor
  -Woverloaded-virtual
  -Wmisleading-indentation
  -Wfloat-conversion
  -Wold-style-cast
  )

set(GCC_WARNINGS
  ${COMMON_WARNINGS}
  -Wsuggest-override
  -Wlogical-op
  -Wuseless-cast
  -Wduplicated-cond
  )

set(CLANG_WARNINGS
  ${COMMON_WARNINGS}
  -Wno-missing-braces
  )

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(ADS_WARNINGS ${GCC_WARNINGS})
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(ADS_WARNINGS ${CLANG_WARNINGS})
endif()

target_compile_options(ads-options-private INTERFACE ${ADS_WARNINGS})
