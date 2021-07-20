if (NOT ADS_ENABLE_COVERAGE)
  return()
endif()

find_program(ADS_LCOV lcov)
find_program(ADS_GENHTML genhtml)
find_program(ADS_LLVM_COV llvm-cov)
find_program(ADS_LLVM_PROFDATA llvm-profdata)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

  set(ADS_CXX_FLAGS_COVERAGE --coverage -fno-inline -O0)
  set(ADS_EXE_LINKER_FLAGS_COVERAGE --coverage)

  add_custom_target(coverage
    COMMAND
      ${ADS_LCOV}
        --zerocounters
        --directory ${PROJECT_BINARY_DIR}
        --rc lcov_branch_coverage=1
    COMMAND
      ${ADS_LCOV}
        --capture --initial
        --directory ${PROJECT_BINARY_DIR}
        --output-file baseline.info
        --rc lcov_branch_coverage=1
    COMMAND
      ${CMAKE_CTEST_COMMAND}
    COMMAND
      ${ADS_LCOV}
        --capture
        --directory ${PROJECT_BINARY_DIR}
        --output-file tests.info
        --rc lcov_branch_coverage=1
    COMMAND
      ${ADS_LCOV}
        --add-tracefile baseline.info
        --add-tracefile tests.info
        --output-file coverage-full.info
        --rc lcov_branch_coverage=1
    COMMAND
      ${ADS_LCOV}
        --extract coverage-full.info
          '${PROJECT_SOURCE_DIR}/src/*'
          '${PROJECT_SOURCE_DIR}/include/*'
        --output-file coverage.info
        --rc lcov_branch_coverage=1
    COMMAND
      ${ADS_GENHTML} coverage.info
        --prefix ${PROJECT_SOURCE_DIR}
        --output-directory html
        --show-details
        --title "ADS coverage"
        --demangle-cpp
        --legend
        --branch-coverage
    DEPENDS ads-suite
  )

elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")

  set(ADS_CXX_FLAGS_COVERAGE
    -fprofile-instr-generate
    -fcoverage-mapping
    -O0
    -fno-inline
    -fno-elide-constructors
  )

  set(ADS_EXE_LINKER_FLAGS_COVERAGE
    -fprofile-instr-generate
  )

  add_custom_target(coverage
    COMMAND
      ${CMAKE_COMMAND} -E env
      LLVM_PROFILE_FILE="ads-%p.profraw" ${CMAKE_CTEST_COMMAND}
    COMMAND
      find ${PROJECT_BINARY_DIR} -name "ads-*.profraw" > profraw-files
    COMMAND
      ${ADS_LLVM_PROFDATA}  merge -sparse
        -input-files profraw-files
        -o coverage.profdata
    COMMAND
      ${ADS_LLVM_COV} show
        -instr-profile coverage.profdata
        -Xdemangler=c++filt
        $<TARGET_FILE:ads-suite>
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        > coverage.txt
    COMMAND
      ${ADS_LLVM_COV} show
        -format=html
        -output-dir=html
        -Xdemangler=c++filt
        -instr-profile coverage.profdata
        $<TARGET_FILE:ads-suite>
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
    DEPENDS ads-suite
  )

endif()

target_compile_options(ads-options-private INTERFACE ${ADS_CXX_FLAGS_COVERAGE})
target_link_options(ads-options-private INTERFACE ${ADS_EXE_LINKER_FLAGS_COVERAGE})
