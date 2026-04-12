include_guard(GLOBAL)

if(ADS_DEV_SANITIZERS)
    string(
        JOIN " "
        SANITIZER_FLAGS
        -fsanitize=address
        -fsanitize=leak
        -fsanitize=pointer-compare
        -fsanitize=pointer-subtract
        -fsanitize=undefined
        -fsanitize-undefined-trap-on-error
        -fno-omit-frame-pointer
    )
endif()

string(
    JOIN " "
    WARNING_FLAGS
    -Wall
    -Wextra
    -pedantic
    # -Wshadow
    -Werror=return-type
    -Wsuggest-override
    -Wold-style-cast
    -Wcast-align
    -Wconversion
    # -Wsign-conversion
    -Wmissing-declarations
    -Wredundant-decls
    -Wmisleading-indentation
    -Wextra-semi
)
