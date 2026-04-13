include_guard(GLOBAL)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Boost 1.70 CONFIG REQUIRED)
find_package(fmt 7.1 CONFIG REQUIRED)

target_link_libraries(
    ads
    PUBLIC #
        LAPACK::LAPACK
        BLAS::BLAS
        Boost::boost
        fmt::fmt
)

if(ADS_USE_MUMPS)
    find_package(MUMPS REQUIRED)
    target_link_libraries(ads PUBLIC MUMPS::MUMPS)
endif()

if(ADS_USE_GALOIS)
    find_package(Galois 6.0 CONFIG REQUIRED)
    target_link_libraries(ads PUBLIC Galois::shmem)
endif()
