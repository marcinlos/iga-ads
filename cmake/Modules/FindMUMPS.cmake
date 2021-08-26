include(FindPackageHandleStandardArgs)

# Find the C dobule header file
find_path(MUMPS_INCLUDE_DIR NAMES dmumps_c.h)

find_library(MUMPS_LIBRARY
  NAMES dmumps dmumps_seq
)

find_package_handle_standard_args(MUMPS
  FOUND_VAR MUMPS_FOUND
  REQUIRED_VARS
    MUMPS_LIBRARY
    MUMPS_INCLUDE_DIR
)

set(MUMPS_LIBRARIES ${MUMPS_LIBRARY})

# mumps_common is a separate library
find_library(MUMPS_COMMON_LIBRARY NAMES mumps_common)

if (MUMPS_COMMON_LIBRARY)
  list(APPEND MUMPS_LIBRARIES ${MUMPS_COMMON_LIBRARY})
endif()

# Find scalapack
find_library(MUMPS_SCALAPACK NAMES scalapack scalapack-openmpi)

if (MUMPS_SCALAPACK)
  list(APPEND MUMPS_LIBRARIES ${MUMPS_SCALAPACK})
endif()

# Find BLAS and LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

list(APPEND MUMPS_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

# Try to find ordering libraries
set(ORDERING_LIBS
  pord
  metis scotch scotcherr esmumps
  parmetis ptscotch ptscotcherr ptesmumps
)

foreach (lib IN LISTS ORDERING_LIBS)
  find_library(MUMPS_${lib}_LIBRARY NAMES ${lib})

  if (MUMPS_${lib}_LIBRARY)
    list(APPEND MUMPS_LIBRARIES ${MUMPS_${lib}_LIBRARY})
  endif()
endforeach()

# We need to find some Fortran components
enable_language(Fortran)

# Link with Fortran standard library
list(APPEND MUMPS_LIBRARIES ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})

# Find MPI
find_package(MPI REQUIRED COMPONENTS CXX Fortran)
list(APPEND MUMPS_LIBRARIES MPI::MPI_CXX MPI::MPI_Fortran)

# Try to enable OMP support to leverage MUMPS parallelism
find_package(OpenMP COMPONENTS Fortran)

if (OpenMP_FOUND)
  list(APPEND MUMPS_LIBRARIES OpenMP::OpenMP_Fortran)
endif()

# Create imported targets
if (MUMPS_FOUND AND NOT TARGET MUMPS::MUMPS)
  add_library(MUMPS::MUMPS INTERFACE IMPORTED)
  set_target_properties(MUMPS::MUMPS PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${MUMPS_LIBRARIES}")
endif()
