include(FindPackageHandleStandardArgs)


find_path(MUMPS_INCLUDE_DIR
  NAMES dmumps_c.h
  PATH_SUFFIXES include mumps)

find_library(MUMPS_LIBRARY
  NAMES dmumps dmumps_seq
  PATH_SUFFIXES lib lib64)

find_package_handle_standard_args(MUMPS
  FOUND_VAR MUMPS_FOUND
  REQUIRED_VARS
    MUMPS_LIBRARY
    MUMPS_INCLUDE_DIR
)

set(MUMPS_LIBRARIES ${MUMPS_LIBRARY})

find_library(MUMPS_COMMON_LIBRARY
  NAMES mumps_common
  PATH_SUFFIXES lib lib64)

if (MUMPS_COMMON_LIBRARY)
  list(APPEND MUMPS_LIBRARIES ${MUMPS_COMMON_LIBRARY})
endif()

# Find scalapack
find_library(MUMPS_SCALAPACK
  NAMES scalapack
  PATH_SUFFIXES lib lib64)

if (MUMPS_SCALAPACK)
  list(APPEND MUMPS_LIBRARIES ${MUMPS_SCALAPACK})
endif()

# Try to find ordering libraries
set(ORDERING_LIBS metis parmetis ptesmumps scotch ptscotch ptscotcherr pord)

foreach (lib IN LISTS ORDERING_LIBS)
  find_library(MUMPS_${lib}_LIBRARY
    NAMES ${lib}
    PATH_SUFFIXES lib lib64)

  if (MUMPS_${lib}_LIBRARY)
    list(APPEND MUMPS_LIBRARIES ${MUMPS_${lib}_LIBRARY})
  endif()
endforeach()


list(APPEND MUMPS_LIBRARIES gfortran)

enable_language(Fortran)

# Find MPI
find_package(MPI REQUIRED)
list(APPEND MUMPS_LIBRARIES MPI::MPI_CXX MPI::MPI_Fortran)

# Try to enable OMP support to leverage MUMPS parallelism
find_package(OpenMP)

if (OpenMP_FOUND)
  list(APPEND MUMPS_LIBRARIES OpenMP::OpenMP_Fortran)
endif()

# Create imported targets
if (MUMPS_FOUND AND NOT TARGET MUMPS::MUMPS)
  add_library(MUMPS::MUMPS UNKNOWN IMPORTED)
  set_target_properties(MUMPS::MUMPS PROPERTIES
    IMPORTED_LOCATION "${MUMPS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${MUMPS_LIBRARIES}")
endif()
