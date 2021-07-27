# IGA-ADS

[![Build status](https://github.com/marcinlos/iga-ads/actions/workflows/build.yml/badge.svg)](https://github.com/marcinlos/iga-ads/actions?query=workflow%3ABuild+branch%3Adevelop)
[![Codecov](https://codecov.io/gh/marcinlos/iga-ads/branch/develop/graph/badge.svg?token=rkZgDGLoTy)](https://codecov.io/gh/marcinlos/iga-ads)
[![License](https://img.shields.io/github/license/marcinlos/iga-ads)](https://raw.githubusercontent.com/marcinlos/iga-ads/develop/LICENSE)
[![Language](https://img.shields.io/badge/C%2B%2B-17-b.svg?logo=cplusplus)](https://isocpp.org/)
[![Latest release](https://img.shields.io/github/v/release/marcinlos/iga-ads)](https://github.com/marcinlos/iga-ads/releases)

IGA-ADS is a C++ framework designed to facilitate creating parallel numerical simulations for
time-dependent and stationary PDEs using isogeometric finite element method.

## Building

### Dependencies

#### Tools
- Modern C++ compiler with C++17 support (GCC >= 7, Clang >= 6)
- [CMake](https://cmake.org/)
  (>= 3.13, 3.20 recommended for [presets](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html))

#### Libraries
- [BLAS](http://www.netlib.org/blas/) and [LAPACK](https://www.netlib.org/lapack/)
- [Boost](https://www.boost.org/) (>= 1.58)
- [{fmt}](https://github.com/fmtlib/fmt) (>= 7.1)

#### Optional
- [Galois](https://iss.oden.utexas.edu/?p=projects/galois) (>= 6.0) - parallel assembly
  (**recommended**)
- [MUMPS](http://mumps.enseeiht.fr/) - used for stationary problems
- [Catch2](https://github.com/catchorg/Catch2) - required to build test suite

### Compilation

To compile the code, in the project root directory execute the following commands:
```bash
mkdir build-dir
cmake -S. -B build-dir [OPTIONS]
cmake --build build-dir --parallel
```
where `OPTIONS` are additional build settings described [below](#build-options).
By default this builds the entire project, including examples and supporting applications.
To build only the library, add `--target ADS` to the last command:
```bash
cmake --build build-dir --parallel --target ADS
```

Once the library is compiled, in can be installed using
```bash
cmake --install build-dir
```
This installs the library, headers and CMake configuration files in a default system location.
To install in a different directory, specify it by using the `--prefix` option as follows:
```bash
cmake --install build-dir --prefix /path/to/desired/directory
```

For more details and additional options, consult
[CMake documentation](https://cmake.org/cmake/help/latest/manual/cmake.1.html).

#### Build options

Options are specified as `-D option=value`, e.g.  `cmake -S . -B build -D ADS_USE_GALOIS=ON`
- `ADS_USE_GALOIS` - decides if parallel executor using Galois framework is compiled. Disabling this
  options also stops example programs using it from being compiled (default: `OFF`)
- `ADS_USE_MUMPS` - decides if MUMPS support is included (default: `OFF`)
- `ADS_BUILD_PROBLEMS` - decides if the example problems are compiled (default: `ON`)
- `ADS_BUILD_TOOLS` - decides if the supporting applications are compiled (default: `ON`)
- [`BUILD_SHARED_LIBS`](https://cmake.org/cmake/help/latest/variable/BUILD_SHARED_LIBS.html) -
  decides if the project is built as a shared library (default: `OFF`)
- [`BUILD_TESTING`](https://cmake.org/cmake/help/latest/module/CTest.html) - decides if the tests
  are compiled (default: `OFF`)

## Using the library

There are three primary ways to use IGA-ADS library in your CMake project.

#### Importing installed package

If IGA-ADS has been built and installed as described [above](#compilation),
it can be imported using `find_package` command.
To import and use the library with some example application,
add the following to you `CMakeLists.txt`:
```cmake
find_package(ADS 0.1.0 REQUIRED)
add_executable(example ...)
target_link_libraries(example PRIVATE ADS::ADS)
```
Note that if the library has been installed in a non-standard location,
it may be necessary to inform CMake about it via
[`CMAKE_PREFIX_PATH`](https://cmake.org/cmake/help/latest/variable/CMAKE_PREFIX_PATH.html)
or
[`ADS_ROOT`](https://cmake.org/cmake/help/latest/variable/PackageName_ROOT.html)
option:
```bash
cmake -D ADS_ROOT=/path/to/install/dir ...
```

#### Including as a subdirectory

In this method, the entire IGA-ADS directory is added as a subdirectory to the project using it.
One easy way to do it is by using git
[submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules):
```bash
git submodule add https://github.com/marcinlos/iga-ads ads
```
The `ads` directory then needs to be added in `CMakeLists.txt` using `add_subdirectory`.
We can also set [build options](#build-options) for IGA-ADS:
```cmake
set(ADS_USE_GALOIS ON)
set(ADS_USE_MUMPS ON)
set(ADS_BUILD_PROBLEMS OFF)
set(ADS_BUILD_TOOLS OFF)

add_subdirectory(ads)
add_executable(example ...)
target_link_libraries(example PRIVATE ADS::ADS)
```

#### Using `FetchContent`

Using [`FetchContent`](https://cmake.org/cmake/help/latest/module/FetchContent.html)
we can automatically download IGA-ADS to the build directory.
```cmake
include(FetchContent)

FetchContent_Declare(ADS
  GIT_REPOSITORY https://github.com/marcinlos/iga-ads
  GIT_TAG develop
)

set(ADS_USE_GALOIS ON)
set(ADS_USE_MUMPS ON)
set(ADS_BUILD_PROBLEMS OFF)
set(ADS_BUILD_TOOLS OFF)

FetchContent_MakeAvailable(ADS)
add_executable(example ...)
target_link_libraries(example PRIVATE ADS::ADS)
```
**Note:** `FetchContent_MakeAvailable` requires CMake >= 3.14, see
[here](https://cmake.org/cmake/help/latest/module/FetchContent.html#fetch-content-canonical-pattern)
for solution working in earlier versions.

## Citation

If you use this code, please cite:
- Marcin Łoś, Maciej Woźniak, Maciej Paszyński, Andrew Lenharth and Keshav Pingali.
  IGA-ADS : Isogeometric Analysis FEM using ADS solver.
  _Computer & Physics Communications_, 217:99-116, 2017.
- Marcin Łoś,  Maciej Paszyński, Adriank Kłusek and Witold Dzwinel.
  Application of fast isogeometric L2 projection solver for tumor growth simulations.
  _Computer Methods in Applied Mechanics and Engineering_, 316:1257-1269, 2017.
- Marcin Łoś,  Adriank Kłusek, M. Amber Hassaan, Keshav Pingali, Witold Dzwinel and Maciej Paszyński.
  Parallel fast isogeometric L2 projection solver with GALOIS system for 3D tumor growth simulations.
  _Computer Methods in Applied Mechanics and Engineering_, 343:1-22, 2019.
