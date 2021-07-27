# IGA-ADS

IGA-ADS is a C++ framework designed to facilitate creating parallel numerical simulations for time-dependent PDEs using isogeometric finite element method.

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
