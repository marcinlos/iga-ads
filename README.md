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

## Compilation

To compile the code, create a directory for the build and execute following commands:

$ cmake <options> ${PROJECT_DIR}
$ make

where ${PROJECT_DIR} is the root directory of the project source. Options allow customizing which
parts of the project are compiled. By default parallelization using Galois is disabled, example
applications are compiled and tests are skipped.

- USE_GALOIS - whether or not parallel executor using Galois framework is compiled. Disabling this options also stops example programs using it from being compiled (default: OFF)
- COMPILE_TESTS - whether the unit tests are compiled. This should be disabled if libunittest is not available (default: OFF)
- SKIP_PROBLEMS - if ON, the example problems are not compiled (default: OFF)

Options are specified as -Doption=value, e.g.
$ cmake -DUSE_GALOIS=ON ..

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
