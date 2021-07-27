# IGA-ADS

IGA-ADS is a C++ framework designed to facilitate creating parallel numerical simulations for time-dependent PDEs using isogeometric finite element method.


## Requirements

1. Dependencies
- LAPACK, BLAS
- Boost, version 1.58 or higher (http://www.boost.org/)
- (optional) Galois framework, version 2.2.1 (http://iss.ices.utexas.edu/?p=projects/galois)
- (optional) libunittest (http://libunittest.sourceforge.net/)

Galois is required for parallelism, libunittest for unit tests only.

2. Tools
- compiler: reasonable C++14 support is required (framework has been tested with GCC 5.3.1)
- build system: CMake 3.1 or higher


2. Compilation

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
