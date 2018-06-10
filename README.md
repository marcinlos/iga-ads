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

- USE_GALOIS - wether or not parallel executor using Galois framework is compiled. Disabling this options also stops example programs using it from being compiled (default: OFF)
- COMPILE_TESTS - wether the unit tests are compiled. This should be disabled if libunittest is not available (default: OFF)
- SKIP_PROBLEMS - if ON, the example problems are not compiled (default: OFF)

Options are specified as -Doption=value, e.g. 
$ cmake -DUSE_GALOIS=ON ..


## Contents

```
 Top-level structure:
  src/ads/      - framework code
  src/problems  - example problem implementations
  test/         - unit tests

Detailed description of the package contents:

CMake build configuration:
  CMakeLists.txt

Simulation infrastructure:
  src/ads/simulation.hpp
  src/ads/simulation/config.hpp
  src/ads/simulation/dimension.cpp
  src/ads/simulation/dimension.hpp
  src/ads/simulation/simulation_1d.cpp
  src/ads/simulation/simulation_1d.hpp
  src/ads/simulation/simulation_2d.cpp
  src/ads/simulation/simulation_2d.hpp
  src/ads/simulation/simulation_3d.cpp
  src/ads/simulation/simulation_3d.hpp
  src/ads/simulation/simulation_base.cpp
  src/ads/simulation/simulation_base.hpp
  src/ads/basis_data.cpp
  src/ads/basis_data.hpp

ADS solver implementation:
  src/ads/solver.hpp

B-spline functions implementation:
  src/ads/bspline/bspline.cpp
  src/ads/bspline/bspline.hpp
  src/ads/bspline/eval.hpp

Executors (single/multithreaded):
  src/ads/executor/galois.cpp
  src/ads/executor/galois.hpp
  src/ads/executor/sequential.hpp

Quadratures:
  src/ads/quad/gauss.hpp
  src/ads/quad/gauss_data.cpp

Computing B-spline basis Gram matrix:
  src/ads/mass_matrix.cpp
  src/ads/mass_matrix.hpp

Computing L2-projection:
  src/ads/projection.hpp

Multi-dimensional arrays:
  src/ads/util/multi_array.hpp
  src/ads/util/multi_array/base.hpp
  src/ads/util/multi_array/ordering/reverse.hpp
  src/ads/util/multi_array/ordering/standard.hpp
  src/ads/util/multi_array/reshape.hpp
  src/ads/util/multi_array/wrapper.hpp

Tensor storage and basic operations:
  src/ads/lin/band_matrix.hpp
  src/ads/lin/band_solve.hpp
  src/ads/lin/tensor.hpp
  src/ads/lin/tensor/base.hpp
  src/ads/lin/tensor/cyclic_transpose.hpp
  src/ads/lin/tensor/defs.hpp
  src/ads/lin/tensor/equality.hpp
  src/ads/lin/tensor/for_each.hpp
  src/ads/lin/tensor/io.hpp
  src/ads/lin/tensor/reshape.hpp
  src/ads/lin/tensor/tensor.hpp
  src/ads/lin/tensor/view.hpp

File output:
  src/ads/output/axis.hpp
  src/ads/output/gnuplot.hpp
  src/ads/output/grid.hpp
  src/ads/output/output_base.hpp
  src/ads/output/output_format.hpp
  src/ads/output/output_manager_base.hpp
  src/ads/output/range.hpp
  src/ads/output/raw.hpp
  src/ads/output/vtk.hpp
  src/ads/output_manager.hpp

Generic utility functions:
  src/ads/util.hpp
  src/ads/util/function_value.hpp
  src/ads/util/function_value/function_value_1d.hpp
  src/ads/util/function_value/function_value_2d.hpp
  src/ads/util/function_value/function_value_3d.hpp
  src/ads/util/io.hpp
  src/ads/util/iter/product.hpp
  src/ads/util/math/vec.hpp
  src/ads/util/math/vec/functions.hpp
  src/ads/util/math/vec/operators.hpp
  src/ads/util/math/vec/vec_2d.hpp
  src/ads/util/math/vec/vec_3d.hpp
  src/ads/util/math/vec/vec_fwd.hpp
  src/ads/util/meta.hpp

Problem implementations:
  src/problems/elasticity/elasticity.hpp
  src/problems/elasticity/main.cpp
  src/problems/flow/environment.hpp
  src/problems/flow/flow.hpp
  src/problems/flow/geometry.hpp
  src/problems/flow/main.cpp
  src/problems/flow/pumps.hpp
  src/problems/heat/heat_1d.hpp
  src/problems/heat/heat_2d.cpp
  src/problems/heat/heat_2d.hpp
  src/problems/heat/heat_3d.cpp
  src/problems/heat/heat_3d.hpp
  src/problems/validation/main.cpp
  src/problems/validation/validation.hpp

Unit tests:
  test/ads/bspline/bspline_test.cpp
  test/ads/bspline/eval_test.cpp
  test/ads/lin/banded_solver_test.cpp
  test/ads/lin/tensor_test.cpp
  test/ads/util/multi_array_test.cpp
```

## Compiling Lapack & BLAS

Both libraries are very hardware-specific. Make sure that if your cluster is not homogenous you use the proper version of them at each node.
If you need to compile them use the instructions from [here](https://stackoverflow.com/questions/23463240/how-to-compile-lapack-so-that-it-can-be-used-correctly-during-installation-of-oc).