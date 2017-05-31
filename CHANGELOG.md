## [...] - ...
- Added sorting functions: sort, sortrow, issorted
- Added a small missing feature in eigs which required the sorting function
- Added diff
- Added support for logical indices in subsref and subsasgn
- Added modulo
- Added any
- Added support for algebra with non-numeric objects
- Improved display of super tiny numbers
- Fixed precision setting for sparse matrices
- Added support for partial precision in toStrings

## [0.1.3] - 2017-02-10
- Added isinteger, isfloat, and basic implementation of sprintf
- Added support for null extraction and assignments in subsref and subsasgn

## [0.1.2] - 2016-12-20
- Updated code to be compatible with latest Eigen, version 3.3
- Windows compatibility

## [0.1.1] - 2016-11-23
- Updated documentation
- Fixed issues when computing few eigenvalues of a small matrix
- Compiled with latest Spectra code
- Added changelog file

## [0.1.0] - 2016-11-19
- First release
- Implements basics matrix algebra for Sparse and dense matrices with arbitrary precision.
- Also implements simple functions like max/min, ... , some trigonometry, eigenvalues decompositions, singular values and basic linear systems solvers.
