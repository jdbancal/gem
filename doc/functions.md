Included functions
==================

Here is a list of the functions currently implemented in the GEM library. The first table shows functions which can be called without reference to a gem/sgem object (this can be used to create a gem instance or to adjust some library parameters for instance).


| function | full support | sparse support | remarks |
|----------|---|---|---|
| gem | ✔ | ✔ | gem object constructor |
| sgem | ✔ | ✔ | sparse gem object constructor |
| gemify | ✔ | ✔ | converts a matrix to gem or sgem object, preserving is sparsity |
| gemRand | ✔ | ✗ | generates high precision random numbers |
| gemRng | ✔ | ✗ | gemRand random seed |
| gemWorkingPrecision | ✔ | ✔ | sets the working precision of the library |
| gemDisplayPrecision | ✔ | ✔ | sets the precision used when displaying gem and sgem objects |
| gemSparseLikeMatlab | ✗ | ✔ | sets whether functions not preserving zeros should be allowed to produce sparse matrices |
| make | ✔ | ✔ | compiles the library |


List of methods:
----------------

The following function can be applied to gem/sgem objects directly.

| methods | full support | sparse support | remarks |
|----------|---|---|---|
| abs  | ✔  | ✔  |   |
| acos | ✔ | ✔ |  |
| acot | ✔ | ✔ |  |
| acsc | ✔ | ✔ |  |
| and | ✔ | ✔ |  |
| angle | ✔ | ✔ |  |
| asec | ✔ | ✔ |  |
| asin | ✔ | ✔ |  |
| atan | ✔ | ✔ |  |
| cat | ✔ | ✔ |  |
| cbrt | ✔ | ✔ |  |
| clean | ✗ | ✔ |  |
| ceil | ✔ | ✔ |  |
| colon | ✔ | ✔ |  |
| complex | ✔ | ✔ |  |
| conj | ✔ | ✔ |  |
| cos | ✔ | ✔ |  |
| cot | ✔ | ✔ |  |
| csc | ✔ | ✔ |  |
| ctranspose | ✔ | ✔ |  |
| diag | ✔ | ✔ |  |
| disp | ✔ | ✔ |  |
| display | ✔ | ✔ |  |
| double | ✔ | ✔ |  |
| eig | ✔ | ✗ |  |
| eigs | ✔ | ✔ |  |
| end | ✔ | ✔ |  |
| eq | ✔ | ✔ |  |
| exp | ✔ | ✔ |  |
| find | ✔ | ✔ |  |
| fix | ✔ | ✔ |  |
| floor | ✔ | ✔ |  |
| full | ✔ | ✔ |  |
| ge | ✔ | ✔ |  |
| gt | ✔ | ✔ |  |
| horzcat | ✔ | ✔ |  |
| imag | ✔ | ✔ |  |
| inv | ✔ | ✔ |  |
| isempty | ✔ | ✔ |  |
| isequal | ✔ | ✔ |  |
| isequaln | ✔ | ✔ |  |
| isfinite | ✔ | ✔ |  |
| ishermitian | ✔ | ✔ |  |
| isinf | ✔ | ✔ |  |
| isnan | ✔ | ✔ |  |
| isnumeric | ✔ | ✔ |  |
| isreal | ✔ | ✔ |  |
| issparse | ✔ | ✔ |  |
| issymmetric | ✔ | ✔ |  |
| kron | ✔ | ✔ | supports more than 2 inputs, as in kron({A,B,C}) |
| ldivide | ✔ | ✔ |  |
| le | ✔ | ✔ |  |
| length | ✔ | ✔ |  |
| loadobj | ✔ | ✔ |  |
| log | ✔ | ✔ |  |
| log10 | ✔ | ✔ |  |
| log2 | ✔ | ✔ |  |
| logical | ✔ | ✔ |  |
| lt | ✔ | ✔ |  |
| max | ✔ | ✔ |  |
| min | ✔ | ✔ |  |
| minus | ✔ | ✔ |  |
| mldivide | ✔ | ✔ |  |
| mpower | ✔ |  | supports two scalars or matrices with +/-1 exponent |
| mrdivide | ✔ | ✔ |  |
| mtimes | ✔ | ✔ |  |
| ne | ✔ | ✔ |  |
| norm | ✔ | ✔ |  |
| not | ✔ | ✔ |  |
| numel | ✔ | ✔ |  |
| or | ✔ | ✔ |  |
| plot | ✔ | ✔ |  |
| plus | ✔ | ✔ |  |
| power | ✔ | ✔ |  |
| precision | ✔ | ✔ |  |
| prod | ✔ | ✔ |  |
| rdivide | ✔ | ✔ |  |
| real | ✔ | ✔ |  |
| reshape | ✔ | ✔ |  |
| round | ✔ | ✔ |  |
| saveobj | ✔ | ✔ |  |
| sec | ✔ | ✔ |  |
| sin | ✔ | ✔ |  |
| size | ✔ | ✔ |  |
| sparse | ✔ | ✔ |  |
| sqrt | ✔ | ✔ |  |
| subsasgn | ✔ | ✔ |  |
| subsref | ✔ | ✔ |  |
| sum | ✔ | ✔ |  |
| svd | ✔ | ✗ |  |
| svds | ✔ | ✔ |  |
| tan | ✔ | ✔ |  |
| times | ✔ | ✔ |  |
| toStrings | ✔ | ✔ |  |
| transpose | ✔ | ✔ |  |
| uminus | ✔ | ✔ |  |
| uplus | ✔ | ✔ |  |
| vertcat | ✔ | ✔ |  |

