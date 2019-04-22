Getting started with the GEM Library
====================================

Here is a short introduction to the usage of GEM library in MATLAB and GNU Octave. This introduction assumes that the gem folder is in MATLAB/GNU Octave's path and that it is working (this can be checked by typing e.g. `gem('pi')` into MATLAB/GNU Octave; if this command yields an error, check that the GEM library is compiled for your system and that its folder is in MATLAB/GNU Octave's path).

## Creating High Precision Matrices

### Dense matrices

A matrix for high precision computation can be created from any MATLAB/GNU Octave matrix by calling the `gem` constructor.
 - `a = gem([1 2; 3 4+5i])` is a 2x2 gem object

    This matrix can be manipulated like a conventional matrix. For instance, `a(2,2)` extracts the (2,2) element of the matrix: `4+5i`.

When constructing a gem object from a double-precision matrix, the first 15 digits of each component are taken into account and completed by trailing zeros (up to the precision of the gem object). This garantees that all numbers are transfered with a predictable way, and that integers remain integers after conversion.

High precision objects can also be constructed by specifying all digits in a string. This allows to define numbers with more than 15 digits precision, as in the following example:
 - `gem('1.93754981629874513245725542343')`

    creates a 1x1 matrix (i.e. a scalar) with 20 non-trivial digits. The remaining digits are set to zeros.

When creating a number from a string which contains more digits than the precision defined by `gemWorkingPrecision`, the precision of the created number is automatically adjusted to hold all specified digits.

All elements of a matrix can be individually initialized to an arbitrary precision by calling the constructor `gem` with a cell array containing the corresponding strings:
 - `gem({'1.321321123123321123' '456.4566544566544564'; '0.789987987789987789' '369.639369366963936'})`


### Sparse matrices

Just like in the case of dense matrices, sparse MATLAB/GNU Octave matrix are easily transfered to sparse gem object by calling the `sgem` constructor:
 - `sgem(speye(3))` is a sparse 3x3 identity matrix

More generally, a high precision object can be constructed from a MATLAB/GNU Octave matrix, whether it is full or sparse, by calling the `sparsify` function. This will create a gem object when given a full matrix, and an sgem object when given a sparse one.

Sgem object can also be constructed by calling the `sparse(i,j,s,m,n)` function with a gem vector `s`:
 - `sparse([1 2], [1 3], gem({'123.45678901234567890' '987.653210987654321'}), 3, 3)` is a sparse matrix with two high precision components


### Special values

The following mathematical constants can be constructed explicitely
 - `gem('pi')` is 3.14159265...
 - `gem('e')` is 2.7182818...
 - `gem('log2')` is 0.69314718...
 - `gem('euler')` is 0.57721566...
 - `gem('catalan')` is 0.91596559...


## Precision

The GEM library works with high precision numbers. The precision of these numbers is defined by the number of digits that are used to describe them (in basis 10). By default, the library takes into account 50 digits. This means that it can distinguish between numbers that differ at the 50th decimal place:
 - In double precision `(pi/10) - (pi/10 - 1e-17)` yields

        0
    but `(gem('pi')/10) - (gem('pi')/10 - 1e-50)` is

        ~1e-50

The precision at which the GEM library works can be adjusted throught the function `gemWorkingPrecision`. Alternatively, it can be passed as a parameter when creating a gem object:
 - `gem(3.141592654,3)` yields just 

        3.14
    (the remaining digits are ignored)

When creating a gem object from a double, at most 15 digits are taken into account. When using strings of digits to specify a number, the precision is automatically adjusted to guarantee that all provided digits are taken into account, unless explicitely stated:
 - `gem('123456789012345678901234567890123456789012345678901234567891') - gem('123456789012345678901234567890123456789012345678901234567890')` produces the result

        1
    even though the relative difference is of the order of 1e-60. In this case, the actual precision of the numbers is higher than the current working precision.

 - `gem('1234567890',3) - gem('1234567891',3)` produces

        0
    because the three first digits of both numbers are equal.

Note that the accuracy of the result of an operation is not only determined by the number of decimals that the library takes into account to define numbers, but it may also depend on the [number of elementary manipulations](https://en.wikipedia.org/wiki/Numerical_error) that were used to produce the result.

The `precision` function can be used to determine how many digits are used to describe a particular number. For instance,
 - `precision(gem('3.2'))` is

        50
    `precision(gem('123456789012345678901234567890123456789012345678901234567891'))` is

        60
    while `precision(gem('e',5))` is

        5


### Display precision

When displaying a high precision number, only part of the number is usually printed out. This is the default behavior, because long string of digits can quickly become cumbersome. Even though all significant digits are not systematically printed, they are indeed kept in memory and used in the computations. The function `gemDisplayPrecision` can be used to adjust the number of displayed digits. Alternatively, the display precision can also be directly specified as a second argument of the `disp` or `display` functions:
 - `disp(gem('pi'),50)` prints 50 digits:

        3.1415926535897932384626433832795028841971693993751

A *negative* display precision prints out all digits in memory:
 - `disp(gem('pi'),-1)` prints the same 50 digits as above.




## Functions

A list of the functions that can be applied to gem and sgem objects is available [here](functions.md). In general, these functions take the same arguments as their MATLAB/GNU Octave counterparts, and behave in the same way. For instance, the `max` function applied to a vector containing two complex numbers returns the complex number with largest magnitude, or the one with largest angle if both magnitudes are equal.

For some functions, all possible behaviors are not yet implemented. For instance, `mpower` currently only supports powers of +1 and -1. Anyone wishing to have a specific feature is invited to open an [issue](https://github.com/jdbancal/gem/issues) and eventually consider contributing to this open source project by implementing the missing behaviors and submitting a push request (see also [this page](howToContribute.md) for more ways to get involved).

Note that the default behavior of some GEM functions does differ from MATLAB's implementation. This is the case for example for functions which don't preserve the sparsity when applied to sparse objects. This point is explained in the `gemSparseLikeMatlab.m` file. This default MATLAB behavior is easily restored by calling `gemSparseLikeMatlab(1)`.

