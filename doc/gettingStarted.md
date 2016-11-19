Getting started with the GEM Library
====================================

Here is a short introduction to the usage of GEM library in matlab. This introduction assumes that the gem folder is in matlab's path and that it is working (this can be checked by typing e.g. `gem('pi')` into matlab; if this command yields an error, check that the gem library is compiled for your system and that its folder is in matlab's path).

## Precision

The GEM library can work with high precision numbers. This precision of a number is defined by the number of digits that are used to describe numbers (in basis 10). By default, the library takes into account 50 digits. This means that it can distinguish between numbers that differ at the 50th decimal place:
 - `(pi/10) - (pi/10 - 1e-17)` yields

        0
    but `(gem('pi')/10) - (gem('pi')/10 - 1e-50)` is

        ~1e-50

The precision at which the GEM library works can be adjusted throught the function `gemWorkingPrecision`. Alternatively, it can be passed as an parameter when a gem object is created:
 - `gem(3.141592654,3)` yields just 

        3.14
    (the remaining digits are ignored)

Note that when creating a gem object from a double, at most 15 digits are taken into account. When using strings of digits to specify a number, the precision is automatically adjusted to guarantee that all provided digits are taken into account, unless explicitely stated:
 - `gem('123456789012345678901234567890123456789012345678901234567891') - gem('123456789012345678901234567890123456789012345678901234567890')` produces the result

        1
    even though the difference is of the order of 1e-60. In this case, the actual precision of the numbers is higher than the current working precision.

 - `gem('1234567890',3) - gem('1234567891',3)` produces

        0
    because the three first digits of both numbers are equal.

Note that the accuracy of the result of an operation is not only determined by the number of decimals that the library takes into account to define numbers, but it may also depend on the [number of elementary manipulations](https://en.wikipedia.org/wiki/Numerical_error) that were used to produce the result.

The `precision` function can be used to determine how many digits are used to describe a particular number. For instance,
 - `precision(gem('3.2'))` is

        50
    while `precision(gem('e',5))` is

        5


### Display precision

When displaying a high precision number, only part of the number is usually printed out. This is the default behavior, because long string of digits can quickly become cumbersome Note however that it is not because additional digits are not printed that they are not kept in memory. The function `gemDisplayPrecision` can be used to adjust the number of displayed digits. Alternatively, the display precision can also be directly specified as a second argument of the `display` function:
 - `display(gem('pi'),50)` prints 50 digits:

        3.1415926535897932384626433832795028841971693993751

A *negative* display precision prints out all digits in memory:
 - `display(gem('pi'),-1)` prints the same 50 digits as above.




## Functions

A list of the functions that can be applied to gem and sgem is available [here](functions.md). In general, these functions take the same arguments as their matlab counterparts, and produce the analog result. For instance, the `max` function applied to complex numbers returns the complex number with largest magnitude, or the one with largest angle if both magnitudes are equal.

For some functions, all possible behaviors are not yet implemented. For instance, mpower currently only supports powers of +1 and -1. Anyone wishing to have such feature is invited to mention so and eventually consider contributing to this open source project by implementing the missing behaviors and submitting a push request (see also [this page](howToContribute.md) for more ways to get involved).

The GEM functions do differ from Matlab's implementation for some functions which don't preserve the sparsity when applied to sparse objects. This point is explained in the `gemSparseLikeMatlab.m` file. This behavior is easily circumvented by calling `gemSparseLikeMatlab(1)`.

