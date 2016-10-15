Getting Started with the GEM Library
====================================

Here is a short introduction to the usage of GEM library in matlab. This introduction assumes that the gem folder is in matlab's path and that it is working (this can be checked by typing e.g. `gem('pi')`; if this command yields an error, check that the gem library is compiled for your system and that its folder is in matlab's path).

### Precision

The GEM library works with numbers that can have a large precision. This precision is defined by the number of digits that are used to describe numbers. By default, the library takes 50 digits into account. This means that the library can distinguish between numbers that differ at the 50th decimal place:

`gem('12345678901234567890123456789012345678901234567891') - gem('12345678901234567890123456789012345678901234567890')` produces the result
> 1

The precision at which the GEM library works can be adjusted throught the function `gemWorkingPrecision`.

Note that the actual precision of the result of an operation is not only determined by the number of decimals that the library takes into account to define numbers, but it may also depend on the [number of elementary manipulations](https://en.wikipedia.org/wiki/Numerical_error) that were used to produce the result.

If you enter only one of the two gem numbers above in the matlab prompt, you will see that only part of the number is printed out. This is the default behavior, because long string of digits can quickly become cumbersome (it is not because additional digits are not printed that they are not in memory). The function `gemDisplayPrecision` can be used to adjust the number of digits printed. Alternatively, the display precision can also be directly specified as a second argument of the `display` function:

`display(gem('pi'),50)` prints 50 digits
> 3.1415926535897932384626433832795028841971693993751

A /negative/ display precision prints out all digits in memory:

    `display(gem('pi'),-1)` prints the same 50 digits as above.




### ...


 - `gem(2)`, `gem(1.23)` create 50-digits precision representations of the numbers 2 and 1.23. When translating a number from double form, exactly 15 digits are taken into account.
 - `gem('1.23456789123456789+2i')` creates a 50-digits representation of the number provided in text form (all digits within the working precision are taken into account
 - `gemWorkingPrecision(100)` updates the working precision to 100 digits
 - `eig(gemRand(100,100))` : computes the eigenvalues of a random 100x100 matrix
 - `sum(gem([1:100000]).^8)-5e39` gives 111111111177777777773111111111333333333330000
 - `notAnInteger = exp(sqrt(gem(163))*gem('pi')); display(notAnInteger, -1)` gives 262537412640768743.9999999999992500725971981856889 (a precision of -1 displays all available digits)
 - `sgem(eye(3))` creates a high precision sparse representation of the 3x3 identity matrix
 - `a=1./gem([1:7]); save('filename','a'); load('filename');` saves and loads a gem object


