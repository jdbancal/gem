Gmp Eigen Matrix Library
========================

The [GEM library](https://github.com/jdbancal/gem/releases) is a variable precision library for [MATLAB](http://www.mathworks.com/products/matlab/) and [GNU Octave](https://www.gnu.org/software/octave/). It provides an open source solution for basic high precision computations in standard numerical computing environments.

The library implements two data types:
 - **gem** for high precision dense matrices
 - **sgem** for high precision sparse matrices

and overloads [a number of MATLAB/GNU Octave functions](doc/functions.md). The full library can be downloaded [here](https://github.com/jdbancal/gem/releases).

The GEM library is coded in C++ and MATLAB/GNU Octave. It currently relies on [GMP](https://gmplib.org/) for high precision arithmetic (through [MPFR](http://www.mpfr.org/) and [MPFR C++](http://www.holoborodko.com/pavel/mpfr/)), and on [Eigen](http://eigen.tuxfamily.org/) and [Spectra](http://yixuan.cos.name/spectra/) for matrix manipulations.

At the moment, priority is given to functionality over performance: the code comes with no garantee of optimality. Nevertheless, appreciable improvement is already available compared to MATLAB's builtin types performance. With 100 digits of precision, for instance, a 100x100 matrix is transfered from _double to high precision_ 10x faster, from _high precision to double_ format 250x faster, and its _column-wise minimum_ is computed 25x faster. _Multiplication_ of two 100x100 dense matrices with 100-digits precision is 10 times faster with gem objects as compared to MATLAB 2016a's vpa type. For a matrix of size 1000x1000 these ratios become respectively 14x, 1500~20000x, 500x, 10x (the GEM library's multithreading capabilities were disactivated when performing these comparisons).

The GEM library also implements sparse matrix multiplication, which is not available with MATLAB's vpa type.


Usage examples
--------------
Here are a few simple examples of GEM library usage.

 - `toStrings(gem('pi',100))` prints the 100 first decimals of pi 
 - `gem(2)`, `gem(1.23)`, `gem([1 2 3])` create 50-digits precision representations of the numbers 2 and 1.23, and of the vector [1 2 3]. When translating a number from double form, exactly 15 digits are taken into account.
 - `gem('1.23456789123456789+2i')` creates a 50-digits representation of the number provided in text form (all digits within the working precision are taken into account for numbers specified in string format)
 - `gemWorkingPrecision(100)` changes the working precision to 100 digits

Once a high precision matrix is created, it can be manipulated by calling usual MATLAB/GNU Octave functions.

 - `eig(gemRand(100,100))` : computes the eigenvalues of a random 100x100 matrix
 - `eigs(gemRand(100,100),1)` : computes the largest eigenvalue of a random 100x100 matrix
 - `sum(gem([1:100000]).^8)-5e39-2e24` gives 111111111177777777771111111111333333333330000 (all digits can be displayed with the help of the function `gemDisplayPrecision`)
 - `notAnInteger = exp(sqrt(gem(163))*gem('pi')); display(notAnInteger, -1)` gives 262537412640768743.9999999999992500725971981856889 (the number of digits displayed can also be specified on a case by case fashion as a parameter to the display fuction; a precision of -1 displays all available digits)
 - `sgem(speye(3))` creates a high precision sparse representation of the 3x3 identity matrix
 - `a=1./gem([1:7]); save('filename','a'); load('filename');` saves and loads a gem object
 - Watch how how 50 digits of the golden ratio are being calculated by the continued fraction formula with `x= gem(1); for i=1:118 x=1/(1+x); disp(num2str(1+x,50)); end; disp(num2str((1+sqrt(gem(5)))/2, 50))`

Check out [getting started with the GEM library](doc/gettingStarted.md) for more usage examples.


Installation
------------

The library comes pre-compiled for linux, macos and windows (64bits). It is therefore straightforward to use : after [downloading](https://github.com/jdbancal/gem/releases) the latest release, just add the gem folder into MATLAB/GNU Octave's path (this can be done by running the command `path(path,'/path_to_the_gem_folder')`), and it is ready to use.

If you experience trouble with the provided binaries, because for instance you are using a different platform (32 bits?), a significantly older version of the operating system than the one on which the provided binaries were compiled, a significantly older version of MATLAB or GNU Octave, you may use the provided script to compile it. For more details on this, please refer to the [compilation instructions](doc/compilationInstructions.md).


License
-------

The GEM library is free and open source. It can be used for both open-source and proprietary application. Therefore, it is also free for academic use. Anyone can [contribute](doc/howToContribute.md) on the [github page](https://github.com/jdbancal/gem). The source code is distributed under a MPL2 license. See [COPYING.md](COPYING.md) for more details.

