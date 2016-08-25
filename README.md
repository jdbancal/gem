Gmp Eigen Matrix Library
========================

The GEM library is a variable precision library for [matlab](http://www.mathworks.com/products/matlab/). It provides a simple way to perform matlab computations with a higher precision.

The library provides two data types:
 - 'gem' for high precision dense matrices
 - 'sgem' for high precision sparse matrices
and overloads numerous matlab functions.

The library is coded in C++. It relies on [GMP](https://gmplib.org/) for the high precision arithmetic (through [MPFR C++](http://www.holoborodko.com/pavel/mpfr/) and [MPFR](http://www.mpfr.org/)) and on [Eigen](http://eigen.tuxfamily.org/) for matrix manipulations.

At the moment, priority is given to functionality over performance, but the library still performs faster than matlab's builtin vpa type. For instance, multiplication of two 100x100 dense matrices with 100-digits precision is 10 times faster with gem objects compared to matlab 2016a's vpa type.

The GEM library is free and open source. Anyone can contribute on the [gitlab page](https://gitlab.com/jdbancal/gem). It is distributed under a ... license.


Usage example
-------------

