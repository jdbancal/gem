Gmp Eigen Matrix Library
========================

The GEM library is an open source variable precision library for [matlab](http://www.mathworks.com/products/matlab/). It provides a simple way to perform matlab computations with a higher precision.

The library provides two data types:
 - 'gem' for high precision dense matrices
 - 'sgem' for high precision sparse matrices
and overloads numerous matlab functions.

The library is coded in C++. It relies on [GMP](https://gmplib.org/) for the high precision arithmetic (through [MPFR C++](http://www.holoborodko.com/pavel/mpfr/) and [MPFR](http://www.mpfr.org/)) and on [Eigen](http://eigen.tuxfamily.org/) for matrix manipulations.

At the moment, priority is given to functionality over performance, but the library still performs faster than matlab's builtin vpa type. For instance, multiplication of two 100x100 dense matrices with 100-digits precision is 10 times faster with gem objects compared to matlab 2016a's vpa type.

The GEM library is free and open source. Anyone can contribute on the [gitlab page](https://gitlab.com/jdbancal/gem). It is distributed under a ... license.


Usage example
-------------



Installation (compilation)
--------------------------

Here are the instructions for installing gem on *ubuntu* (please update this file on [gitlab](https://gitlab.com/jdbancal/gem) with instructions if you find out how to install it on other platforms, unfortunately I don't have any of those to play with).

1. Check out this repository in the folder of your choice
2. Download the latest version of Eigen on [eigen.tuxfamily.org](eigen.tuxfamily.org) and place it into the src folder.
3. Install the gmp and mpfrc++ libraries with the command
'''sudo apt-get install libmpfrc++-dev libgmp-dev
4. Start matlab and go to the gem folder.
5. Type 'make'. This will launch the compilation of the library. If everything goes fine, the program will return 'Compilation successful'.
6. Add the gem folder to your matlab path. You can now safely play with high precision objects :-)



