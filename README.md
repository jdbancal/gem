Gmp Eigen Matrix Library
========================

The GEM library is an open source variable precision library for [matlab](http://www.mathworks.com/products/matlab/). It provides a simple way to perform matlab computations with a high precision.

500-fold overheads when performing simple comparison operations (c.f. below) and the absence of support for high-precision sparse matrices are hardly admissible constraints for building innovative high precision algorithms, just like proprietary code. The GEM library is an attempt to address these problems by providing high precision capabilities for basic matlab functions in an open source way.

The library provides two data types:
 - 'gem' for high precision dense matrices
 - 'sgem' for high precision sparse matrices
and overloads numerous matlab functions.

The library is coded in C++. It relies on [GMP](https://gmplib.org/) for the high precision arithmetic (through [MPFR C++](http://www.holoborodko.com/pavel/mpfr/) and [MPFR](http://www.mpfr.org/)) and on [Eigen](http://eigen.tuxfamily.org/) for matrix manipulations.

At the moment, priority is given to functionality over performance. Nevertheless, appreciable perforance improvement is already availble compared to matlab's builtin vpa type. With 100 digits of precision, for instance, a 100x100 matrix is transfered from double to high precision 10x faster, from high precision to double format 250x faster, and its column-wise minimum is computed 25x faster. Multiplication of two 100x100 dense matrices with 100-digits precision is 10 times faster with gem objects compared to matlab 2016a's vpa type. For a matrix of size 1000x100 these ratios become respectively 14x, 1500~20000x, 500x, 10x (with multithreading desactivated).

The GEM library is free and open source. It is therefore also free for academic use. Anyone can contribute on the [gitlab page](https://gitlab.com/jdbancal/gem) (see below for more details). It is distributed under a ... license.


Usage examples
--------------

 - gem(2), gem(1.23) create 50-digits precision representations of the numbers 2 and 1.23. When translating a number from double form, exactly 15 digits are taken into account.
 - gem('1.23456789123456789+2i') creates a 50-digits representation of the number provided in text form (all digits within the working precision are taken into account
 - gemWorkingPrecision(100) updates the working precision to 100 digits
 - eig(gemRand(100,100)) : computes the eigenvalues of a random 100x100 matrix
 - notAnInteger = exp(sqrt(gem(163))*gem('pi')); display(notAnInteger, -1) gives 262537412640768743.9999999999992500725971981856889 (a precision of -1 displays all available digits)
 - sgem(eye(3)) creates a high precision sparse representation of the 3x3 identity matrix
 - a=1./gem([1:7]); save('filename','a'); load('filename'); saves and loads a gem object


Installation
------------

This library comes pre-compiled for ubuntu 64bits. It is therefore straightforward to use : just add the gem folder into matlab's path, and it is ready to use. This can be done by running the command '''path(path,'/path_to_the_gem_folder/gem').

If you are using a different platform (32 bits, mac os or windows), or in case you use an older version of linux/matlab than the one on which this was compiled, you need to compile it. Below are instructions


Installation (compilation)
--------------------------

Here are the instructions for compiling GEM on *ubuntu* (please update this file on [gitlab](https://gitlab.com/jdbancal/gem) with instructions if you find out how to install it on other platforms, unfortunately I don't have any of those to play with).

1. Check out this repository in the folder of your choice
2. Download the latest version of Eigen on [eigen.tuxfamily.org](eigen.tuxfamily.org) and place it into the src folder.
3. Install the gmp and mpfrc++ libraries with the command
'''sudo apt-get install libmpfrc++-dev libgmp-dev
4. Start matlab and go to the gem folder.
5. Type 'make'. This will launch the compilation of the library. If everything goes fine, the program will return 'Compilation successful'.
6. Add the gem folder to your matlab path. You can now safely play with high precision objects :-)



How to contribute
-----------------

Here is a detailed overview of the steps to follow if you want to add one function to the GEM library:

 - First of all, find out:
    - How many parameters does the function depends on? Which of these parameters can be gem objects, which ones must be indices? If there are several input gem parameters, do they all need to be of the same size, or can the function mix matrices and scalars?
    - How many outputs does this function produce? Which ones are gem objects or indices?
    - Does the function preserve sparsity? (i.e. are there sparse inputs which can produce sparse outputs? This is not the case of the cosine functions, for instance.)
 - Identify an existing gem function which has identical or similar properties (e.g. tan(x) is similar to sin(x), plus(x,y) is not similar to find(x))
 - Copy this function and rename it to you new function. Do this in the @gem folder, src/gem_mex.cpp file, src/gem.hpp file, src/gem.cpp file first.
 - Now you can modify these new pieces of code to do what you wish.
 - Compile your code and test it.
 - Now you can implement you function for sparse matrices as well by copying and modifying the files/functions in the @sgem folder, src/sgem_mex.cpp file, src/sgem.hpp file, src/sgem.cpp file. If your function never preserves sparsity, you only need to perform modification at the matlab code in the @sgem folder (i.e. no modification in the src folder at this stage).
 - Compile and check your code.
 - Make sure that your code has a minimal amount of help information and that it contains a few helpful comments that explain what is happening.
 - When all is fine, send a pull request on gitlab to add your new feature to the library!



Included functions
------------------


| function  |   |   |   |   |
|---|---|---|---|---|
| 1  |   |   |   |   |
| 2  |   |   |   |   |
| 3  |  rwe |   |   |   |



Design considerations
---------------------

- The library uses the type 'mpreal' provided by the mpfrc++ library. Therefore, it deal with complex numbers by itself. In particular, all interactions with the Eigen library involves purely real numbers.

Note that relying on std::complex has been shown to lead to problems, because several algorithms assume that 'std::complex' comes with double precision, hence leading to a loss of precision (c.f. comment from April 20, 2016 on http://www.holoborodko.com/pavel/mpfr/).

- Truly sparse operations are implemented only if there is a chance that the result of the operation applied to a sparse matrix produces a sparse result. This diverges from matlab's default behavior. However, matlab's default behavior can be restored through the function 'gemSparseLikeMatlab'.

This means that sin(x) has a sparse implementation, but not cos(x) (sin(0) = 0, but cos(0) is not 0). Also the matrix inverse function inv(X) admits a sparse implementation, even though the inverse of most sparse matrices is not sparse. This is because there exist sparse matrices X whose inverse is also sparse (e.g. X = eye).




