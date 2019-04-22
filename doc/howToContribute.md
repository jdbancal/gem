How to contribute
-----------------

Here is a detailed overview of the steps to be followed by anyone who would like to contribute to the GEM library.

Improvements to the code should be performed on a fork of the github repository (see [here](https://docs.gitlab.com/ce/workflow/forking_workflow.html) for information on git related collaborative programming). In case you are wondering what kind of features could be still worth implementing in the library, a few ideas are listed at the [bottom](#Desired Features) of this file.

Assuming you want to implement a particular function and already forked a version of the source code to work on, here is how you can proceed:
 - First of all, identify some general feature of the function you are considering:
    - How many parameters does the function depends on? Which of these parameters can be gem objects, which ones must be indices? If there are several input gem parameters, do they all need to be of the same size, or can the function mix matrices and scalars?
    - How many outputs does this function produce? Which ones are gem objects or indices?
    - Does the function preserve sparsity? (i.e. are there sparse inputs which can produce sparse outputs? This is not the case of the cosine functions, for instance.)
 - Identify an existing gem function which has identical or similar properties (e.g. `tan(x)` is similar to `sin(x)`, `plus(x,y)` is _not_ similar to `find(x)`), and read the code related to this function (both MATLAB/GNU Octave and c++ code).
 - Copy this function and rename it to your new function. Do this in the `@gem` folder, `src/gem_mex.cpp` file, `src/gem.hpp` file, `src/gem.cpp` file first.
 - Now you can modify these new pieces of code to do what you wish.
 - Compile your code and test it.
 - It should now be straightforward to implement your function for sparse matrices as well by copying and modifying the files/functions in the `@sgem` folder, `src/sgem_mex.cpp` file, `src/sgem.hpp` file, `src/sgem.cpp` file. If your function never preserves sparsity, no modification of the src folder is required at this stage, though: you only need to perform modification of the matlab code in the `@sgem` folder.
 - Compile and check your code.
 - Make sure that your code has a minimal amount of help information (including e.g. a usage example) and that it contains a few helpful comments that explain what is happening.
 - When all is fine, send a pull request on github to add your new feature to the library!


Design considerations
---------------------

- By default, truly sparse operations are implemented only if there is a chance that the result of the operation applied to a sparse matrix produces a sparse result. This diverges from matlab's default behavior. However, matlab's default behavior can be restored through the function 'gemSparseLikeMatlab'.

This means that `sin(x)` has a sparse implementation, but not `cos(x)` (`sin(0) = 0`, but `cos(0)` is not `0`). The implementation of the cosine function for sparse objects thus starts by a casting of the input into a full matrix, before calling the full implementation of the cosine function (see `@sgem/cos.m` for an example). Also the matrix inverse function `inv(X)` admits a sparse implementation, even though the inverse of most sparse matrices is not sparse. This is because there exist sparse matrices X whose inverse is also sparse (e.g. `X = eye(d)`).

- The library uses the data type *mpreal* provided by the mpfrc++ library. Therefore, it deals with complex numbers by itself, storing the real and imaginary components of complex objects in two different matrices. In particular, all interactions with the Eigen library involves purely real numbers.

Note that relying on std::complex is risky, because several algorithms assume that `std::complex` comes with double precision, hence leading to a loss of precision (c.f. comment from April 20, 2016 on http://www.holoborodko.com/pavel/mpfr/). Since the complex algebra is contained in the matrix algebra, complex operations can be implemented as matrix operations.

- The fact that the library interfaces with MATLAB/GNU Octave requires a dynamical allocation of objects. This is taken care of by a class adapted from [this c++ class interface](https://fr.mathworks.com/matlabcentral/fileexchange/38964-example-matlab-class-wrapper-for-a-c++-class). One consequence is that most methods in the c++ code have two variants :
  - `method` which performs the required computation and returns the result in a static variable (this variable is cleared once the library returns to matlab).
  - `method_new` which performs the same operation but returns the result in a dynamic variable (this variable remains in memory between to calls from matlab).

  The code of both functions must be identical, except for the variable definition.

- The MATLAB/GNU Octave is responsible of performing all parameters checks before calling the c++ library. So for instance it must check that the size of two matrices are compatible before asking the c++ library to multiply these matrices. No such checks should be expected to be performed on the c++ side. This is meant to lighten the c++ code, which is already complex enough as is.


Desired Features
----------------

Here is a partial list of features/functions that would be nice to add to the library.

 - improve the function `sprintf` to provide all digits
 - Implement the matrix `mpower` function for powers different from +/-1.
 - Implement the matrix exponential function `expm`
 - Parallelize the for loops appearing in simple functions such as `sin`.
 - `triu`, `tril`
 - Beta special functions: `beta`, `betainv`, `betainc`, `betaincinv`
 - For more ways to contribute, check if there is anything more [here](http://github.com/jdbancal/gem/issues)


