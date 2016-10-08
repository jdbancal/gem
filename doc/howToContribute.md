How to contribute
-----------------

Here is a detailed overview of the steps to follow if you want to add one function to the GEM library. In case you are done and wonder what else could be worth implementing, you can have a look at the list at the [bottom](#Desired Features) of this file.

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


Design considerations
---------------------

- Truly sparse operations are implemented only if there is a chance that the result of the operation applied to a sparse matrix produces a sparse result. This diverges from matlab's default behavior. However, matlab's default behavior can be restored through the function 'gemSparseLikeMatlab'.

This means that sin(x) has a sparse implementation, but not cos(x) (sin(0) = 0, but cos(0) is not 0). Also the matrix inverse function inv(X) admits a sparse implementation, even though the inverse of most sparse matrices is not sparse. This is because there exist sparse matrices X whose inverse is also sparse (e.g. X = eye).

- The library uses the type *mpreal* provided by the mpfrc++ library. Therefore, it deal with complex numbers by itself. In particular, all interactions with the Eigen library involves purely real numbers.

Note that relying on std::complex has been shown to lead to problems, because several algorithms assume that 'std::complex' comes with double precision, hence leading to a loss of precision (c.f. comment from April 20, 2016 on http://www.holoborodko.com/pavel/mpfr/).


Desired Features
----------------

Here is a list of some features/functions that would be nice to add to the library. (Yes, this section is not called _todo_, because no contributor is forced to do anything, the whole project relies on free contributions.)

 - Some sorting functions such as `sort`, `sortrow` and `unique`
 - Singular value decomposition `svd`, this would allow computing the 2-norm of a matrix as well as `cond`
 - Add eigenvalue decomposition for sparse matrices. This is most likely going to happen throught the function `eigs` (use [this library](https://github.com/yixuan/arpack-eigen)?)
 - Add linear system solvers (`\` operator) -- something like this is already done in the function `inv`.
 - Implement the matrix `mpow` function for powers different from +/-1.
 - Add a function that checks whether a matrix is `symmetric` or `hermitian`. Then, allow functions such as `eig`, `eigs`, `svd`, `inv` and `\` to adjust their algorithm choice accordingly.
 - Parallelize the for loops appearing in simple functions such as `sin`.


