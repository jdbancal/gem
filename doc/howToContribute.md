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


