How to compile the GEM library
==============================

This document describes steps which lead to a successful compilation of the GEM library on various platforms.

Here we assume that the target environment is MATLAB, but the same instructions also apply to GNU Octave.


Steps to compile the GEM library on *Ubuntu* :
----------------------------------------------

The following steps work equally well for MATLAB and GNU Octave.

1. Download the latest library repository with the following command: `git clone https://www.github.com/jdbancal/gem`. This creates a folder called `gem`.
2. Inside this folder, download the submodules with the command `git submodule update --init`.
3. Install the gmp, mpfr and mpfrc++ libraries with the command
`sudo apt-get install libmpfrc++-dev libgmp-dev`
4. Start MATLAB or GNU Octave and go to the gem folder.
5. Type `make`. This will launch the compilation of the library. If everything goes fine, the program will conclude with 'Compilation successful'. Otherwise, a message should inform you about what is missing.
6. Add the gem folder to your MATLAB/GNU Octave path. You can now perform your favorite computation in high precision!


Steps to compile the GEM library on *Windows* (64 bits) :
---------------------------------------------------------

The following was only tested on MATLAB.

1. Download the latest library repository either with git (see the first Ubuntu instruction above) of from the [release page](https://github.com/jdbancal/gem/releases).
2. Install msys on your system from [http://mingw.org/wiki/msys](http://mingw.org/wiki/msys)
3. Launch the MinGW Installation Manager and install packages `mingw-developer-toolkit` and `msys-base` 
4. Download GMP from [https://gmplib.org/#DOWNLOAD](https://gmplib.org/#DOWNLOAD) and place it into gem's 'external' folder
5. Download MPFR from [http://www.mpfr.org/mpfr-current/#download](http://www.mpfr.org/mpfr-current/#download) and place it into gem's 'external' folder
6. Download MPFR-C++ from [http://www.holoborodko.com/pavel/mpfr/](http://www.holoborodko.com/pavel/mpfr/
) and place it into gem's 'external' folder
7. Create the folder 'external/staticLibraries'
8. From within msys (launch C:\MinGW\msys\1.0\msys.bat), go into the external/gmp folder and run the following commands:
    - ``./configure --disable-shared --enable-static CFLAGS=-fPIC --with-pic --prefix=`pwd`/../staticLibraries``
    - `make`
    - `make check`
    - `make install`
9. From within msys, go into the external/mpfr folder and run the following commands:
    - ``./configure --disable-shared --enable-static CFLAGS=-fPIC --with-pic --prefix=`pwd`/../staticLibraries --with-gmp=`pwd`/../staticLibraries``
    - `make`
    - `make check`
    - `make install`
10. Install the TDM64-GCC compiler from [http://tdm-gcc.tdragon.net/download](http://tdm-gcc.tdragon.net/download), making sure you include the component `Components/gcc/openmp` in the TDM-GCC Setup. Configure it for MATLAB by typing `setenv('MW_MINGW64_LOC','C:\TDM-GCC-64')` and `mex -setup cpp` in MATLAB (see [here](https://fr.mathworks.com/help/matlab/matlab_external/compiling-c-mex-files-with-mingw.html) for more details).
11. Download the latest Eigen source code on [eigen.tuxfamily.org](http://eigen.tuxfamily.org) and place it into gem's 'external' folder.
12. Download the latest version of Spectra on [http://yixuan.cos.name/spectra/download.html](http://yixuan.cos.name/spectra/download.html) and place it into gem's 'external' folder.
13. In MATLAB, go the the gem folder and type `make(1,0)` to compile the library.
14. Note: the compiled require the `libgomp` library, located in `C:\TDM-GCC-64\bin\libgomp_64-1.dll`. To make your compilation portable, copy this file into the `gem` folder.


Steps to compile the GEM library on *macOS* (tested on el capitain 10.11.6) :
-----------------------------------------------------------------------------

### Instructions for GNU Octave
These instructions were tested for a version of GNU Octave installed as an [App Bundle](https://octave-app.org/Download.html).

1. Download the latest library repository with the following command: `git clone https://www.github.com/jdbancal/gem`. This creates a folder called `gem`.
2. Inside this folder, download the submodules with the command `git submodule update --init`.
3. Download GMP from [https://gmplib.org/#DOWNLOAD](https://gmplib.org/#DOWNLOAD) and place it into gem's 'external' folder
4. Download MPFR from [http://www.mpfr.org/mpfr-current/#download](http://www.mpfr.org/mpfr-current/#download) and place it into gem's 'external' folder
5. Download MPFR-C++ from [http://www.holoborodko.com/pavel/mpfr/](http://www.holoborodko.com/pavel/mpfr/
) and place it into gem's 'external' folder
6. Start GNU Octave and go to the gem folder.
7. Type `make(1,0)`. This will launch the compilation of the library. If everything goes fine, the program will conclude with 'Compilation successful'. Otherwise, a message should inform you about what is missing.
8. Add the gem folder to your GNU Octave path. You can now perform your favorite computation in high precision!


### Instructions for MATLAB

Tentative roadmap:
1. Install Xcode
2. Install homebrew and most of the packages mentioned in https://www.topbug.net/blog/2013/04/14/install-and-use-gnu-command-line-tools-in-mac-os-x/
3. Install gcc throught homebrew with `brew install gcc --without-muiltilib` (or maybe without the option...)
4. Install maybe also `gmp` `mpfr` and `libmpc` with brew...
4. In /usr/local/bin, add reference to gcc and g++ with `ln -s gcc-7 gcc` and `ln -s g++-7 g++`
5. Launch MATLAB from the terminal, so that /usr/local/bin is the first folder looked for when running `unix(‘which g++’)` in MATLAB
6. Type `make` from  within the gem folder and follow the instructions
7. Add the gem folder to your MATLAB path.

These additional links may help
 - Instructions to be able to run `mex -setup` : http://www.mathworks.com/matlabcentral/answers/246507-why-can-t-mex-find-a-supported-compiler-in-matlab-r2015b-after-i-upgraded-to-xcode-7-0#comment_392485

