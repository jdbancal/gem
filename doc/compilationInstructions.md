How to compile the GEM library
==============================

This document describes steps which lead to a successful compilation of the GEM library on various platforms.


Steps to compile the GEM library on *ubuntu* :
----------------------------------------------

1. Download the [latest library release](https://github.com/jdbancal/gem/releases) in the folder of your choice
2. Download the latest Eigen source code on [eigen.tuxfamily.org](http://eigen.tuxfamily.org) and place it into gem's src folder.
3. Download the latest version of Spectra on [http://yixuan.cos.name/spectra/download.html](http://yixuan.cos.name/spectra/download.html) and place it into gem's src folder.
4. Install the gmp and mpfrc++ libraries with the command
`sudo apt-get install libmpfrc++-dev libgmp-dev`
5. Start matlab and go to the gem folder.
6. Type `make`. This will launch the compilation of the library. If everything goes fine, the program will conclude with 'Compilation successful'. Otherwise, a message should inform you about what is missing.
7. Add the gem folder to your matlab path. You can now perform your favorite computation in high precision !-)


Steps to compile the GEM library on *windows* (64 bits) :
---------------------------------------------------------

1. Install msys on your system from [http://mingw.org/wiki/msys](http://mingw.org/wiki/msys)
2. Launch the MinGW Installation Manager and install packages `mingw-developer-toolkit` and `msys-base` 
3. Download GMP from [https://gmplib.org/#DOWNLOAD](https://gmplib.org/#DOWNLOAD) and place it into gem's src folder
4. Download MPFR from [http://www.mpfr.org/mpfr-current/#download](http://www.mpfr.org/mpfr-current/#download) and place it into gem's src folder
5. Download MPFR-C++ from [http://www.holoborodko.com/pavel/mpfr/](http://www.holoborodko.com/pavel/mpfr/
) and place it into gem's src folder
6. Create the folder src/staticLibraries
7. From within msys (launch C:\MinGW\msys\1.0\msys.bat), go into the src/gmp folder and run the following commands:
    - `./configure --disable-shared --enable-static CFLAGS=-fPIC --prefix=`pwd`/../staticLibraries`
    - `make`
    - `make check`
    - `make install`
8. From within msys, go into the src/mpfr folder and run the following commands:
    - `./configure --disable-shared --enable-static CFLAGS=-fPIC --prefix=`pwd`/../staticLibraries --with-gmp=`pwd`/../staticLibraries`
    - `make`
    - `make check`
    - `make install`
9. Install the TDM64-GCC compiler from [http://tdm-gcc.tdragon.net/download](http://tdm-gcc.tdragon.net/download), making sure you include the component `Components/gcc/openmp` in the TDM-GCC Setup. Configure it for matlab by typing `setenv('MW_MINGW64_LOC','C:\TDM-GCC-64')` and `mex -setup cpp` in matlab (see (here)[https://fr.mathworks.com/help/matlab/matlab_external/compiling-c-mex-files-with-mingw.html] for more details).
10. In the file `make.m`, set the variabe `useSharedGmpAndMpfr` to 0 (it takes value 1 by default)
11. Follow the instructions above for the installation on ubuntu (except for point number 4)
12. Note: the compiled require the `libcomp` library, located in `C:\TDM-GCC-64\bin\libcomp-1.dll`.by default: copy this file into the `gem` folder


Steps to compile the GEM library on *MacOs* (tested on el capitain 10.11.6) :
-----------------------------------------------------------------------------

Tentative roadmap:
1. Install Xcode
2. Install homebrew and most of the packages mentioned in https://www.topbug.net/blog/2013/04/14/install-and-use-gnu-command-line-tools-in-mac-os-x/
3. Install gcc throught homebrew with `brew install gcc --without-muiltilib` (or maybe without the option...)
4. Install maybe also `gmp` `mpfr` and `libmpc` with brew...
4. In /usr/local/bin, add reference to gcc and g++ with `ln -s gcc-7 gcc` and `ln -s g++-7 g++`
5. Launch matlab from the terminal, so that /usr/local/bin is the first folder looked for when running `unix(‘which g++’)` in matlab
6. Type `make` from  within the gem folder and follow the instructions
7. Add the gem folder to your matlab path.

These additional links may help
 - Instructions to be able to run `mex -setup` : http://www.mathworks.com/matlabcentral/answers/246507-why-can-t-mex-find-a-supported-compiler-in-matlab-r2015b-after-i-upgraded-to-xcode-7-0#comment_392485

