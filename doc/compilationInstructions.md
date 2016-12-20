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

1. Install msys on your system
2. Install the TDM-GCC compiler from [http://tdm-gcc.tdragon.net/download](http://tdm-gcc.tdragon.net/download) and configure it for matlab (see (here)[https://fr.mathworks.com/help/matlab/matlab_external/compiling-c-mex-files-with-mingw.html] for more details).
3. Download GMP from [https://gmplib.org/#DOWNLOAD](https://gmplib.org/#DOWNLOAD) and place it into gem's src folder
4. Download MPFR from [http://www.mpfr.org/mpfr-current/#download](http://www.mpfr.org/mpfr-current/#download) and place it into gem's src folder
5. Download MPFR-C++ from [http://www.holoborodko.com/pavel/mpfr/
](http://www.holoborodko.com/pavel/mpfr/
) and place it into gem's src folder
6. Create the folder src/staticLibraries
7. From within msys, go into the src/gmp folder and run the following commands:
    - `./configure --disable-shared --enable-static CFLAGS=-fPIC --prefix=`pwd`/../staticLibraries`
    - `make`
    - `make check`
    - `make install`
8. From within msys, go into the src/mpfr folder and run the following commands:
    - `./configure --disable-shared --enable-static CFLAGS=-fPIC --prefix=`pwd`/../staticLibraries --with-gmp=`pwd`/../staticLibraries`
    - `make`
    - `make check`
    - `make install`
9. In the file `make.m`, set the variabe `useSharedGmpAndMpfr` to 0 (it takes value 1 by default)
10. Follow the instructions above for the installation on ubuntu (except for point number 4)


Steps to compile the GEM library on *MacOs* and other platform :
----------------------------------------------------------------

Compilation on *other platforms* : unfortunately, I don't have other operating systems on which to run tests. Since the compilation commands are relatively simple (not much configuration involved), it should be rather straightforward to compile the GEM library with something similar to the above commands. The best way to start is to make sure that a running version of [gcc](https://gcc.gnu.org/) is installed on your machine. Once this is done, you can try running the compilation commands one by one, and see what happens. If you succeed, tell us how you did it!


