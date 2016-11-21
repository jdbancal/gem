function make(parallelization)
% make(parallelization)
%
% Run this file from matlab to compile the C++ part of the GEM library. Use
% parallelization = 0 to disable openmp, parallelization = 1 otherwise.
%
% Note : the code relies on the eigen, spectra, gmp and mpfrc++ libraries.
% - The Eigen library (version 3.2 (!)) should be dowloaded from http://eigen.tuxfamily.org 
%   and placed in the gem/src folder.
% - The Spectra library should be downloaded from 
%   http://yixuan.cos.name/spectra/index.html and places in the gem/src 
%   folder.
% On ubuntu, the two remaining libraries can be installed by running the 
% following command in the terminal :
%   sudo apt-get install libgmp-dev libmpfr-dev libmpfrc++-dev
%
% It is also possible to perform a compilation with gmp and mpfr source 
% codes instead of relying on precompiled gmp and mpfr libraries. In this
% case proceed as follow :
% - Download the gmp source code should from https://gmplib.org/#DOWNLOAD
%   and unpack it into the 'src' folder (this creates the folder src/gmp*)
% - Download the mpfr source code from http://www.mpfr.org/mpfr-current/#download
%   and unpack it into the 'src' folder (this creates the folder src/mpfr*)
% - Download the mpfrc++ source from http://www.holoborodko.com/pavel/mpfr/
%   and place it into the 'src' folder (this creates the folder src/mpfrc++*)
% - Set the value of 'useSharedGmpAndMpfr' to 0 in this file
% - Run this file from matlab


%% Settings
if nargin < 1
    parallelization = 1;
end

% Set the following variable to 0 in order to compile a library which does
% not depend on shared gmp and mpfr libraries:
useSharedGmpAndMpfr = 1;


%% We perform some checks
% This file needs to be run from within its folder
tmp = mfilename('fullpath');
tmp = tmp(1:find(tmp=='/',1,'last')-1);
if ~isequal(tmp, pwd)
    error('Unable to compile the gem library from a different folder.');
end

% We extract the matlab version
v = version;
vDots = find(v=='.');
v = str2num(v(1:vDots(2)-1));


% We check that the eigen and spectra source code were downloaded
list = dir('src');
eigenFound = 0;
spectraFound = 0;
for i=1:length(list)
    eigenFound = eigenFound + isequal(list(i).name,'eigen');
    spectraFound = spectraFound + isequal(list(i).name,'spectra');
end
if eigenFound == 0
    error('The eigen folder was not found in the gem/src folder. Please download the eigen library from http://eigen.tuxfamily.org');
end
if spectraFound == 0
    error('The spectra folder was not found in the gem/src folder. Please download it from http://yixuan.cos.name/spectra/index.html');
end


% We set the flags needed for parallelization (if requested)
if parallelization == 1
    parallelizationFlag = '-pthread -fopenmp';
else
    parallelizationFlag = '';
end


% additional flags
optimizationFlag = '-O3';


%% Now we proceed to compilation
if useSharedGmpAndMpfr == 1
    flags(1) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Isrc/spectra/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG src/utils.cpp -o src/utils.o']);
    flags(2) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Isrc/spectra/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG src/gem.cpp -o src/gem.o']);
    flags(3) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Isrc/spectra/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG src/gem_mex.cpp -o src/gem_mex.o']);
    flags(4) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Isrc/spectra/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG src/sgem.cpp -o src/sgem.o']);
    flags(5) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Isrc/spectra/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG src/sgem_mex.cpp -o src/sgem_mex.o']);
    flags(6) = unix(['g++ -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/glnxa64/mexFunction.map" src/gem_mex.o src/gem.o src/sgem.o src/utils.o   -lmpfr  -lgmp   -Wl,-rpath-link,', matlabroot, '/bin/glnxa64 -L"', matlabroot, '/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -o gem_mex.', mexext]);
    flags(7) = unix(['g++ -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/glnxa64/mexFunction.map" src/sgem_mex.o src/sgem.o src/gem.o src/utils.o   -lmpfr  -lgmp   -Wl,-rpath-link,', matlabroot, '/bin/glnxa64 -L"', matlabroot, '/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -o sgem_mex.', mexext]);
    if sum(flags) == 0
        display('Compilation successful.');
    else
        error('Compilation error.');
    end
else
    % Here we compile a version of the library which does not rely on gmp 
    % and mpfr shared libraries. It requires additional source code :
    % - download the gmp source code from https://gmplib.org/#DOWNLOAD
    %   and unpack it into the 'src' folder (this creates the folder src/gmp*)
    % - download the mpfr source code from http://www.mpfr.org/mpfr-current/#download
    %   and unpack it into the 'src' folder (this creates the folder src/mpfr*)
    % - download the mpfrc++ source code from http://www.mpfr.org/mpfr-current/#download
    %   and unpack it into the 'src' folder (this creates the folder src/mpfrc++*)
    % - Download the spectra source from http://yixuan.cos.name/spectra/index.html
    %   and place it into the 'src' folder (this creates the folder src/spectra)

    % We check that the gmp, mpfr and mpfrc++ sources are available
    gmpFolder = '';
    mpfrFolder = '';
    mpfrcppFolder = '';
    spectraFolder = '';
    for i=1:length(list)
        if isequal(lower(list(i).name(1:min(end,3))),'gmp')
            gmpFolder = list(i).name;
        end
        if isequal(lower(list(i).name(1:min(end,4))),'mpfr') && ~isequal(lower(list(i).name(1:min(end,7))),'mpfrc++')
            mpfrFolder = list(i).name;
        end
        if isequal(lower(list(i).name(1:min(end,7))),'mpfrc++')
            mpfrcppFolder = list(i).name;
        end
        if isequal(lower(list(i).name(1:min(end,7))),'spectra')
            spectraFolder = list(i).name;
        end
    end
    if isempty(gmpFolder)
        error('The gmp source code was not found. Please download it from https://gmplib.org and place it into the gem/src folder');
    end
    if isempty(mpfrFolder)
        error('The mpfr source code was not found. Please download it from http://www.mpfr.org and place it into the gem/src folder');
    end
    if isempty(mpfrcppFolder) || ~isequal(exist(['src/', mpfrcppFolder, '/mpreal.h']), 2)
        error('The mpfrc++ source code was not found. Please download it from http://www.holoborodko.com/pavel/mpfr and place it into the gem/src folder');
    end
    if isempty(spectraFolder)
        error('The spectra source code was not found. Please download it from http://yixuan.cos.name/spectra/index.html and place it into the gem/src folder');
    end
    
    % We compile the libraries if needed
    cd src
    if exist('staticLibraries') ~= 7
        unix('mkdir staticLibraries');
    end
    if exist('staticLibraries/lib/libgmp.a') ~= 2
        cd(gmpFolder);
        unix('./configure --disable-shared --enable-static CFLAGS=-fPIC --prefix=`pwd`/../staticLibraries');
        unix('make');
        unix('make check');
        unix('make install');
        cd ..
    end
    if exist('staticLibraries/lib/libmpfr.a') ~= 2
        cd(mpfrFolder);
        unix('./configure --disable-shared --enable-static CFLAGS=-fPIC --prefix=`pwd`/../staticLibraries --with-gmp=`pwd`/../staticLibraries');
        unix('make');
        unix('make check');
        unix('make install');
        cd ..
    end

    % Now gmp and mpfr should be ready, we can compile gem
    flags(1) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -I', mpfrcppFolder, ' -IstaticLibraries/include -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG utils.cpp -o utils.o']);
    flags(2) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -I', mpfrcppFolder, ' -IstaticLibraries/include -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG gem.cpp -o gem.o']);
    flags(3) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -I', mpfrcppFolder, ' -IstaticLibraries/include -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG gem_mex.cpp -o gem_mex.o']);
    flags(4) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -I', mpfrcppFolder, ' -IstaticLibraries/include -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG sgem.cpp -o sgem.o']);
    flags(5) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -I', mpfrcppFolder, ' -IstaticLibraries/include -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG sgem_mex.cpp -o sgem_mex.o']);
    flags(6) = unix(['g++ -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/glnxa64/mexFunction.map" gem_mex.o gem.o sgem.o utils.o  staticLibraries/lib/libmpfr.a staticLibraries/lib/libgmp.a  -Wl,-rpath-link,', matlabroot, '/bin/glnxa64 -L"', matlabroot, '/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -o ../gem_mex.', mexext]);
    flags(7) = unix(['g++ -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/glnxa64/mexFunction.map" sgem_mex.o sgem.o gem.o utils.o  staticLibraries/lib/libmpfr.a staticLibraries/lib/libgmp.a  -Wl,-rpath-link,', matlabroot, '/bin/glnxa64 -L"', matlabroot, '/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -o ../sgem_mex.', mexext]);
    cd ..
    if sum(flags) == 0
        display('Compilation successful.');
    else
        error('Compilation error.');
    end
end

return;

% The following line launches the compilation with default options (and no
% parallelization) :
%mex -largeArrayDims -Isrc/eigen -Isrc/eigen/unsupported -Isrc/spectra/include -I/usr/include -lmpfr -lgmp src/gem_mex.cpp src/gem.cpp src/utils.cpp

