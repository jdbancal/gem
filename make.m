function make(parallelization, useSharedGmpAndMpfr)
% make([parallelization], [useSharedGmpAndMpfr])
%
% Run this file from matlab to compile the C++ part of the GEM library. Use
% parallelization = 0 to disable openmp, parallelization = 1 otherwise. Use
% useSharedGmpAndMpfr = 1 if you want to use the GMP and MPFR libraries
% installed on your system, use useSharedGmpAndMpfr = 0 if you want to use
% the manually downloaded versions.
%
% Note : the code relies on the eigen, spectra, gmp, mpfr and mpfrc++ libraries.
% - The Eigen and Spectra libraries are included as git submodules. They
%   should be downloaded with the command
%     git submodule update --init
% On ubuntu, the three remaining libraries can be installed by running the 
% following command in the terminal :
%     sudo apt-get install libgmp-dev libmpfr-dev libmpfrc++-dev
%
% It is also possible to perform a compilation with gmp and mpfr source 
% codes instead of relying on precompiled gmp and mpfr libraries (this
% should be the way to go for other systems than ubuntu). In this case, 
% proceed as follow :
% - Download the gmp source code should from https://gmplib.org/#DOWNLOAD
%   and unpack it into the 'external' folder (this creates the folder external/gmp*)
% - Download the mpfr source code from http://www.mpfr.org/mpfr-current/#download
%   and unpack it into the 'external' folder (this creates the folder external/mpfr*)
% - Download the mpfrc++ source from http://www.holoborodko.com/pavel/mpfr/
%   and place it into the 'external' folder (this creates the folder external/mpfrc++*)
% - Call this script from matlab with the option 'useSharedGmpAndMpfr' set
%   to 0


%% Settings
if nargin < 1
    parallelization = 1;
end

% Set the following variable to 0 in order to compile a library which does
% not depend on shared gmp and mpfr libraries:
if nargin < 2
    useSharedGmpAndMpfr = 1;
end

warning('You are trying to compile the GEM library, but binaries may be available for your platform on https://www.github.com/jdbancal/gem/releases');


%% We perform some checks
% This file needs to be run from within its folder
tmp = mfilename('fullpath');
if ispc
    tmp = tmp(1:find(tmp=='\',1,'last')-1);
else
    tmp = tmp(1:find(tmp=='/',1,'last')-1);
end
if ~isequal(tmp, pwd)
    error('Unable to compile the gem library from a different folder.');
end

% We extract the matlab version
v = version;
vDots = find(v=='.');
v = str2num(v(1:vDots(2)-1));


% We check that the eigen and spectra source code were downloaded
list = dir('external');
eigenFolder = '';
spectraFolder = '';
for i=1:length(list)
    if isequal(lower(list(i).name(1:min(end,5))),'eigen')
        eigenFolder = list(i).name;
    end
    if isequal(lower(list(i).name(1:min(end,7))),'spectra')
        spectraFolder = list(i).name;
    end
end
if isempty(eigenFolder)
    error('The eigen folder was not found in the gem/external folder. Please run the command `git submodule update --init`.');
end
if isempty(spectraFolder)
    error('The spectra source code was not found in the gem/external folder. Please run the command `git submodule update --init`');
end
eigenFolder = ['external/', eigenFolder];
spectraFolder = ['external/', spectraFolder];


% We set the flags needed for parallelization (if requested and if possible)
if (parallelization == 1) && ~ismac
    parallelizationFlag = '-pthread -fopenmp';
else
    parallelizationFlag = '';
end


% additional flags
optimizationFlag = '-O3';

% Check if we are running octave
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;


%% Now we proceed to compilation
if useSharedGmpAndMpfr == 1
    if (ismac) || (~isunix)
        % If we pass here, we expect to be on a non-linux machine, but
        % environments are only well defined on linux... so we suggest to
        % go proceed another way
        warning('Trying to compile the GEM library with shared GMP and MPFR libraries, but it is not clear whether these libraries are available on the system. You may want to try compiling without shared libraries first.');
    end
    
    if isOctave
        flags(1) = unix(['g++ -c -Wdate-time  -I"', matlabroot, '/include/octave-', version, '/octave" -I/usr/include/hdf5/serial -I/usr/include/mpi -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', spectraFolder, '/include -I/usr/include -fstack-protector-strong -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG src/utils.cpp -o src/utils.o']);
        flags(2) = unix(['g++ -c -Wdate-time  -I"', matlabroot, '/include/octave-', version, '/octave" -I/usr/include/hdf5/serial -I/usr/include/mpi -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', spectraFolder, '/include -I/usr/include -fstack-protector-strong -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG src/gem.cpp -o src/gem.o']);
        flags(3) = unix(['g++ -c -Wdate-time  -I"', matlabroot, '/include/octave-', version, '/octave" -I/usr/include/hdf5/serial -I/usr/include/mpi -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', spectraFolder, '/include -I/usr/include -fstack-protector-strong -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG src/gem_mex.cpp -o src/gem_mex.o']);
        flags(4) = unix(['g++ -c -Wdate-time  -I"', matlabroot, '/include/octave-', version, '/octave" -I/usr/include/hdf5/serial -I/usr/include/mpi -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', spectraFolder, '/include -I/usr/include -fstack-protector-strong -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG src/sgem.cpp -o src/sgem.o']);
        flags(5) = unix(['g++ -c -Wdate-time  -I"', matlabroot, '/include/octave-', version, '/octave" -I/usr/include/hdf5/serial -I/usr/include/mpi -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', spectraFolder, '/include -I/usr/include -fstack-protector-strong -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG src/sgem_mex.cpp -o src/sgem_mex.o']);
        flags(6) = unix(['g++ -Wl,--no-undefined -fstack-protector-strong -Wl,-z,relro -shared src/gem_mex.o src/gem.o src/sgem.o src/utils.o -lmpfr -lgmp -L/usr/lib/x86_64-linux-gnu/octave/', version, ' -L/usr/lib/x86_64-linux-gnu -loctinterp -loctave -lstdc++ ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -o gem_mex.', mexext]);
        flags(7) = unix(['g++ -Wl,--no-undefined -fstack-protector-strong -Wl,-z,relro -shared src/sgem_mex.o src/sgem.o src/gem.o src/utils.o -lmpfr -lgmp -L/usr/lib/x86_64-linux-gnu/octave/', version, ' -L/usr/lib/x86_64-linux-gnu -loctinterp -loctave -lstdc++ ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -o sgem_mex.', mexext]);
    else
        flags(1) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE -DMEX_DOUBLE_HANDLE  -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG src/utils.cpp -o src/utils.o']);
        flags(2) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE -DMEX_DOUBLE_HANDLE  -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG src/gem.cpp -o src/gem.o']);
        flags(3) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE -DMEX_DOUBLE_HANDLE  -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG src/gem_mex.cpp -o src/gem_mex.o']);
        flags(4) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE -DMEX_DOUBLE_HANDLE  -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG src/sgem.cpp -o src/sgem.o']);
        flags(5) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE -DMEX_DOUBLE_HANDLE  -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG src/sgem_mex.cpp -o src/sgem_mex.o']);
        flags(6) = unix(['g++ -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/', lower(computer), '/mexFunction.map" src/gem_mex.o src/gem.o src/sgem.o src/utils.o   -lmpfr  -lgmp   -Wl,-rpath-link,', matlabroot, '/bin/', lower(computer), ' -L"', matlabroot, '/bin/', lower(computer), '" -lmx -lmex -lmat -lm -lstdc++ ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -o gem_mex.', mexext]);
        flags(7) = unix(['g++ -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/', lower(computer), '/mexFunction.map" src/sgem_mex.o src/sgem.o src/gem.o src/utils.o   -lmpfr  -lgmp   -Wl,-rpath-link,', matlabroot, '/bin/', lower(computer), ' -L"', matlabroot, '/bin/', lower(computer), '" -lmx -lmex -lmat -lm -lstdc++ ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -o sgem_mex.', mexext]);
    end
else
    % Here we compile a version of the library which does not rely on the
    % gmp and mpfr shared libraries. It requires additional source code :
    % - download the gmp source code from https://gmplib.org/#DOWNLOAD
    %   and unpack it into the 'external' folder (this creates the folder external/gmp*)
    % - download the mpfr source code from http://www.mpfr.org/mpfr-current/#download
    %   and unpack it into the 'external' folder (this creates the folder external/mpfr*)
    % - download the mpfrc++ source code from http://www.holoborodko.com/pavel/mpfr/
    %   and unpack it into the 'external' folder (this creates the folder external/mpfrc++*)

    % We check that the gmp, mpfr and mpfrc++ sources are available
    gmpFolder = '';
    mpfrFolder = '';
    mpfrcppFolder = '';
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
    end
    if isempty(gmpFolder)
        error('The gmp source code was not found. Please download it from https://gmplib.org and place it into the gem/external folder');
    end
    if isempty(mpfrFolder)
        error('The mpfr source code was not found. Please download it from http://www.mpfr.org and place it into the gem/external folder');
    end
    if isempty(mpfrcppFolder) || ~isequal(exist(['external/', mpfrcppFolder, '/mpreal.h']), 2)
        error('The mpfrc++ source code was not found. Please download it from http://www.holoborodko.com/pavel/mpfr and place it into the gem/external folder');
    end
    
    % We compile the libraries if needed
    cd external
    if exist('staticLibraries') ~= 7
        unix('mkdir staticLibraries');
    end
    if exist('staticLibraries/lib/libgmp.a') ~= 2
        if ispc
            error('Please download msys (e.g. from https://sourceforge.net/projects/mingw-w64/files/External%20binary%20packages%20%28Win64%20hosted%29/MSYS%20%2832-bit%29/MSYS-20111123.zip/download) and run the instructions in this section of the make.m file from within msys. Then, try to launch make.m again.');
            % The instructions to be run on msys are the 8 following 'unix' commands
            % (to be executed in their respective folders)
        end
        cd(gmpFolder);
        unix('./configure --disable-shared --enable-static CFLAGS=-fPIC --with-pic --prefix=`pwd`/../staticLibraries');
        unix('make');
        unix('make check');
        unix('make install');
        cd ..
    end
    if exist('staticLibraries/lib/libmpfr.a') ~= 2
        cd(mpfrFolder);
        unix('./configure --disable-shared --enable-static CFLAGS=-fPIC --with-pic --prefix=`pwd`/../staticLibraries --with-gmp=`pwd`/../staticLibraries');
        unix('make');
        unix('make check');
        unix('make install');
        cd ..
    end
    cd ..

    % Now gmp and mpfr should be ready, we can compile gem
    cd src
    % We adjust all paths
    eigenFolder = ['../', eigenFolder];
    spectraFolder = ['../', spectraFolder];
    gmpFolder = ['../external/', gmpFolder];
    mpfrFolder = ['../external/', mpfrFolder];
    mpfrcppFolder = ['../external/', mpfrcppFolder];
    if ismac
        % Octave on OSX returns true for both ismac and isunix, so we start
        % by the maxos case
        if isOctave
            flags = 0;
            try
                % We just do a generic build
                eval(['mex -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', spectraFolder, '/include -I', mpfrcppFolder, ' -I../external/staticLibraries/include -L../external/staticLibraries/lib -lmpfr -lgmp ', parallelizationFlag, ' ', optimizationFlag, ' gem_mex.cpp gem.cpp sgem.cpp utils.cpp']);
                eval(['mex -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', spectraFolder, '/include -I', mpfrcppFolder, ' -I../external/staticLibraries/include -L../external/staticLibraries/lib -lmpfr -lgmp ', parallelizationFlag, ' ', optimizationFlag, ' sgem_mex.cpp gem.cpp sgem.cpp utils.cpp']);
                unix('mv *.', mexext, ' ..');
            catch me
                flags = 1;
            end
        else
            flags(1) = unix(['/usr/bin/xcrun -sdk macosx10.12 clang++ -c -DMATLAB_MEX_FILE  -I../external/eigen -I../external/eigen/unsupported -Impfrc++-3 -I../external/staticLibraries/include -I../external/spectra/include -I/usr/include  -I"/Applications/MATLAB_R2016a.app/extern/include" -I"/Applications/MATLAB_R2016a.app/simulink/include" -fno-common -arch x86_64 -mmacosx-version-min=10.12 -fexceptions -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.12.sdk -fobjc-arc -std=c++11 -stdlib=libc++ -O2 -fwrapv -DNDEBUG utils.cpp -o utils.o']);
            flags(2) = unix(['/usr/bin/xcrun -sdk macosx10.12 clang++ -c -DMATLAB_MEX_FILE  -I../external/eigen -I../external/eigen/unsupported -Impfrc++-3 -I../external/staticLibraries/include -I../external/spectra/include -I/usr/include  -I"/Applications/MATLAB_R2016a.app/extern/include" -I"/Applications/MATLAB_R2016a.app/simulink/include" -fno-common -arch x86_64 -mmacosx-version-min=10.12 -fexceptions -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.12.sdk -fobjc-arc -std=c++11 -stdlib=libc++ -O2 -fwrapv -DNDEBUG gem.cpp -o gem.o']);
            flags(3) = unix(['/usr/bin/xcrun -sdk macosx10.12 clang++ -c -DMATLAB_MEX_FILE  -I../external/eigen -I../external/eigen/unsupported -Impfrc++-3 -I../external/staticLibraries/include -I../external/spectra/include -I/usr/include  -I"/Applications/MATLAB_R2016a.app/extern/include" -I"/Applications/MATLAB_R2016a.app/simulink/include" -fno-common -arch x86_64 -mmacosx-version-min=10.12 -fexceptions -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.12.sdk -fobjc-arc -std=c++11 -stdlib=libc++ -O2 -fwrapv -DNDEBUG gem_mex.cpp -o gem_mex.o']);
            flags(4) = unix(['/usr/bin/xcrun -sdk macosx10.12 clang++ -c -DMATLAB_MEX_FILE  -I../external/eigen -I../external/eigen/unsupported -Impfrc++-3 -I../external/staticLibraries/include -I../external/spectra/include -I/usr/include  -I"/Applications/MATLAB_R2016a.app/extern/include" -I"/Applications/MATLAB_R2016a.app/simulink/include" -fno-common -arch x86_64 -mmacosx-version-min=10.12 -fexceptions -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.12.sdk -fobjc-arc -std=c++11 -stdlib=libc++ -O2 -fwrapv -DNDEBUG sgem.cpp -o sgem.o']);
            flags(5) = unix(['/usr/bin/xcrun -sdk macosx10.12 clang++ -c -DMATLAB_MEX_FILE  -I../external/eigen -I../external/eigen/unsupported -Impfrc++-3 -I../external/staticLibraries/include -I../external/spectra/include -I/usr/include  -I"/Applications/MATLAB_R2016a.app/extern/include" -I"/Applications/MATLAB_R2016a.app/simulink/include" -fno-common -arch x86_64 -mmacosx-version-min=10.12 -fexceptions -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.12.sdk -fobjc-arc -std=c++11 -stdlib=libc++ -O2 -fwrapv -DNDEBUG sgem_mex.cpp -o sgem_mex.o']);
            flags(6) = unix(['/usr/bin/xcrun -sdk macosx10.12 clang++ -Wl,-twolevel_namespace -undefined error -arch x86_64 -mmacosx-version-min=10.12 -Wl,-syslibroot,/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.12.sdk -framework Cocoa -bundle  -Wl,-exported_symbols_list,"/Applications/MATLAB_R2016a.app/extern/lib/maci64/mexFunction.map" -stdlib=libc++ -O -Wl,-exported_symbols_list,"/Applications/MATLAB_R2016a.app/extern/lib/maci64/mexFunction.map" gem_mex.o gem.o sgem.o utils.o ../external/staticLibraries/lib/libmpfr.a  -lmpfr  -lgmp   -L"/Applications/MATLAB_R2016a.app/bin/maci64" -lmx -lmex -lmat -o ../gem_mex.mexmaci64']);
            flags(7) = unix(['/usr/bin/xcrun -sdk macosx10.12 clang++ -Wl,-twolevel_namespace -undefined error -arch x86_64 -mmacosx-version-min=10.12 -Wl,-syslibroot,/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.12.sdk -framework Cocoa -bundle  -Wl,-exported_symbols_list,"/Applications/MATLAB_R2016a.app/extern/lib/maci64/mexFunction.map" -stdlib=libc++ -O -Wl,-exported_symbols_list,"/Applications/MATLAB_R2016a.app/extern/lib/maci64/mexFunction.map" sgem_mex.o gem.o sgem.o utils.o ../external/staticLibraries/lib/libmpfr.a  -lmpfr  -lgmp   -L"/Applications/MATLAB_R2016a.app/bin/maci64" -lmx -lmex -lmat -o ../sgem_mex.mexmaci64']);
        end
    elseif isunix
        if isOctave
            flags(1) = unix(['g++ -c -Wdate-time  -I"', matlabroot, '/include/octave-', version, '/octave" -I/usr/include/hdf5/serial -I/usr/include/mpi -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', mpfrcppFolder, ' -I../external/staticLibraries/include -I', spectraFolder, '/include -I/usr/include -fstack-protector-strong -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG utils.cpp -o utils.o']);
            flags(2) = unix(['g++ -c -Wdate-time  -I"', matlabroot, '/include/octave-', version, '/octave" -I/usr/include/hdf5/serial -I/usr/include/mpi -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', mpfrcppFolder, ' -I../external/staticLibraries/include -I', spectraFolder, '/include -I/usr/include -fstack-protector-strong -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG gem.cpp -o gem.o']);
            flags(3) = unix(['g++ -c -Wdate-time  -I"', matlabroot, '/include/octave-', version, '/octave" -I/usr/include/hdf5/serial -I/usr/include/mpi -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', mpfrcppFolder, ' -I../external/staticLibraries/include -I', spectraFolder, '/include -I/usr/include -fstack-protector-strong -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG gem_mex.cpp -o gem_mex.o']);
            flags(4) = unix(['g++ -c -Wdate-time  -I"', matlabroot, '/include/octave-', version, '/octave" -I/usr/include/hdf5/serial -I/usr/include/mpi -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', mpfrcppFolder, ' -I../external/staticLibraries/include -I', spectraFolder, '/include -I/usr/include -fstack-protector-strong -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG sgem.cpp -o sgem.o']);
            flags(5) = unix(['g++ -c -Wdate-time  -I"', matlabroot, '/include/octave-', version, '/octave" -I/usr/include/hdf5/serial -I/usr/include/mpi -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', mpfrcppFolder, ' -I../external/staticLibraries/include -I', spectraFolder, '/include -I/usr/include -fstack-protector-strong -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG sgem_mex.cpp -o sgem_mex.o']);
            flags(6) = unix(['g++ -Wl,--no-undefined -fstack-protector-strong -Wl,-z,relro -shared gem_mex.o gem.o sgem.o utils.o ../external/staticLibraries/lib/libmpfr.a ../external/staticLibraries/lib/libgmp.a -L/usr/lib/x86_64-linux-gnu/octave/', version, ' -L/usr/lib/x86_64-linux-gnu -loctinterp -loctave -lstdc++ ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -o ../gem_mex.', mexext]);
            flags(7) = unix(['g++ -Wl,--no-undefined -fstack-protector-strong -Wl,-z,relro -shared sgem_mex.o sgem.o gem.o utils.o ../external/staticLibraries/lib/libmpfr.a ../external/staticLibraries/lib/libgmp.a -L/usr/lib/x86_64-linux-gnu/octave/', version, ' -L/usr/lib/x86_64-linux-gnu -loctinterp -loctave -lstdc++ ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -o ../sgem_mex.', mexext]);
        else
            flags(1) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE -DMEX_DOUBLE_HANDLE  -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', mpfrcppFolder, ' -I../external/staticLibraries/include -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG utils.cpp -o utils.o']);
            flags(2) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE -DMEX_DOUBLE_HANDLE  -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', mpfrcppFolder, ' -I../external/staticLibraries/include -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG gem.cpp -o gem.o']);
            flags(3) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE -DMEX_DOUBLE_HANDLE  -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', mpfrcppFolder, ' -I../external/staticLibraries/include -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG gem_mex.cpp -o gem_mex.o']);
            flags(4) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE -DMEX_DOUBLE_HANDLE  -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', mpfrcppFolder, ' -I../external/staticLibraries/include -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG sgem.cpp -o sgem.o']);
            flags(5) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE -DMEX_DOUBLE_HANDLE  -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', mpfrcppFolder, ' -I../external/staticLibraries/include -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG sgem_mex.cpp -o sgem_mex.o']);
            flags(6) = unix(['g++ -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/', lower(computer), '/mexFunction.map" gem_mex.o gem.o sgem.o utils.o  ../external/staticLibraries/lib/libmpfr.a ../external/staticLibraries/lib/libgmp.a  -Wl,-rpath-link,', matlabroot, '/bin/', lower(computer), ' -L"', matlabroot, '/bin/', lower(computer), '" -lmx -lmex -lmat -lm -lstdc++ ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -o ../gem_mex.', mexext]);
            flags(7) = unix(['g++ -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/', lower(computer), '/mexFunction.map" sgem_mex.o sgem.o gem.o utils.o  ../external/staticLibraries/lib/libmpfr.a ../external/staticLibraries/lib/libgmp.a  -Wl,-rpath-link,', matlabroot, '/bin/', lower(computer), ' -L"', matlabroot, '/bin/', lower(computer), '" -lmx -lmex -lmat -lm -lstdc++ ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -o ../sgem_mex.', mexext]);
        end
    elseif ispc
        if isOctave
            flags = 0;
            try
                % We try to do a generic build
                eval(['mex -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', spectraFolder, '/include -I', mpfrcppFolder, ' -I../external/staticLibraries/include -L../external/staticLibraries/lib -lmpfr -lgmp ', parallelizationFlag, ' ', optimizationFlag, ' gem_mex.cpp gem.cpp sgem.cpp utils.cpp']);
                eval(['mex -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', spectraFolder, '/include -I', mpfrcppFolder, ' -I../external/staticLibraries/include -L../external/staticLibraries/lib -lmpfr -lgmp ', parallelizationFlag, ' ', optimizationFlag, ' sgem_mex.cpp gem.cpp sgem.cpp utils.cpp']);
                unix('mv *.', mexext, ' ..');
            catch me
                flags = 1;
            end
        else
            flags(1) = unix(['g++ -c -m64 -DMATLAB_MEX_FILE  -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', mpfrcppFolder, ' -I../external/staticLibraries/include -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG utils.cpp -o utils.o']);
            flags(2) = unix(['g++ -c -m64 -DMATLAB_MEX_FILE  -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', mpfrcppFolder, ' -I../external/staticLibraries/include -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG gem.cpp -o gem.o']);
            flags(3) = unix(['g++ -c -m64 -DMATLAB_MEX_FILE  -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', mpfrcppFolder, ' -I../external/staticLibraries/include -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG gem_mex.cpp -o gem_mex.o']);
            flags(4) = unix(['g++ -c -m64 -DMATLAB_MEX_FILE  -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', mpfrcppFolder, ' -I../external/staticLibraries/include -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG sgem.cpp -o sgem.o']);
            flags(5) = unix(['g++ -c -m64 -DMATLAB_MEX_FILE  -I', eigenFolder, ' -I', eigenFolder, '/unsupported -I', mpfrcppFolder, ' -I../external/staticLibraries/include -I', spectraFolder, '/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -DNDEBUG sgem_mex.cpp -o sgem_mex.o']);
            flags(6) = unix(['g++ -m64 -s -Wl,--no-undefined  -shared gem_mex.o gem.o sgem.o utils.o  ../external/staticLibraries/lib/libmpfr.a ../external/staticLibraries/lib/libgmp.a -L"', matlabroot, '/bin/', lower(computer), '" -llibmx -llibmex -llibmat -lm -llibmwlapack -llibmwblas ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -L"', matlabroot, '\extern\lib\win64\mingw64" -o ../gem_mex.', mexext]);
            flags(7) = unix(['g++ -m64 -s -Wl,--no-undefined  -shared sgem_mex.o sgem.o gem.o utils.o  ../external/staticLibraries/lib/libmpfr.a ../external/staticLibraries/lib/libgmp.a -L"', matlabroot, '/bin/', lower(computer), '" -llibmx -llibmex -llibmat -lm -llibmwlapack -llibmwblas ', parallelizationFlag, ' ', optimizationFlag, ' -DEIGEN_NO_DEBUG -L"', matlabroot, '\extern\lib\win64\mingw64" -o ../sgem_mex.', mexext]);
        end
    end
    cd ..
end

% Was the compilation successful?
if sum(flags) == 0
    disp('Compilation successful.');
else
    error('Compilation error.');
end

return;

% The following line displays the default instructions for compilation 
% with default options (and no parallelization) :
%mex -n -largeArrayDims -Isrc/eigen -Isrc/eigen/unsupported -Isrc/spectra/include -I/usr/include -lmpfr -lgmp src/gem_mex.cpp src/gem.cpp src/sgem.cpp src/utils.cpp

