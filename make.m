function make(parallelization)
% make(parallelization)
%
% Run this file from matlab to compile the C++ part of the GEM library. Use
% parallelization = 1 to enable openmp, parallelization = 0 otherwise.
%
% Note : the code relies on the gmp, mpfrc++ and eigen libraries. The Eigen
% library should be dowloaded from http://eigen.tuxfamily.org and placed in
% the gem/src folder. On ubuntu, the two other libraries can be installed 
% by running the following command in the terminal :
%   sudo apt-get install libmpfrc++-dev libgmp-dev
%
% It is also possible to perform a compilation with gmp and mpfr source 
% codes instead of relying on precompiled gmp and mpfr libraries. In this
% case proceed as follow :
% - Download the gmp source code should from https://gmplib.org/#DOWNLOAD
%   and unpack it into the 'src' folder (this creates the folder src/gmp*)
% - Download the mpfr source code from http://www.mpfr.org/mpfr-current/#download
%   and unpack it into the 'src' folder
% - Set the value of 'useSharedGmpAndMpfr' to 0 in this file
% - Run this file from matlab

v = version;
vDots = find(v=='.');
v = str2num(v(1:vDots(2)-1));

if nargin < 1
    parallelization = 1;
end

% Set the following variable to 0 in order to compile a library which does
% not depend on shared gmp and mpfr libraries:
useSharedGmpAndMpfr = 0;

% We check that the eigen source code was downloaded
list = dir('src');
eigenFound = 0;
for i=1:length(list)
    eigenFound = eigenFound + isequal(list(i).name,'eigen');
end
if eigenFound == 0
    error('The eigen library was not found in the gem/src folder. Please download it from http://eigen.tuxfamily.org');
end

if useSharedGmpAndMpfr == 1
    if parallelization == 1
        % Compilation with parallelization
        cd src
        flags(1) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp utils.cpp -o utils.o']);
        flags(2) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp gem.cpp -o gem.o']);
        flags(3) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp gem_mex.cpp -o gem_mex.o']);
        flags(4) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp sgem.cpp -o sgem.o']);
        flags(5) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp sgem_mex.cpp -o sgem_mex.o']);
        flags(6) = unix(['g++ -pthread -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/glnxa64/mexFunction.map" gem_mex.o gem.o sgem.o utils.o   -lmpfr  -lgmp   -Wl,-rpath-link,', matlabroot, '/bin/glnxa64 -L"', matlabroot, '/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -O3 -DEIGEN_NO_DEBUG -fopenmp -o ../gem_mex.', mexext]);
        flags(7) = unix(['g++ -pthread -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/glnxa64/mexFunction.map" sgem_mex.o sgem.o gem.o utils.o   -lmpfr  -lgmp   -Wl,-rpath-link,', matlabroot, '/bin/glnxa64 -L"', matlabroot, '/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -O3 -DEIGEN_NO_DEBUG -fopenmp -o ../sgem_mex.', mexext]);
        cd ..
        if sum(flags) == 0
            display('Compilation successful.');
        else
            error('Compilation error.');
        end
    else
        % Compilation without parallelization
        cd src
        flags(1) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -I/usr/include  -I"', matlabroot, '/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG utils.cpp -o utils.o']);
        flags(2) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -I/usr/include  -I"', matlabroot, '/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG gem.cpp -o gem.o']);
        flags(3) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -I/usr/include  -I"', matlabroot, '/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG gem_mex.cpp -o gem_mex.o']);
        flags(4) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -I/usr/include  -I"', matlabroot, '/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG sgem.cpp -o sgem.o']);
        flags(5) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -I/usr/include  -I"', matlabroot, '/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG sgem_mex.cpp -o sgem_mex.o']);
        flags(6) = unix(['g++ -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/glnxa64/mexFunction.map" gem_mex.o gem.o sgem.o utils.o   -lmpfr  -lgmp   -Wl,-rpath-link,', matlabroot, '/bin/glnxa64 -L"', matlabroot, '/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -O3 -DEIGEN_NO_DEBUG -o ../gem_mex.', mexext]);
        flags(7) = unix(['g++ -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/glnxa64/mexFunction.map" sgem_mex.o sgem.o gem.o utils.o   -lmpfr  -lgmp   -Wl,-rpath-link,', matlabroot, '/bin/glnxa64 -L"', matlabroot, '/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -O3 -DEIGEN_NO_DEBUG -o ../sgem_mex.', mexext]);
        cd ..
        if sum(flags) == 0
            display('Compilation successful.');
        else
            error('Compilation error.');
        end
    end
else
    % Here we compile a version of the library which does not rely on gmp 
    % and mpfr shared libraries. It requires additional source code :
    % - download the gmp source code from https://gmplib.org/#DOWNLOAD
    %   and unpack it into the 'src' folder (this creates the folder src/gmp*)
    % - download the mpfr source code from http://www.mpfr.org/mpfr-current/#download
    %   and unpack it into the 'src' folder

    % We check that the gmp and mpfr sources have been downloaded
    gmpFolder = '';
    mpfrFolder = '';
    for i=1:length(list)
        if isequal(list(i).name(1:min(end,3)),'gmp')
            gmpFolder = list(i).name;
        end
        if isequal(list(i).name(1:min(end,4)),'mpfr')
            mpfrFolder = list(i).name;
        end
    end
    if isempty(gmpFolder)
        error('The gmp source code was not found. Please download it from https://gmplib.org');
    end
    if isempty(mpfrFolder)
        error('The mpfr source code was not found. Please download it from http://www.mpfr.org');
    end
    
    % we compile the libraries if needed
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
    if parallelization == 1
        % Compilation with parallelization
        flags(1) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -IstaticLibraries/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp utils.cpp -o utils.o']);
        flags(2) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -IstaticLibraries/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp gem.cpp -o gem.o']);
        flags(3) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -IstaticLibraries/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp gem_mex.cpp -o gem_mex.o']);
        flags(4) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -IstaticLibraries/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp sgem.cpp -o sgem.o']);
        flags(5) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -IstaticLibraries/include -I/usr/include  -I"', matlabroot, '/extern/include" -I"', matlabroot, '/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp sgem_mex.cpp -o sgem_mex.o']);
        flags(6) = unix(['g++ -pthread -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/glnxa64/mexFunction.map" gem_mex.o gem.o sgem.o utils.o  staticLibraries/lib/libmpfr.a staticLibraries/lib/libgmp.a  -Wl,-rpath-link,', matlabroot, '/bin/glnxa64 -L"', matlabroot, '/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -O3 -DEIGEN_NO_DEBUG -fopenmp -o ../gem_mex.', mexext]);
        flags(7) = unix(['g++ -pthread -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/glnxa64/mexFunction.map" sgem_mex.o sgem.o gem.o utils.o  staticLibraries/lib/libmpfr.a staticLibraries/lib/libgmp.a  -Wl,-rpath-link,', matlabroot, '/bin/glnxa64 -L"', matlabroot, '/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -O3 -DEIGEN_NO_DEBUG -fopenmp -o ../sgem_mex.', mexext]);
        cd ..
        if sum(flags) == 0
            display('Compilation successful.');
        else
            error('Compilation error.');
        end
    else
        % Compilation without parallelization
        flags(1) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -IstaticLibraries/include -I/usr/include  -I"', matlabroot, '/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG utils.cpp -o utils.o']);
        flags(2) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -IstaticLibraries/include -I/usr/include  -I"', matlabroot, '/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG gem.cpp -o gem.o']);
        flags(3) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -IstaticLibraries/include -I/usr/include  -I"', matlabroot, '/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG gem_mex.cpp -o gem_mex.o']);
        flags(4) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -IstaticLibraries/include -I/usr/include  -I"', matlabroot, '/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG sgem.cpp -o sgem.o']);
        flags(5) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Ieigen -Ieigen/unsupported -Ieigen/unsupported/test/mpreal -IstaticLibraries/include -I/usr/include  -I"', matlabroot, '/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG sgem_mex.cpp -o sgem_mex.o']);
        flags(6) = unix(['g++ -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/glnxa64/mexFunction.map" gem_mex.o gem.o sgem.o utils.o  staticLibraries/lib/libmpfr.a staticLibraries/lib/libgmp.a  -Wl,-rpath-link,', matlabroot, '/bin/glnxa64 -L"', matlabroot, '/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -O3 -DEIGEN_NO_DEBUG -o ../gem_mex.', mexext]);
        flags(7) = unix(['g++ -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/glnxa64/mexFunction.map" sgem_mex.o sgem.o gem.o utils.o  staticLibraries/lib/libmpfr.a staticLibraries/lib/libgmp.a  -Wl,-rpath-link,', matlabroot, '/bin/glnxa64 -L"', matlabroot, '/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -O3 -DEIGEN_NO_DEBUG -o ../sgem_mex.', mexext]);
        cd ..
        if sum(flags) == 0
            display('Compilation successful.');
        else
            error('Compilation error.');
        end
    end

end

return;

% The following line launches the compilation with default options (and no
% parallelization) :
%mex -largeArrayDims -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include -lmpfr -lgmp src/gem_mex.cpp src/gem.cpp src/utils.cpp

