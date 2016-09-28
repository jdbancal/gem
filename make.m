function make(parallelization)
% make(parallelization)
%
% Run this file from matlab to compile the C++ part of the GEM library. Use
% parallelization = 1 to enable openmp, parallelization = 0 otherwise.
%
% Note : the code relies on the gmp, mpfrc++ and eigen libraries. The Eigen
% library should be dowloaded from http://eigen.tuxfamily.org and placed in
% the gem or gem/src folder. On ubuntu, the two other libraries can be 
% installed by running the following command :
%   sudo apt-get install libmpfrc++-dev libgmp-dev

v = version;
vDots = find(v=='.');
v = str2num(v(1:vDots(2)-1));

if nargin < 0
    parallelization = 1;
end

if parallelization == 1
    % Compilation with parallelization
	flags(1) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"', matlabroot, '/extern/include" -I"/opt/MATLAB/R2016a/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/utils.cpp -o utils.o']);
	flags(2) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"', matlabroot, '/extern/include" -I"/opt/MATLAB/R2016a/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/gem.cpp -o gem.o']);
	flags(3) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"', matlabroot, '/extern/include" -I"/opt/MATLAB/R2016a/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/gem_mex.cpp -o gem_mex.o']);
	flags(4) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"', matlabroot, '/extern/include" -I"/opt/MATLAB/R2016a/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/sgem.cpp -o sgem.o']);
	flags(5) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"', matlabroot, '/extern/include" -I"/opt/MATLAB/R2016a/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/sgem_mex.cpp -o sgem_mex.o']);
	flags(6) = unix(['g++ -pthread -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/glnxa64/mexFunction.map" gem_mex.o gem.o sgem.o utils.o   -lmpfr  -lgmp   -Wl,-rpath-link,', matlabroot, '/bin/glnxa64 -L"', matlabroot, '/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -O3 -DEIGEN_NO_DEBUG -fopenmp -o gem_mex.', mexext]);
	flags(7) = unix(['g++ -pthread -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/glnxa64/mexFunction.map" sgem_mex.o sgem.o gem.o utils.o   -lmpfr  -lgmp   -Wl,-rpath-link,', matlabroot, '/bin/glnxa64 -L"', matlabroot, '/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -O3 -DEIGEN_NO_DEBUG -fopenmp -o sgem_mex.', mexext]);
    if sum(flags) == 0
        display('Compilation successful.');
    else
        error('Compilation error.');
    end
else
    % Compilation without parallelization
	flags(1) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"', matlabroot, '/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG src/utils.cpp -o utils.o']);
	flags(2) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"', matlabroot, '/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG src/gem.cpp -o gem.o']);
	flags(3) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"', matlabroot, '/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG src/gem_mex.cpp -o gem_mex.o']);
	flags(4) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"', matlabroot, '/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG src/sgem.cpp -o sgem.o']);
	flags(5) = unix(['g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"', matlabroot, '/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG src/sgem_mex.cpp -o sgem_mex.o']);
	flags(6) = unix(['g++ -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/glnxa64/mexFunction.map" gem_mex.o gem.o sgem.o utils.o   -lmpfr  -lgmp   -Wl,-rpath-link,', matlabroot, '/bin/glnxa64 -L"', matlabroot, '/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -O3 -DEIGEN_NO_DEBUG -o gem_mex.', mexext]);
	flags(7) = unix(['g++ -Wl,--no-undefined  -shared -Wl,--version-script,"', matlabroot, '/extern/lib/glnxa64/mexFunction.map" sgem_mex.o sgem.o gem.o utils.o   -lmpfr  -lgmp   -Wl,-rpath-link,', matlabroot, '/bin/glnxa64 -L"', matlabroot, '/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -O3 -DEIGEN_NO_DEBUG -o sgem_mex.', mexext]);
    if sum(flags) == 0
        display('Compilation successful.');
    else
        error('Compilation error.');
    end
else
    % This line launches the compilation with default options
    mex -largeArrayDims -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include -lmpfr -lgmp src/gem_mex.cpp src/gem.cpp src/utils.cpp
end

