function make
% Run this file from matlab to compile C++ part of the GEM library
% Compilation :
%   mex -I../eigen -I../eigen/unsupported -I../eigen/unsupported/test/mpreal -I/usr/include gem_mex.cpp -lmpfr -lgmp
% after running
%   sudo apt-get install libmpfrc++-dev libgmp-dev

v = version;
vDots = find(v=='.');
v = str2num(v(1:vDots(2)-1));

if v >= 9
    % Version 2016a without parallelization
%    /usr/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"/opt/MATLAB/R2016a/extern/include" -I"/opt/MATLAB/R2016a/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -O -DNDEBUG /home/bancal/Sync/MangoAssembled/Matlab_pour_lrs/Eigen/gem/src/gem_mex.cpp -o /tmp/mex_5920715111474_21812/gem_mex.o
%    /usr/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"/opt/MATLAB/R2016a/extern/include" -I"/opt/MATLAB/R2016a/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -O -DNDEBUG /home/bancal/Sync/MangoAssembled/Matlab_pour_lrs/Eigen/gem/src/gem.cpp -o /tmp/mex_5920715111474_21812/gem.o
%    /usr/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"/opt/MATLAB/R2016a/extern/include" -I"/opt/MATLAB/R2016a/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -O -DNDEBUG /home/bancal/Sync/MangoAssembled/Matlab_pour_lrs/Eigen/gem/src/utils.cpp -o /tmp/mex_5920715111474_21812/utils.o
%    /usr/bin/g++ -pthread -Wl,--no-undefined  -shared -O -Wl,--version-script,"/opt/MATLAB/R2016a/extern/lib/glnxa64/mexFunction.map" /tmp/mex_5920715111474_21812/gem_mex.o /tmp/mex_5920715111474_21812/gem.o /tmp/mex_5920715111474_21812/utils.o   -lmpfr  -lgmp   -Wl,-rpath-link,/opt/MATLAB/R2016a/bin/glnxa64 -L"/opt/MATLAB/R2016a/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -o gem_mex.mexa64
	flags(1) = unix('/usr/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"/opt/MATLAB/R2016a/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG src/utils.cpp -o utils.o');
	flags(2) = unix('/usr/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"/opt/MATLAB/R2016a/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG src/gem.cpp -o gem.o');
	flags(3) = unix('/usr/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"/opt/MATLAB/R2016a/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG src/gem_mex.cpp -o gem_mex.o');
	flags(4) = unix('/usr/bin/g++ -Wl,--no-undefined  -shared -Wl,--version-script,"/opt/MATLAB/R2016a/extern/lib/glnxa64/mexFunction.map" gem_mex.o gem.o sgem.o utils.o   -lmpfr  -lgmp   -Wl,-rpath-link,/opt/MATLAB/R2016a/bin/glnxa64 -L"/opt/MATLAB/R2016a/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -O3 -DEIGEN_NO_DEBUG -o gem_mex.mexa64');
	flags(5) = unix('/usr/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"/opt/MATLAB/R2016a/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG src/sgem.cpp -o sgem.o');
	flags(6) = unix('/usr/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"/opt/MATLAB/R2016a/extern/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG src/sgem_mex.cpp -o sgem_mex.o');
	flags(7) = unix('/usr/bin/g++ -Wl,--no-undefined  -shared -Wl,--version-script,"/opt/MATLAB/R2016a/extern/lib/glnxa64/mexFunction.map" sgem_mex.o sgem.o gem.o utils.o   -lmpfr  -lgmp   -Wl,-rpath-link,/opt/MATLAB/R2016a/bin/glnxa64 -L"/opt/MATLAB/R2016a/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -O3 -DEIGEN_NO_DEBUG -o sgem_mex.mexa64');
    if sum(flags) == 0
        display('Compilation successful.');
    else
        error('Compilation error.');
    end
elseif v >= 9
    % Version 2016a with parallelization
	flags(1) = unix('/usr/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"/opt/MATLAB/R2016a/extern/include" -I"/opt/MATLAB/R2016a/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/utils.cpp -o utils.o');
	flags(2) = unix('/usr/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"/opt/MATLAB/R2016a/extern/include" -I"/opt/MATLAB/R2016a/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/gem.cpp -o gem.o');
	flags(3) = unix('/usr/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"/opt/MATLAB/R2016a/extern/include" -I"/opt/MATLAB/R2016a/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/gem_mex.cpp -o gem_mex.o');
	flags(4) = unix('/usr/bin/g++ -pthread -Wl,--no-undefined  -shared -Wl,--version-script,"/opt/MATLAB/R2016a/extern/lib/glnxa64/mexFunction.map" gem_mex.o gem.o sgem.o utils.o   -lmpfr  -lgmp   -Wl,-rpath-link,/opt/MATLAB/R2016a/bin/glnxa64 -L"/opt/MATLAB/R2016a/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -O3 -DEIGEN_NO_DEBUG -fopenmp -o gem_mex.mexa64');
	flags(5) = unix('/usr/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"/opt/MATLAB/R2016a/extern/include" -I"/opt/MATLAB/R2016a/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/sgem.cpp -o sgem.o');
	flags(6) = unix('/usr/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"/opt/MATLAB/R2016a/extern/include" -I"/opt/MATLAB/R2016a/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/sgem_mex.cpp -o sgem_mex.o');
	flags(7) = unix('/usr/bin/g++ -pthread -Wl,--no-undefined  -shared -Wl,--version-script,"/opt/MATLAB/R2016a/extern/lib/glnxa64/mexFunction.map" sgem_mex.o sgem.o gem.o utils.o   -lmpfr  -lgmp   -Wl,-rpath-link,/opt/MATLAB/R2016a/bin/glnxa64 -L"/opt/MATLAB/R2016a/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -O3 -DEIGEN_NO_DEBUG -fopenmp -o sgem_mex.mexa64');
    if sum(flags) == 0
        display('Compilation successful.');
    else
        error('Compilation error.');
    end
elseif v >= 7.12
    % Version 2015b with parallelization
    % Compile -O3 -DEIGEN_NO_DEBUG and -fopenmp options, takes 14 seconds instead of 10 seconds
	flags(1) = unix('/usr/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"/opt/MATLAB/R2015b/extern/include" -I"/opt/MATLAB/R2015b/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/utils.cpp -o utils.o');
	flags(2) = unix('/usr/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"/opt/MATLAB/R2015b/extern/include" -I"/opt/MATLAB/R2015b/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/gem.cpp -o gem.o');
	flags(3) = unix('/usr/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"/opt/MATLAB/R2015b/extern/include" -I"/opt/MATLAB/R2015b/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/gem_mex.cpp -o gem_mex.o');
	%flags(4) = unix('/usr/bin/g++ -pthread -Wl,--no-undefined  -shared -Wl,--version-script,"/opt/MATLAB/R2015b/extern/lib/glnxa64/mexFunction.map" gem_mex.o gem.o sgem.o utils.o /usr/lib/x86_64-linux-gnu/libmpfr.a /usr/lib/x86_64-linux-gnu/libgmp.a -Wl,-rpath-link,/opt/MATLAB/R2015b/bin/glnxa64 -L"/opt/MATLAB/R2015b/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -O3 -DEIGEN_NO_DEBUG -fopenmp -o gem_mex.mexa64');
    flags(4) = unix('/usr/bin/g++ -pthread -Wl,--no-undefined  -shared -Wl,--version-script,"/opt/MATLAB/R2015b/extern/lib/glnxa64/mexFunction.map" gem_mex.o gem.o sgem.o utils.o  -lmpfr  -lgmp   -Wl,-rpath-link,/opt/MATLAB/R2015b/bin/glnxa64 -L"/opt/MATLAB/R2015b/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -O3 -DEIGEN_NO_DEBUG -fopenmp -o gem_mex.mexa64');
	flags(5) = unix('/usr/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"/opt/MATLAB/R2015b/extern/include" -I"/opt/MATLAB/R2015b/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/sgem.cpp -o sgem.o');
	flags(6) = unix('/usr/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include  -I"/opt/MATLAB/R2015b/extern/include" -I"/opt/MATLAB/R2015b/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -Wno-deprecated -pthread -std=c++11 -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/sgem_mex.cpp -o sgem_mex.o');
	flags(7) = unix('/usr/bin/g++ -pthread -Wl,--no-undefined  -shared -Wl,--version-script,"/opt/MATLAB/R2015b/extern/lib/glnxa64/mexFunction.map" sgem_mex.o sgem.o gem.o utils.o   -lmpfr  -lgmp   -Wl,-rpath-link,/opt/MATLAB/R2015b/bin/glnxa64 -L"/opt/MATLAB/R2015b/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -O3 -DEIGEN_NO_DEBUG -fopenmp -o sgem_mex.mexa64');
    if sum(flags) == 0
        display('Compilation successful.');
    else
        error('Compilation error.');
    end
elseif v >= 7.11
    % Version 2010b with parallelization
    flags(1) = unix('g++ -c  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include -I/opt/MATLABR2010b/extern/include -I/opt/MATLABR2010b/simulink/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -pthread -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/utils.cpp');
    flags(2) = unix('g++ -c  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include -I/opt/MATLABR2010b/extern/include -I/opt/MATLABR2010b/simulink/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -pthread -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/gem.cpp');
    flags(3) = unix('g++ -c  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include -I/opt/MATLABR2010b/extern/include -I/opt/MATLABR2010b/simulink/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -pthread -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/gem_mex.cpp');
    flags(4) = unix('g++ -pthread -shared -Wl,--version-script,/opt/MATLABR2010b/extern/lib/glnxa64/mexFunction.map -Wl,--no-undefined gem_mex.o gem.o sgem.o utils.o  -lmpfr -lgmp -Wl,-rpath-link,/opt/MATLABR2010b/bin/glnxa64 -L/opt/MATLABR2010b/bin/glnxa64 -lmx -lmex -lmat -lm -O3 -DEIGEN_NO_DEBUG -fopenmp -o gem_mex.mexa64');
    flags(5) = unix('g++ -c  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include -I/opt/MATLABR2010b/extern/include -I/opt/MATLABR2010b/simulink/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -pthread -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/sgem.cpp');
    flags(6) = unix('g++ -c  -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include -I/opt/MATLABR2010b/extern/include -I/opt/MATLABR2010b/simulink/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -pthread -O3 -DEIGEN_NO_DEBUG -DNDEBUG -fopenmp src/sgem_mex.cpp');
    flags(7) = unix('g++ -pthread -shared -Wl,--version-script,/opt/MATLABR2010b/extern/lib/glnxa64/mexFunction.map -Wl,--no-undefined sgem_mex.o sgem.o gem.o utils.o  -lmpfr -lgmp -Wl,-rpath-link,/opt/MATLABR2010b/bin/glnxa64 -L/opt/MATLABR2010b/bin/glnxa64 -lmx -lmex -lmat -lm -O3 -DEIGEN_NO_DEBUG -fopenmp -o sgem_mex.mexa64');
    if sum(flags) == 0
        display('Compilation successful.');
    else
        error('Compilation error.');
    end
else
    % Compile with default options
    mex -largeArrayDims -Isrc/eigen -Isrc/eigen/unsupported -Ieigen -Ieigen/unsupported -I/usr/include -lmpfr -lgmp src/gem_mex.cpp src/gem.cpp src/utils.cpp
end

return;

% Here is a step by step version : 
mex -I../eigen -I../eigen/unsupported -I/usr/include -c utils.cpp
mex -I../eigen -I../eigen/unsupported -I/usr/include -c gem.cpp
mex -I../eigen -I../eigen/unsupported -I/usr/include -c gem_mex.cpp
mex -I../eigen -I../eigen/unsupported -I/usr/include -lmpfr -lgmp gem_mex.o utils.o gem.o
