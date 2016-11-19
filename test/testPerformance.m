% Here we perform some performance test

% First we choose the precision
nbDigits = 100;
digits(nbDigits);
gemWorkingPrecision(nbDigits);

% Now we define some random matrices
size = 100;
a1 = rand(size,size);
a2 = rand(size,size);

% We test the conversion speed from double to high precision
tic;b1 = vpa(a1);toc
tic;b2 = vpa(a2);toc
tic;c1 = gem(a1);toc
tic;c2 = gem(a2);toc

% We test the conversion speed from high precision to double
tic;double(b1);toc
tic;double(b2);toc
tic;double(c1);toc
tic;double(c2);toc

% We test the minimization function
tic;min(b1);toc
tic;min(b2);toc
tic;min(c1);toc
tic;min(c2);toc

% We test the matrix product
tic; b1*b2; toc
tic; b2*b1; toc
tic; c1*c2; toc
tic; c2*c1; toc
tic; a1*a2; toc
tic; a2*a1; toc
