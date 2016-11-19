% Here we test the performance of matrix products

% We vary the precision of one of the matrices
A = gem(rand(500));
co = 0;
precisions = [1 2 3 4 5 7 9 10 20 30 40 50 60 70 80 90 100];
for precision = precisions
    co = co + 1
    tmp = [];
    for i = 1:10
        B = gem(round(rand(500)), precision);
        tic;
        C = A*B;
        tmp(i) = toc;
    end
    timing(co) = mean(tmp);
end

