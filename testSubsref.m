% We test various subsref strategies


A = gemRand(501,1001);
A(find(A < 1/10)) = 0;
B = sgem(A);

global test
co1 = 0;
t1 = [];
t2 = [];
for i = 1:10:501
    co1 = co1 + 1;
    co2 = 0;
    cont = 1;
    for j = 1:10:501
        if cont == 1
            co2 = co2 + 1;
            test = 0;
            tic;
            tmp=B([1:i],[1:j]);
            t1(co1,co2) = toc;
            test = 1;
            tic;
            tmp=B([1:i],[1:j]);
            t2(co1,co2) = toc;
            if t2(co1,co2) < t1(co1,co2)
                cont = 0;
            end
        end
    end;
end;


A = gemRand(1001,501);
A(find(A < 1/10)) = 0;
B = sgem(A);

global test
co1 = 0;
t3 = [];
t4 = [];
for i = 1:10:501
    co1 = co1 + 1;
    co2 = 0;
    cont = 1;
    for j = 1:10:501
        if cont == 1
            co2 = co2 + 1;
            test = 0;
            tic;
            tmp=B([1:i],[1:j]);
            t3(co1,co2) = toc;
            test = 1;
            tic;
            tmp=B([1:i],[1:j]);
            t4(co1,co2) = toc;
            if t4(co1,co2) < t3(co1,co2)
                cont = 0;
            end
        end
    end;
end;

A = gemRand(501,1001);
A(find(A < 1/10)) = 0;
B = sgem(A);

global test
co1 = 0;
t5 = [];
t6 = [];
for i = 1:10:501
    co1 = co1 + 1;
    co2 = 0;
    cont = 1;
    for j = 1:10:501
        if cont == 1
            co2 = co2 + 1;
            test = 0;
            ri = randperm(i);
            rj = randperm(j);
            tic;
            tmp=B(ri,rj);
            t5(co1,co2) = toc;
            test = 1;
            tic;
            tmp=B(ri,rj);
            t6(co1,co2) = toc;
            if t6(co1,co2) < t5(co1,co2)
                cont = 0;
            end
        end
    end;
end;


A = gemRand(501,1001);
A(find(A < 1/100)) = 0;
B = sgem(A);

global test
co1 = 0;
t7 = [];
t8 = [];
for i = 1:10:501
    co1 = co1 + 1;
    co2 = 0;
    cont = 1;
    for j = 1:10:501
        if cont == 1
            co2 = co2 + 1;
            test = 0;
            tic;
            tmp=B([1:i],[1:j]);
            t7(co1,co2) = toc;
            test = 1;
            tic;
            tmp=B([1:i],[1:j]);
            t8(co1,co2) = toc;
            if t8(co1,co2) < t7(co1,co2)
                cont = 0;
            end
        end
    end;
end;
