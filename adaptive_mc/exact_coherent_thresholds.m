%% Compare SVT and our matrix completion algorithm on matrices with
%% coherent row space. 

clear all;
r = 5;
ps = 0.05:0.05:0.9;
iters = 100;
n = 100;

a1 = [];
for p=ps,
    s = 0;
    for i=1:iters,
        X = coherent_matrix(n,n,r,1);
        [xhat, obs] = mc_exact_adapt(X, p);
        if norm(xhat - X, 'fro') < 0.001,
            s = s+1;
        end;
    end;
    fprintf('p=%0.2f s=%0.2f\n', p, s/iters);
    a1 = [a1 s/iters];
end;

a2 = [];
for p=ps,
    s = 0;
    for i=1:iters,
        X = coherent_matrix(n,n,r,2);
        [xhat, obs] = mc_exact_adapt(X, p);
        if norm(xhat - X, 'fro') < 0.001,
            s = s+1;
        end;
    end;
    fprintf('p=%0.2f s=%0.2f\n', p, s/iters);
    a2 = [a2 s/iters];
end;


s1 = [];
for p=ps,
    s = 0;
    for i=1:iters,
        X = coherent_matrix(n,n,r,1);
        [xhat, obs] = mc_svt(X, p);
        if norm(xhat - X, 'fro') < 0.001,
            s = s+1;
        end;
    end;
    fprintf('p=%0.2f s=%0.2f\n', p, s/iters);
    s1 = [s1 s/iters];
end;

s2 = [];
for p=ps,
    s = 0;
    for i=1:iters,
        X = coherent_matrix(n,n,r,2);
        [xhat, obs] = mc_svt(X, p);
        if norm(xhat - X, 'fro') < 0.001,
            s = s+1;
        end;
    end;
    fprintf('p=%0.2f s=%0.2f\n', p, s/iters);
    s2 = [s2 s/iters];
end;

save('./data/exact_coherent_threshold.m', 'ps', 'a1', 'a2', 's1', 's2');