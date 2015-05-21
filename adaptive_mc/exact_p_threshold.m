%% Probability of success curves for adaptive matrix completion
%% algorithm for different values of matrix size n.

clear all;
r = 10;
ps = 0.01:0.01:0.4;
iters = 200;


n = 250;
m250 = [];
for p=ps,
    s = 0;
    for i=1:iters,
        X = incoherent_matrix(n,n,r,1);
        [xhat, obs] = mc_exact_adapt(X, p);
        if norm(xhat - X, 'fro') < 0.001,
            s = s+1;
        end;
    end;
    fprintf('p=%0.2f s=%0.2f\n', p, s/iters);
    m250 = [m250 s/iters];
end;

n = 500;
m500 = [];
for p=ps,
    s = 0;
    for i=1:iters,
        X = incoherent_matrix(n,n,r,1);
        [xhat, obs] = mc_exact_adapt(X, p);
        if norm(xhat - X, 'fro') < 0.001,
            s = s+1;
        end;
    end;
    fprintf('p=%0.2f s=%0.2f\n', p, s/iters);
    m500 = [m500 s/iters];
end;

n = 750;
m750 = [];
for p=ps,
    s = 0;
    for i=1:iters,
        X = incoherent_matrix(n,n,r,1);
        [xhat, obs] = mc_exact_adapt(X, p);
        if norm(xhat - X, 'fro') < 0.001,
            s = s+1;
        end;
    end;
    fprintf('p=%0.2f s=%0.2f\n', p, s/iters);
    m750 = [m750 s/iters];
end;

n = 1000;
m1000 = [];
for p=ps,
    s = 0;
    for i=1:iters,
        X = incoherent_matrix(n,n,r,1);
        [xhat, obs] = mc_exact_adapt(X, p);
        if norm(xhat - X, 'fro') < 0.001,
            s = s+1;
        end;
    end;
    fprintf('p=%0.2f s=%0.2f\n', p, s/iters);
    m1000 = [m1000 s/iters];
end;

save('data/exact_p_threshold.mat', 'ps', 'm250', 'm500', 'm750', 'm1000');