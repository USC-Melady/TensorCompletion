%% Probability of success curves for adaptive matrix completion
%% algorithm for different values of coherence \mu. 

clear all;
n = 500;
r = 10;
iters = 100;
ps = 0.02:0.02:0.6;

mu = 1;
mu1 = [];
for p=ps,
    s = 0;
    for i=1:iters,
        X = incoherent_matrix(n,n,r,mu);
        [xhat, obs] = mc_exact_adapt(X, p);
        if norm(xhat - X, 'fro') < 0.001,
            s = s+1;
        end;
    end;
    fprintf('p=%0.2f s=%0.2f\n', p, s/iters);
    mu1 = [mu1 s/iters];
end;

mu = 2;
mu2 = [];
for p=ps,
    s = 0;
    for i=1:iters,
        X = incoherent_matrix(n,n,r,mu);
        [xhat, obs] = mc_exact_adapt(X, p);
        if norm(xhat - X, 'fro') < 0.001,
            s = s+1;
        end;
    end;
    fprintf('p=%0.2f s=%0.2f\n', p, s/iters);
    mu2 = [mu2 s/iters];
end;

mu = 3;
mu3 = [];
for p=ps,
    s = 0;
    for i=1:iters,
        X = incoherent_matrix(n,n,r,mu);
        [xhat, obs] = mc_exact_adapt(X, p);
        if norm(xhat - X, 'fro') < 0.001,
            s = s+1;
        end;
    end;
    fprintf('p=%0.2f s=%0.2f\n', p, s/iters);
    mu3 = [mu3 s/iters];
end;

mu = 4;
mu4 = [];
for p=ps,
    s = 0;
    for i=1:iters,
        X = incoherent_matrix(n,n,r,mu);
        [xhat, obs] = mc_exact_adapt(X, p);
        if norm(xhat - X, 'fro') < 0.001,
            s = s+1;
        end;
    end;
    fprintf('p=%0.2f s=%0.2f\n', p, s/iters);
    mu4 = [mu4 s/iters];
end;

save('./data/exact_mu_threshold.mat', 'ps', 'mu1', 'mu2', 'mu3', 'mu4');