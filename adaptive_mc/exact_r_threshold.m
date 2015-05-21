%% Probability of success curves for adaptive matrix completion
%% algorithm for different values of rank r. 

clear all;
n = 500;
r = 10;

ps = 0.02:0.02:0.6;
iters = 100;

r10 = [];
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
    r10 = [r10 s/iters];
end;

r = 20;
r20 = [];
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
    r20 = [r20 s/iters];
end;

r = 30;
r30 = [];
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
    r30 = [r30 s/iters];
end;

r = 40;
r40 = [];
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
    r40 = [r40 s/iters];
end;

save('./data/exact_r_threshold.mat', 'ps', 'r10', 'r20', 'r30', 'r40');