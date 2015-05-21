%% Probability of success curves for SVT algorithm for different
%% values of n. 

clear all;

r = 5;
iters = 50;

ps = 0.05:0.05:0.6;

n = 100;
m100 = [];
for p=ps,
    s = 0;
    for i=1:iters,
        X = incoherent_matrix(n,n,r,1);
        [xhat, obs] = mc_svt(X, p);
        if 1/(n^2) * norm(xhat - X, 'fro')^2 < 10^(-5),
            s = s+1;
        end;
    end;
    fprintf('p=%0.2f s=%0.2f\n', p, s/iters);
    m100 = [m100 s/iters];
end;

n = 150;
m150 = [];
for p=ps,
    s = 0;
    for i=1:iters,
        X = incoherent_matrix(n,n,r,1);
        [xhat, obs] = mc_svt(X, p);
        if 1/(n^2)*norm(xhat - X, 'fro')^2 < 10^(-5),
            s = s+1;
        end;
    end;
    fprintf('p=%0.2f s=%0.2f\n', p, s/iters);
    m150 = [m150 s/iters];
end;


n=200;
m200 = [];
for p=ps,
    s = 0;
    for i=1:iters,
        X = incoherent_matrix(n,n,r,1);
        [xhat, obs] = mc_svt(X, p);
        if 1/(n^2) * norm(xhat - X, 'fro')^2 < 10^(-5),
            s = s+1;
        end;
    end;
    fprintf('p=%0.2f s=%0.2f\n', p, s/iters);
    m200 = [m200 s/iters];
end;

save('./data/exact_svt_threshold.mat', 'ps', 'm100', 'm150', 'm200');