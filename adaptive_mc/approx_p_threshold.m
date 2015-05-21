%% Generate error curves for 4 matrix sizes as a function of p (the
%% fraction of entries sampled per column).
%% This is for the low-rank-approximation algorithm.

clear all;

r = 10;
ps = 0.01:0.01:0.4;
iters = 10;

n = 250;
m250 = [];
for p=ps,
    errs = [];
    for i=1:iters;
        %% X = incoherent_approx_matrix_random(n,n,r);
        X = incoherent_approx_matrix(n,n,r,1);
        [xhat, obs] = mc_approx_adapt(X,p,r);
        [u,s,v] = svds(X, r);
        Xstar = u*s*v';
        errs = [errs (norm(X-xhat, 'fro') - norm(X-Xstar, 'fro'))/norm(X, ...
                                                          'fro')];
    end;
    fprintf('p=%0.2f err=%0.3f\n', p, mean(errs));
    m250 = [m250, mean(errs)];
end;

n = 500;
m500 = [];
for p=ps,
    errs = [];
    for i=1:iters;
        %% X = incoherent_approx_matrix_random(n,n,r);
        X = incoherent_approx_matrix(n,n,r,1);
        [xhat, obs] = mc_approx_adapt(X,p,r);
        [u,s,v] = svds(X, r);
        Xstar = u*s*v';
        errs = [errs (norm(X-xhat, 'fro') - norm(X-Xstar, 'fro'))/norm(X, ...
                                                          'fro')];
    end;
    fprintf('p=%0.2f err=%0.3f\n', p, mean(errs));
    m500 = [m500, mean(errs)];
end;

n = 750;
m750 = [];
for p=ps,
    errs = [];
    for i=1:iters;
        %% X = incoherent_approx_matrix_random(n,n,r);
        X = incoherent_approx_matrix(n,n,r,1);
        [xhat, obs] = mc_approx_adapt(X,p,r);
        [u,s,v] = svds(X, r);
        Xstar = u*s*v';
        errs = [errs (norm(X-xhat, 'fro') - norm(X-Xstar, 'fro'))/norm(X, ...
                                                          'fro')];
    end;
    fprintf('p=%0.2f err=%0.3f\n', p, mean(errs));
    m750 = [m750, mean(errs)];
end;

n = 1000;
m1000 = [];
for p=ps,
    errs = [];
    for i=1:iters;
        X = incoherent_approx_matrix(n,n,r,1);
        [xhat, obs] = mc_approx_adapt(X,p,r);
        [u,s,v] = svds(X, r);
        Xstar = u*s*v';
        errs = [errs (norm(X-xhat, 'fro') - norm(X-Xstar, 'fro'))/norm(X, ...
                                                          'fro')];
    end;
    fprintf('p=%0.2f err=%0.3f\n', p, mean(errs));
    m1000 = [m1000, mean(errs)];
end;

save('./data/approx_p_threshold.mat', 'ps', 'm250', 'm500', 'm750', ...
     'm1000');

