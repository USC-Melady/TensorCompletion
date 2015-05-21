%% Generate error curves for 4 target ranks as a function of p (the
%% fraction of entries sampled per column).
%% This is for the low-rank-approximation algorithm.

clear all;

n = 500;
iters = 100;
ps = 0.02:0.02:0.6;

r = 10;
m10 = [];
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
    m10 = [m10, mean(errs)];
end;


r = 20;
m20 = [];
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
    m20 = [m20, mean(errs)];
end;

r = 30;
m30 = [];
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
    m30 = [m30, mean(errs)];
end;

save('./data/approx_r_threshold.mat', 'ps', 'm10', 'm20', 'm30');