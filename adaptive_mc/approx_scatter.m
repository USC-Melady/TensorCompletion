%% Generate data for a scatter plot of error vs number of
%% observations for the low-rank approximation algorithm. 

clear all;

n = 200;
ps = 0.01:0.01:0.1;
iters = 10;
r = 10;

eadapt = [];
madapt = [];
esoft = [];
msoft = [];
for p=ps,
    fprintf('p=%0.2f\n', p);
    for i=1:iters,
        X = incoherent_approx_matrix(n,n,r,1);
        [u,s,v] = svds(X,r);
        Xstar = u*s*v';
        [xhat, obs] = mc_approx_adapt(X,p,r);
        eadapt = [eadapt abs(norm(X - xhat, 'fro') - norm(X - Xstar, ...
                                                          'fro'))/norm(X, ...
                                                          'fro')];
        madapt = [madapt nnz(obs)/(n*n)];
        [xhat, obs] = mc_soft_threshold(X,p,r);
        esoft = [esoft abs(norm(X - xhat, 'fro') - norm(X - Xstar, ...
                                                          'fro'))/norm(X, ...
                                                          'fro')];
        msoft = [msoft nnz(obs)/(n*n)];
    end;
end;