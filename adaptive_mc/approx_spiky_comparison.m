%% Compare the adaptive algorithm with soft-thresholding on
%% matrices with non-uniform column norm. 

clear all;

ps = 0.01:0.01:0.5;
iters = 100;
n = 500;
r = 10;

adapt = [];
soft = [];
for p=ps,
    sa = [];
    st = [];
    for i=1:iters,
        X = incoherent_approx_matrix(n,n,r,1);
        [xadapt, obs] = mc_approx_adapt(X, p, r);
        [xsoft, obs] = mc_soft_threshold(X, p, r);
        [u,s,v] = svds(X,r);
        Xstar = u*s*v';
        sa = [sa abs(norm(X - xadapt, 'fro') - norm(X - Xstar, ...
                                                 'fro'))/norm(X, 'fro')];
        st = [st abs(norm(X - xsoft, 'fro') - norm(X - Xstar, ...
                                                'fro'))/norm(X, 'fro')];
    end;
    fprintf('p=%0.2f sa=%0.3f st=%0.3f\n', p, mean(sa), mean(st));
    adapt = [adapt mean(sa)];
    soft = [soft mean(st)];
end;

save('./data/approx_spiky_comparison.mat', 'ps', 'adapt', 'soft');