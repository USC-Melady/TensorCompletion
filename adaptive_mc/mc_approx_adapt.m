function [Xhat Obs probs] = mc_approx_adapt(X, p, r)
%% Adaptive Matrix Approximation algorithm from JMLR submission. 
    d = size(X,1);
    n = size(X,2);

    Xhat = zeros(d,n);
    Obs = zeros(d,n);

    frac = 1/10;
    norms = [];
    for i=1:n,
        inds = unique(randi(d, [max(int32(d*p*frac), 1), 1]));
        Obs(inds,i) = 1;
        xto = X(inds, i);
        norms = [norms (p*frac)^(-1)*norm(xto, 'fro')^2];
    end;

    probs = norms/(sum(norms));
    for i=1:n,
        p2i = n*(1-frac)*p*probs(i);
        if p2i < 1,
            inds = unique(randi(d, [int32(d*p2i), 1]));
            Obs(inds, i) = 1;
            Xhat(inds, i) = 1/p2i*X(inds,i);
        else
            inds = (1:1:d)';
            Obs(inds, i) = 1;
            Xhat(inds, i) = X(inds, i);
        end;
    end;
    [u,s,v] = svds(Xhat, r);
    Xhat = u*s*v';