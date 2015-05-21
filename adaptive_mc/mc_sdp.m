function [Xhat Obs] = mc_sdp(X, p, r)
%% SDP algorithm of Negahban and Wainwright. Requires cvx. 
%% lambda should be 4 L sigma \sqrt{n/m \log n}
n = size(X, 2);
d = size(X, 1);
inds = unique(randi(n*d, [int32(n*d*p), 1]));
obs = zeros(d*n,1);
obs(inds) = 1;
Obs = reshape(obs, [d,n]);

Xtilde = X.*Obs;

n = nnz(Obs);
inds = find(Obs);
[n1,n2] = size(X);
lambda = 4 * sqrt(0.5*(n1+n2)/n * log(0.5*(n1+n2)));
cvx_begin sdp
  variable X(n1,n2)
  minimize ( 1.0/(2*n) * square_pos(norm(X(inds) - Xtilde(inds),2)) + lambda * norm_nuc(X) )
  subject to
    norm(vec(X), Inf) <= 1;
cvx_end

[u,s,v] = svds(X, r);
Xhat = u*s*v';