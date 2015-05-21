function [Xhat, Obs] = mc_soft_threshold(X, p, r)
%% Soft thresholding algorithm used by Koltchinskii, Lounici, and Tsybakov.
  n = size(X,2);
  d = size(X, 1);

  inds = unique(randi(n*d, [int32(n*d*p), 1]));
  obs = zeros(d*n,1);
  obs(inds) = 1;
  Obs = reshape(obs, [d,n]);

  Xtilde = 1/p*X.*Obs;

  [u,s,v] = svds(Xtilde,r);
  Xhat = u*s*v';
