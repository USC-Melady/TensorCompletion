function [M] = incoherent_approx_matrix_random(d,n,r)
%% Generate a rank r+1 matrix from random gaussian entries.
%% The idea here is that the energy in the lower ranks is
%% concentrated. 
M = normrnd(0, 1, [d,n]);
[U,S,V] = svds(M, d);
s = diag(repmat(1, 1, d));
M = U*s*V';