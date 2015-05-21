function [M] = incoherent_matrix(d,n,r,mu)
%% Generate a low rank matrix with incoherent row and column
%% space. 
l = floor(d/(mu*r));
U = zeros(d,r);
for i=1:r,
    v = zeros(d,1);
    v((i-1)*l+1:(i)*l) = 1;
    U(:,i) = v;
end;
V = normrnd(0, 1, n, r);
M = U*V';
