function [M] = coherent_matrix(d,n,r,mu)
%% Generate a low rank matrix with coherent row space. 
l = floor(d/(mu*r));
U = zeros(d,r);
for i=1:r,
    v = zeros(d,1);
    v((i-1)*l+1:(i)*l) = 1;
    U(:,i) = v;
end;
y = randsample(1:n, r);
V = zeros(n,r);
for i=1:r,
    V(y(i), i) = 1;
end;
M = U*V';
