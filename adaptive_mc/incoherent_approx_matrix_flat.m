function [M] = incoherent_approx_matrix_flat(d,n,r,mu)
%% Generate a matrix with incoherent columns (not column space) and
%% with column norms distributed uniformly from [0.9, 1.1].
l = floor(d/(mu*r));
U = zeros(d,r);
for i=1:r,
    v = zeros(d,1);
    v((i-1)*l+1:(i)*l) = 1;
    U(:,i) = v;
end;

M = zeros(d,n);

%% distribute lengths uniformly between 0.9 and 1.1
lengths = abs(unifrnd(0.9,1.1,[n,1]));

for i=1:n,
    sgn = 2*binornd(1, 0.5, [r, 1])-1;
    M(:,i) = lengths(i)*U*sgn;
end;
M = M + normrnd(0, 1.0/(d*n), [d,n]);