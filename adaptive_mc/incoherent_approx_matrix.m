function [M] = incoherent_approx_matrix(d,n,r,mu)
%% Generate a matrix with incoherent columns (not column space) and
%% with column norms distributed non-uniformly (According to
%% log-normal distribution).

l = floor(d/(mu*r));
U = zeros(d,r);
for i=1:r,
    v = zeros(d,1);
    v((i-1)*l+1:(i)*l) = 1;
    U(:,i) = v;
end;

M = zeros(d,n);

%% distribute lengths log normally;
lengths = normrnd(0,1,[n,1]);
lengths = exp(lengths);

%% distribute lengths uniformly on 0.1, 10.0;
% lengths = unifrnd(0.1, 10.0, [n,1]);

for i=1:n,
    sgn = 2*binornd(1, 0.5, [r, 1])-1;
    M(:,i) = lengths(i)*U*sgn;
end;
M = M + normrnd(0, 1.0/(d*n), [d,n]);