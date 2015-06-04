clear;
clc;
sz = [20,10,15];
n = prod(sz);
Omega = randsample(1:n, ceil(0.8 * n))';
M = rand(sz);
b = M(Omega);
sub_idx = cell(1,3);
for j = 1:3
    sub_idx{j} = 1:sz(j);
end
lambda = 5e2;
rho = 1e-2;
maxiter = 500;
tol= 1e-6 ; 

[ Z,X,obj] = admm_tensor_solver( sz,Omega, b, sub_idx, lambda, rho, maxiter, tol );
% [ X,Z,obj] = admm_low_n_rank(sz,Omega,b,lambda, rho, maxiter);
plot(obj);