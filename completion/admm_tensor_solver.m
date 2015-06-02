function [ Z,X,obj] = admm_tensor_solver( sz,Omega, b, sub_idx, lambda, rho, maxiter, tol )
%ADMM_SOLVER_TENSOR : tensor completion algorithm using ADM
%generalizeion using tensor package

global VERBOSE
if isempty(VERBOSE)
    % -- feel free to change these 'verbosity' parameters
    % VERBOSE = false;
    VERBOSE = 1;    % a little bit of output
    % VERBOSE = 2;    % even more output
end

if nargin < 8 || isempty(tol)
    tol = 1e-4;
end
if nargin < 7 || isempty(maxiter)
    maxiter = 500;
end

obj = zeros(maxiter+1,1);
[m,N] = size(sub_idx); %number of components, dimension of a tensor

% initialize
Y = cell(m,N); % dual variable , N is three for now
X = cell(m,N); % auxiliary variables for each part
Z = tenzeros(sz); % auxiliary varible for complete
P = tenzeros(sz); % Submatrix Projection


Z_old = Z;

for i = 1:m
    sub_sz = zeros(1,N);
    for n = 1:N
       sub_sz(n) =length(sub_idx{i,n});
    end
    
    for n = 1:N
      X{i,n} = tenzeros(sub_sz);
      Y{i,n} = tenzeros(sub_sz);
    end
    
    P(sub_idx{1}, sub_idx{2}, sub_idx{3}) = P(sub_idx{1}, sub_idx{2}, sub_idx{3}) +1;% extract subtensor
end
completable_ind = find(P>=1);
intersect_ind = intersect(Omega, completable_ind);
ind_1 = setdiff(Omega, intersect_ind); % observed, not completable
ind_2 = setdiff(completable_ind, intersect_ind); % unobserved, completable


if VERBOSE==1, fprintf('\nIteration:   '); end

for iter = 1: maxiter
    if VERBOSE==1, fprintf('\b\b\b\b%4d',iter);  end

   % solve for each X_i separallel
   for i = 1:m
        for n = 1:N
            Z_sub = Z(sub_idx{1}, sub_idx{2}, sub_idx{3});
            % mode n unfolding of Z
            Z_sub_n = tenmat(Z_sub,n);
            Y_n = tenmat(Y{i,n},n);
            X_sub_n = shrink(Z_sub_n.data - 1/rho* Y_n.data, 1/rho); 
            X{i,n} = tensor (X_sub_n, sz);
        end
    end
% solve Z: closed form solution for each i
    tmp = tenzeros(sz);
    for i = 1:m
        for n = 1:N
           tmp(sub_idx{1}, sub_idx{2}, sub_idx{3}) = tmp(sub_idx{1}, sub_idx{2}, sub_idx{3}) + rho * X{i,n} + Y{i,n};      
        end
    end
    tmp(Omega) = tmp(Omega) +  lambda * b;
    % branching
    Z = tenzeros(sz);% unobserved and imcompletable: zero
    Z(intersect_ind)  = tmp(intersect_ind)./(lambda + N* rho * P(intersect_ind));% observed and completable   
    Z(ind_1) = tmp(ind_1)/(lambda);
    Z(ind_2) = tmp(ind_2)/(N*rho);
     
% update Y
    for i = 1:m
        for n = 1:n
        Y{i,n} = Y{i,n}  + rho * (X{i,n} - Z(sub_idx{1}, sub_idx{2}, sub_idx{3}) );
        end
    end
    % objective funciton 
    % obj(iter+1) = eval_objective(X,Z,Omega,b, lambda);
    
    % convergence criteria
    stop_criterion = norm ( Z(:) - Z_old(:))  /  norm(Z_old);
    if (stop_criterion< tol)
        break
    end
    obj(iter)= stop_criterion;
    Z_old = Z;
end

if VERBOSE==1, fprintf('\n'); end
end

