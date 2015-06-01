function [ Z,X, obj ] = admm_solver_tensor( sz,Omega, b, submat_idx, lambda, rho, maxiter, tol )
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
[m,N] = size(submat_idx); %number of components, dimension of a tensor

% initialize
Y = cell(m,N); % dual variable 
X = cell(m,N); % auxiliary variables for each part
Z = tenzeros(sz); % auxiliary varible for complete
P = tenzeros(sz); % Submatrix Projection


Z_old = Z;
for i = 1:m
    subs = [];
    subs_len = zeros(1,N);
    for n = 1:N
     subs = [subs,submat_idx{i,n}];
     subs_len(n) =length(submat_idx{i,n});
    end
    
    for n = 1:N
     X{i,n} = tenzeros(subs_len);
     Y{i,n} = tenzeros(subs_len);
    end
    
     P(subs) = P(subs) +1;%extract subtensor
end
completable_ind = find(P>=1);
intersect_ind = intersect(Omega, completable_ind);
ind_1 = setdiff(Omega, intersect_ind); % observed, not completable
ind_2 = setdiff(completable_ind, intersect_ind); % unobserved, completable


if VERBOSE==1, fprintf('\nIteration:   '); end

for iter = 1: maxiter
    if VERBOSE==1, fprintf('\b\b\b\b%4d',iter);  end

% solve for each X_i separately
    for i = 1:m
        for n = 1:N
        row_idx = submat_idx{i,n};
        end
        X{i} = Z(row_idx, col_idx);
        X{i} = shrink(X{i} - 1/rho* Y{i}, 1/rho); 
    end
% solve Z: closed form solution for each i
    tmp = zeros(sz);
    for i = 1:m
        row_idx = submat_idx{i,1};
        col_idx = submat_idx{i,2};
        tmp(row_idx,col_idx) = tmp(row_idx, col_idx) + rho * X{i} + Y{i};      
    end
    tmp(Omega) = tmp(Omega) +  lambda * b;
    % branching
    Z = zeros(sz);% unobserved and imcompletable: zero
    Z(intersect_ind)  = tmp(intersect_ind)./(lambda + rho * P(intersect_ind));% observed and completable   
    Z(ind_1) = tmp(ind_1)/(lambda);
    Z(ind_2) = tmp(ind_2)/(rho);
     
% update Y
    for i = 1:m
        row_idx = submat_idx{i,1};
        col_idx = submat_idx{i,2};
        Y{i} = Y{i}  + rho * (X{i} - Z(row_idx, col_idx) );
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

