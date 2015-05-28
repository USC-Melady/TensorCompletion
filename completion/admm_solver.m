function [Z,X, obj] = admm_solver(sz,Omega, b, submat_idx, lambda, rho, maxiter, tol )
%ADMM_SOLVER: ADMM solver for graph-based matrix completion
% Input
% sz: size of the orignal matrix [n1,n2]
% Omega: linear indice of observed entries
% data: observed entries corresponding to position in Omega
% submat_idx: m x 2 cell array of m completable components, 
              % first column row indice, second column column indice
% lambda: positive number, Lagrange multiplier: larger lambda, faster,
            % lower error on Omega
      
% rho: positive number, augmented Lagrange multiplier: smaller rho, lower
            % recover error 
% max_iter: max number of iteration

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
m = size(submat_idx,1); %number of components

% initialize
Y = cell(m,1); % dual variable 
X = cell(m,1); % auxiliary variables for each part
Z = zeros(sz); % auxiliary varible for complete
P = zeros(sz); % Submatrix Projection


Z_old = Z;
for i = 1:m
     row_idx = submat_idx{i,1};
     col_idx = submat_idx{i,2};
     X{i} = zeros(length(row_idx), length(col_idx));
     Y{i} = zeros(length(row_idx), length(col_idx));
     P(row_idx, col_idx) = P(row_idx,col_idx) +1;
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
        row_idx = submat_idx{i,1};
        col_idx = submat_idx{i,2};
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
    Z_old = Z;
end
if VERBOSE==1, fprintf('\n'); end

end

