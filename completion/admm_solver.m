function [ X, Z, obj] = admm_solver( M,Omega, submat_idx, lambda, rho,max_iter )
%ADMM_SOLVER Summary of this function goes here
%   Detailed explanation goes here
% M is the full matrix
% Omega is the linear indexing 
global VERBOSE
if isempty(VERBOSE)
    % -- feel free to change these 'verbosity' parameters
    % VERBOSE = false;
    VERBOSE = 1;    % a little bit of output
    % VERBOSE = 2;    % even more output
end

sz = size(M);
obj = zeros(max_iter,1);
m = size(submat_idx,1); %number of components

% initialize
Y = cell(m,1); % dual variable 
X = cell(m,1); % auxiliary variables for each part
Z = zeros(sz); % auxiliary varible for complete
P = zeros(sz); % Submatrix Projection

for i = 1:m
     row_idx = submat_idx{i,1};
     col_idx = submat_idx{i,2};
     X{i} = zeros(length(row_idx), length(col_idx));
     Y{i} = zeros(length(row_idx), length(col_idx));
     P(row_idx, col_idx) = P(row_idx,col_idx) +1;
end

if VERBOSE==1, fprintf('\nIteration:   '); end

for iter = 1: max_iter
    if VERBOSE==1, fprintf('\b\b\b\b%4d',iter);  end

% solve for each X_i separately
    for i = 1:m
        row_idx = submat_idx{i,1};
        col_idx = submat_idx{i,2};
        tmp = Z(row_idx, col_idx);
        X{i} = shrink(tmp - 1/rho* Y{i}, 1/rho); 
    end
% solve Z: closed form solution for each i
    tmp = zeros(sz);
    for i = 1:m
        row_idx = submat_idx{i,1};
        col_idx = submat_idx{i,2};
        tmp(row_idx,col_idx) = tmp(row_idx, col_idx) + rho * X{i} + Y{i};
        
    end
    tmp = tmp +  lambda * M.*Omega;
    % branching
    Z = zeros(sz);% both zero
    obs_ind = find(Omega==1);
    sub_ind = find(P>=1);
    
    intersect_ind = intersect(obs_ind, sub_ind);
    Z(intersect_ind) = tmp(intersect_ind) / (lambda + rho);
    
    ind_1 = setdiff(obs_ind, intersect_ind);
    Z(ind_1) = tmp(ind_1)/(lambda);

    ind_2 = setdiff(sub_ind, intersect_ind);
    Z(ind_2) = tmp(ind_2)/rho;
     
% update Y
    for i = 1:m
        row_idx = submat_idx{i,1};
        col_idx = submat_idx{i,2};
        Y{i} = Y{i}  + rho * (X{i} - Z(row_idx, col_idx) );
    end
    % objective funciton 
    obj(iter) = eval_objective(X,Y,Z,M,Omega, submat_idx, lambda, rho);

end
if VERBOSE==1, fprintf('\n'); end

end

