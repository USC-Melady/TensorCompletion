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
obj = zeros(max_iter,1);
m = size(submat_idx,1); %number of components

% initialize
Y = cell(m,1); % dual variable 
X = cell(m,1); % auxiliary variables for each part
Z = M; % auxiliary varible for complete
for i = 1:m
     row_idx = submat_idx{i,1};
     col_idx = submat_idx{i,2};
     X{i} = rand(length(row_idx), length(col_idx));
     Y{i} = zeros(length(row_idx), length(col_idx));
end

if VERBOSE==1, fprintf('\nIteration:   '); end

for iter = 1: max_iter
    if VERBOSE==1, fprintf('\b\b\b\b%4d',iter);  end

% solve for each X_i separately
    for i = 1:m
        row_idx = submat_idx{i,1};
        col_idx = submat_idx{i,2};
        tmp = Z(row_idx, col_idx);
        X{i} = shrink(tmp - Y{i}, 1/rho); 
    end
% solve Z: closed form solution
    for i = 1:m
        row_idx = submat_idx{i,1};
        col_idx = submat_idx{i,2};
        tmp = 2* lambda* M.* Omega;
        tmp = tmp (row_idx, col_idx);
        Z(row_idx,col_idx) = tmp +...
                            rho * (X{i} + Y{i} );
    end
     
% update Y
    for i = 1:m
        row_idx = submat_idx{i,1};
        col_idx = submat_idx{i,2};
        Y{i} = Y{i}  + (X{i} - Z(row_idx, col_idx) );
    end
    % objective funciton 
    obj(iter) = eval_objective(X,Y,Z,M,Omega, submat_idx, lambda, rho);

end
if VERBOSE==1, fprintf('\n'); end

end

