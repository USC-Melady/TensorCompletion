function [ obj_val ] = eval_objective( X,Y,Z,M,Omega, submat_idx, lambda, rho )
%EVAL_OBJECTIVE Summary of this function goes here
%   Detailed explanation goes here

m = length(X);

t1= 0;
for i = 1:m
    t1 = t1 +  trace(sqrt(X{i}'*X{i}));
end
t2 = lambda * norm(( Z-M ) .* Omega, 'fro')^2;
t3 =0;

for i = 1:m
    row_idx = submat_idx{i,1};
    col_idx = submat_idx{i,2};
    t3 = t3 + trace(Y{i}' * (X{i}- Z(row_idx, col_idx)));
    t3 = t3 + 0.5* norm(X{i}-Z(row_idx, col_idx),'fro');    
end
t3 = rho*t3;

obj_val = t1+ t2 + t3;

end

