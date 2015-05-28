function [ obj_val ] = eval_objective( X,Z,Omega,data, lambda)
%EVAL_OBJECTIVE : evaluate objective function of 
% lambda/2 \| P_\Omega(Z) - P_\Omega(M)\|_F^2 + \sum_{i=1}^m \|X_i\|_*
% subject to X_i = P_{\Omega_i}(Z)

m = length(X);

t1= 0;
for i = 1:m
    t1 = t1 +  trace_norm (X{i});
end
t2 = 0.5* lambda * norm( Z(Omega)- data , 'fro')^2;


% t3 =0;
% for i = 1:m
%     row_idx = submat_idx{i,1};
%     col_idx = submat_idx{i,2}; 
%     t3 = t3 + trace(Y{i}' * (X{i}- Z(row_idx, col_idx)));
%     t3 = t3 + 0.5*rho* norm(X{i}-Z(row_idx, col_idx),'fro')^2;    
% end

obj_val = t1+ t2;


function val = trace_norm(X)
    [U , S , V] = svd(X);
    val = sum(diag(S));

end
end

