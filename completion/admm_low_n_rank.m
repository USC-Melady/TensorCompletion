function [ X,Z,obj  ] = admm_low_n_rank( sz,Omega, b, lambda, rho, maxiter )
%ADMM_LOW_N_RANK Summary of this function goes here
%   Detailed explanation goes here
N= 3;
X = tenzeros(sz);
Z = cell(1,N);
Y = cell(1,N);
comp_ind = [1:prod(sz)]';
ind_2  = setdiff(comp_ind, Omega);

X_old = X;
obj = zeros(maxiter,1);
for n = 1:N
    Z{n} = tenzeros(sz);
    Y{n} = tenzeros(sz);
end
for iter = 1: maxiter
    fprintf('\b\b\b\b%4d',iter); 
    tmp = tenzeros(sz);
    for n = 1:N
        tmp = tmp + Y{n} + rho * Z{n};
    end
    
    X(Omega) = ( tmp(Omega) + lambda * b ) ./ ( lambda + N *rho) ;
    X(ind_2) =  tmp(ind_2)./ (N*rho);
    
    for n = 1:N
        X_n = tenmat(X,n);
        X_n = X_n.data;
        Y_n = tenmat(Y{n},n);
        Y_n = Y_n.data;
        Z_n = shrink(X_n - 1/rho * Y_n ,1/rho );
        Z{n} = refold(Z_n, n, sz);
    end
    
    
    for n = 1:N
        Y{n} = Y{n} + rho * (Z{n} -X );
    end
    
    obj(iter) = norm(X - X_old) / norm(X_old);
    X_old = X;
end

end

