function M = shrink(Y, tau)
%      [U L V]=mySVD(A);
%     options.tol = 1e-5;  % REducing the accuracy for faster convergence
%     [U, L, V] = svds(Y, min(size(Y)), 'L', options);
%     eig=diag(L)- tau;
%     eig(eig<0)=0;
%     eigM=diag(eig);
%     L(1:size(eigM,1), 1:size(eigM,2))=eigM;
%     M=U*L*V';

% Candeas routine on soft-thresholdin operator
    [U,Sigma,V] = svd(Y ,'econ'); 
    sigma = diag(Sigma); r = sum(sigma > tau);
    U = U(:,1:r); V = V(:,1:r); sigma = sigma(1:r) - tau; 
    M=U*diag(sigma)*V';

end

