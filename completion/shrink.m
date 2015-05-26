function M = shrink(A, s)
% [U L V]=mySVD(A);
options.tol = 1e-5;  % REducing the accuracy for faster convergence
[U, L, V] = svds(A, min(size(A)), 'L', options);
eig=diag(L)-s;
eig(eig<0)=0;
eigM=diag(eig);
L(1:size(eigM,1), 1:size(eigM,2))=eigM;
M=U*L*V';
end

