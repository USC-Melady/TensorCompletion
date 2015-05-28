
clear;
clc;
%% Setup a matrix
randn('state',2009);
rand('state',2009);

% if you're daring, try with complex numbers:
COMPLEX = false;

n1 = 150; n2 = 300; r = 10;
if COMPLEX
    M = (randn(n1,r)+1i*randn(n1,r))*(randn(r,n2)+1i*randn(r,n2))/2;
else
    M = randn(n1,r)*randn(r,n2);
end

df = r*(n1+n2-r);
oversampling = 5; 
m = min(5*df,round(.99*n1*n2) ); 
p  = m/(n1*n2);

Omega = randsample(n1*n2,m);  % this requires the stats toolbox
% a workaround, if you don't have the stats toolbox, is this:
%   Omega = randperm(n1*n2); Omega = Omega(1:m);

data = M(Omega);
% add in noise, if desired
sigma = 0;
% sigma = .05*std(data);
data = data + sigma*randn(size(data));

fprintf('Matrix completion: %d x %d matrix, rank %d, %.1f%% observations\n',...
    n1,n2,r,100*p);
fprintf('\toversampling degrees of freedom by %.1f; noise std is %.1e\n',...
    m/df, sigma );
if ~isreal(M), disp('Matrix is complex'); end
%% Set parameters and solve

tau = 5*sqrt(n1*n2); 
delta = 1.2/p;    
%{
 if n1 and n2 are very different, then
   tau should probably be bigger than 5*sqrt(n1*n2)

 increase tau to increase accuracy; decrease it for speed

 if the algorithm doesn't work well, try changing tau and delta
   i.e. if it diverges, try a smaller delta (e.g. delta < 2 is a 
   safe choice, but the algorithm may be slower than necessary).
%}
maxiter = 500; 
tol = 1e-4;
%% Approximate minimum nuclear norm solution by SVT algorithm
% Note: SVT, as called below, is setup for noiseless data 
%   (i.e. equality constraints).

fprintf('\nSolving by SVT...\n');
tic
[U,S,V,numiter] = SVT([n1 n2],Omega,data,tau,delta,maxiter,tol);
toc
    
X = U*S*V';
    
% Show results
fprintf('The recovered rank is %d\n',length(diag(S)) );
fprintf('The relative error on Omega is: %d\n', norm(data-X(Omega))/norm(data))
fprintf('The relative recovery error is: %d\n', norm(M-X,'fro')/norm(M,'fro'))
fprintf('The relative recovery in the spectral norm is: %d\n', norm(M-X)/norm(M))
%% ADMM
submat_idx = cell(1,2);
submat_idx{1,1}= 1:n1;
submat_idx{1,2}= 1:n2;

RelError = [];
rho = 1e-2;
lambda = 1e3;
max_iter = 100;
[X, ~, obj] = admm_solver([n1,n2], Omega, data , submat_idx, lambda, rho,max_iter );

% Show results
fprintf('rho %d\n',rho);
fprintf('The relative error on Omega is: %d\n', norm(data-X(Omega))/norm(data))
fprintf('The relative recovery error is: %d\n', norm(M-X,'fro')/norm(M,'fro'))
fprintf('The relative recovery in the spectral norm is: %d\n', norm(M-X)/norm(M))
