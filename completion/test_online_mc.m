clear; 
clc;
%% gen data

n1 = 100; n2 = 200;
sz = [n1, n2];
X = rand(sz);
% low rank projection
r = 20;
n = prod(sz); % total number of entries
p = 0.65; %log(n)* r/n;
[U, S, V] = svd(X);
S(r+1:end,r+1:end)=0;
X = U*S*V';

%(1) uniform random
% k = ceil (ratio * n );
% Omega_ind = randsample(n, k );  %observed entries
% Omega = zeros(n1,n2);
% Omega(Omega_ind) = 1;
% [I,J] = ind2sub(sz,Omega_ind);

%(2) block overlapp
Omega = zeros(n1, n2);

n1_tmp = 20;
Omega_ind_tmp = randsample(n1_tmp *n2, ceil(p *n1_tmp *n2 ));
Omega_tmp = zeros(n1_tmp, n2);
Omega_tmp(Omega_ind_tmp) = 1;
Omega(1:n1_tmp, :) = Omega_tmp;

n2_tmp = 30;
Omega_ind_tmp = randsample(n1*n2_tmp, ceil (p * n1*n2_tmp) );
Omega_tmp = zeros(n1, n2_tmp);
Omega_tmp(Omega_ind_tmp) = 1;
Omega(:, 1:n2_tmp) =  Omega(:, 1:n2_tmp) + Omega_tmp;
[I, J] = find(Omega>0);
Omega(Omega>1)= 1;

% simulate observed entries 
% % non-overlapping
% A = eye(10);
% B = ones(5);
% tmp = kron(A,B) ;
% Omega = zeros(sz);
% Omega(1:50,51:end) = tmp;
% Omega(51:end, 1:50)= tmp;

% overlapping
% A = eye(10);
% B = ones(8);
% tmp = kron(A,B);
% Omega = zeros(sz);
% Omega(1:80,21:end) = tmp;
% Omega(21:end,1:80) = Omega(21:end,1:80)+  tmp;
% Omega(Omega>1) =1;
% 

%% trim degree

tic;
% construct adjacent matrix
G = sparse([I,J+n1],[J+n1,I], 1, n1+n2, n1+n2);
G = full(G);


v_s = []; % set of removed vertices
degree = sum(G);
v = find(degree < r );
v_s=[v_s,v];
while( ~isempty(v))
    % remove vertex
    G(v,:) = 0;
    G(:,v) = 0; 
    % check degree
    degree = sum(G);
    v = setdiff(v_s, find(degree < r ));    
    v_s = [v_s, v];
%     disp(length(v_s));
end

if(length(v_s) ==size(G,1))
    fprintf('No completable submatrix\n');
end

%% find strongly-connected component

graph = sparse(G);
[K,C] = graphconncomp(graph); % BUG: too many singletons

%% optimization

df = r*(n1+n2-r);
oversampling = 5; 
m = min(5*df,round(.99*n1*n2) ); 
p  = m/(n1*n2);

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
maxiter = 100; 
tol = 1e-4;

X_c = zeros(sz);
submat_idx = [];
for i = 1:K
    blk_idx = find(C==i);
    %ignore singleton
    n_i = length(blk_idx);
    if (n_i ==1)
        continue;
    end  
    row_idx = {blk_idx(blk_idx<=n1)};
    col_idx = {blk_idx(blk_idx>n1)-n1};
    submat_idx = [submat_idx;row_idx, col_idx];
end

m = size(submat_idx,1);
fprintf('%d connnected components found\n', m);

% optmization
for i = 1: m
    row_idx =  submat_idx{i,1};
    col_idx =  submat_idx{i,2};
    n1_i = length(row_idx);
    n2_i = length(col_idx);
    fprintf('component %d X %d \n', n1_i,n2_i);
    Omega_i=  Omega( row_idx, col_idx);
    Omega_i_ind = find(Omega_i==1); % must be linear index
    Xi = X(row_idx,col_idx);
    data_i = Xi(Omega_i_ind);
    % solve svt for each block
    [U,S,V,~] = SVT([n1_i, n2_i],Omega_i_ind,data_i,tau,delta,maxiter,tol);
    Xi_c = U*S*V';
    X_c(row_idx, col_idx) = Xi_c;   
end
runtime_1= toc;

%% compare with full matrix completion
Omega_ind = find(Omega==1);
data = X(Omega_ind);
[U,S,V,~] = SVT(sz,Omega_ind,data,tau,delta,maxiter,tol);
X_c2 = U*S*V';
runtime_2 = toc;


[U,S,V,~] = SVT(sz,Omega_ind,data,tau,delta,maxiter,tol);
X_c3 = U*S*V';
runtime_3 = toc;

%% evaluate
rmse_1 = eval_RMSE( X, X_c, submat_idx );
fprintf('RMSE submatrix: %d, run time: %d \n',rmse_1, runtime_1 );

rmse_2 = eval_RMSE( X, X_c2, submat_idx );
fprintf('RMSE full: %d, run_time: %d \n', rmse_2, runtime_2);

rmse_3 = eval_RMSE( X, X_c3, submat_idx );
fprintf('RMSE full: %d, run_time: %d \n', rmse_3, runtime_3);

% TD: overlapp components
% TD: compare with adaptive mc
% TD: generalize to tensor


