clear; 
clc;
global VERBOSE;
VERBOSE =1;
%% gen data
n1 = 300; n2 = 500;
sz = [n1, n2];
M = rand(sz);
% low rank projection
r = 20;
n = prod(sz); % total number of entries
p = 0.65; %log(n)* r/n;
[U, S, V] = svd(M);
S(r+1:end,r+1:end)=0;
M = U*S*V';

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

%% construct Omega_i
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

%% admm
lambda = 1e1;
rho = 1e-2;
max_iter = 500;

Omega_ind = find(Omega>0);
data = M(Omega_ind);
[Z, X, obj] = admm_solver(sz,Omega_ind,data, submat_idx, lambda, rho,max_iter );
