clear; clc;
%% gen data
sz = [100, 100];
X = rand(sz);

n = prod(sz); % total number of entries
m = 3; % total number of blocks
k = ceil(n * 0.8);


% simulate observed entries
A = eye(10);
B = ones(5);
tmp = kron(A,B) ;
Omega = zeros(sz);

% non-overlapping
Omega(1:50,51:end) = tmp;
Omega(51:end, 1:50)= tmp;

%% trim degree
r = 5;
v = [];
while( length(v)~= length(Omega) )
    degree = sum(Omega);
    v = find(degree >= r );
    Omega = Omega( v, v);
end


%% find strongly-connected component

graph = sparse(Omega);
[K,C] = graphconncomp(graph);

%% optimization
n1 = sz(1);
n2 = sz(2);

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
maxiter = 500; 
tol = 1e-4;

X_c = zeros(sz);
for i = 1:K
    blk_idx = find(C==i);
    n_i = length(blk_idx);
    Omega_i=  Omega(blk_idx, blk_idx);
    Omega_i = find(Omega_i==1); % must be linear index
    Xi = X(blk_idx,blk_idx);
    data_i = Xi(Omega_i);
    % solve svt for each block
    [U,S,V,numiter] = SVT(n_i,Omega_i,data_i,tau,delta,maxiter,tol);
    Xi_c = U*S*V';
    X_c(blk_idx, blk_idx) = Xi_c;   
end

%% evaluate
fprintf('The relative recovery error is: %d\n', norm(X-X_c,'fro')/norm(X,'fro'))

%% compare with full matrix completion
Omega = find(Omega==1);
data = X(Omega);
[U,S,V,numiter] = SVT(sz,Omega,data,tau,delta,maxiter,tol);
X_c2 = U*S*V'
fprintf('The relative recovery error is: %d\n', norm(X-X_c2,'fro')/norm(X,'fro'))
