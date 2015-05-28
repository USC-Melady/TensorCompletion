function [submat_idx] = find_comp( n, G , seed, r, d)
%FIND_COMP Summary of this function goes here
%   Detailed explanation goes here
% G: adjacent matrix of a graph
% seed: vertex index for the seed
% r: radius
% d: depth

n1 = n(1);
n2 = n(2);
% perform  BFS
[dist ,disc_time, predecessor]  = bfs(G,seed);
% construct subgraph
G_sub_ind  = find(dist <=d);
G_full = full(G);
G_sub = G_full(G_sub_ind, G_sub_ind);

% find the maximal r-conected subgraph
[K,C] = find_comp_simple(G_sub, r);


%%

%% 
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


end

