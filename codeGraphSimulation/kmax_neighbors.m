% Finds the number of "kmin"-neighbors (k links away at a minimum) for every node
% If nodes are k-links away due to loops (so they appear as m-neighbours, m<k), they are not counted
% INPUTS: adjacency matrix, node index, k - number of links
% OUTPUTS: vector of "kmin"-neighbors indices
% GB, May 16, 2011

function kneigh = kmax_neighbors(adj,ind,k)

close_neighbors=[];

adjk = adj;
for i=1:k
  close_neighbors = [close_neighbors find(adjk(ind,:)>0)];
  adjk = adjk*adj;
end

% kneigh = setdiff(find(adjk(ind,:)>0),[close_neighbors ind]);
kneigh = unique([close_neighbors ind]);