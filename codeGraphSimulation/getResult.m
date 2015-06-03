function [ pre1,rec1,pre2,rec2 ] = getResult( G, d, r ,n, nRuns)
%GETRESULT Summary of this function goes here
%   Detailed explanation goes here

N = size(G,1);
[Is,~] = find(G);
m = numel(Is);

pre1 = zeros(nRuns,1);
rec1 = zeros(nRuns,1);
pre2 = zeros(nRuns,1);
rec2 = zeros(nRuns,1);


for irun = 1:nRuns
    % sample a vertex from the first cluster.
    ms = randperm(m,1);
    while Is(ms)>n
        ms = randperm(m,1);
    end
    v = Is(ms);
    kneigh = kmax_neighbors(G,v,d);    
    [ pre1(irun), rec1(irun) ] = getPreRec( n,N,kneigh );
      
    % Skrinking
    % find subgraph
    Gs = G(kneigh,kneigh);
    
    keep = skrinkR( Gs,r );
    keep = kneigh(keep);
    [ pre2(irun), rec2(irun) ] = getPreRec( n,N,keep );
end
end

