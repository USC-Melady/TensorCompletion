function [ G, subIndexes ] = generateGraph( n,rho,p )
%GENERATEGRAPH Summary of this function goes here
%   Detailed explanation goes here

nOverlap = ceil(n*rho);
nSize = 2*n-nOverlap;

M = zeros(nSize);

tmp = rand(n);
tmp(tmp>p)=0;
M(1:n,1:n) = M(1:n,1:n) +tmp;

tmp = rand(n);
tmp(tmp>p)=0;
M(end-n+1:end,end-n+1:end) = M(end-n+1:end,end-n+1:end) +tmp;
M = sign(M);

G = [zeros(nSize), M ;M',zeros(nSize)];
subIndexes =cell(2,2);
subIndexes{1,1} = 1:n;
subIndexes{1,2} = 1:n;

subIndexes{1,1} = (nSize-n+1):nSize;
subIndexes{1,2} = (nSize-n+1):nSize;

end

