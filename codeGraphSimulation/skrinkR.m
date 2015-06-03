function [ keep ] = skrinkR( G,r )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
keep = 1:size(G);
deg = sum(G);
[~,~,vdeg] = find(deg);
while min(vdeg)< r
    keep = find(deg>=r);
    Gt = zeros(size(G));
    Gt(keep,keep) = G(keep,keep);
    G = Gt;
    deg = sum(G);
    [~,~,vdeg] = find(deg);
end

end

