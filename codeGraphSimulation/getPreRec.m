function [ pre, rec ] = getPreRec( n,N,kneigh )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = N/2;
kneigh(kneigh>N) = kneigh(kneigh>N)-N;

hit = numel(find(kneigh<=n));

pre = hit/numel(kneigh);
rec = hit/2/n;

end

