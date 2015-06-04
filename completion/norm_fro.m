function [ norm_val ] = norm_fro( X )
%NORM_FRO Summary of this function goes here
%   Detailed explanation goes here


R = size(X,3);
norm_val = 0;
for r = 1:R
    norm_val =  norm_val +  norm(X(:,:,r),'fro')^2;
end

norm_val = sqrt(norm_val);
end

