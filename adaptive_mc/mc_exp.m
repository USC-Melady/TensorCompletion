function [ms vs] = mc_exp(n1, n2, incs, samples)
%% Generate error rate for sdp algorithm. 
sigma = 0.1;
ms = [];
vs = [];
for i=1:length(incs),
    sube = [];
    M = incoherent_matrix(n1,n2,incs(i));
    for j=1:5,
        R = sigma/sqrt(n1*n2)*normrnd(0, 1, [n1 n2]);
        samples(i)/(n1*n2)
        B = binornd(1, samples(i)/(n1*n2), [n1 n2])
        X = mc_sdp((M+R).*B, sigma^4, norm(vec(M),Inf));
        error = square_pos(norm(X - M));
        sube = [sube error];
    end;
    ms = [ms mean(sube)];
    vs = [vs var(sube)];
end;