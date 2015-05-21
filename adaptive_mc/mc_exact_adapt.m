function [Xhat Obs] = mc_exact_adapt(X,p)
%% Adaptive matrix completion algorithm from NIPS paper and JMLR
%% submission. 
    tau = 0.00001;
    d = size(X,1);
    n = size(X,2);

    Xhat = zeros(d,n);

    Obs = zeros(d,n);
    
    U = 0;

    inds = unique(randi(d, [int32(d*p), 1]));
    for i=1:n,
        xto = X(inds,i);
        Obs(inds,i) = 1;
        if U == 0,
            if norm(xto) > tau,
                U = X(:,i);
                U = orthogonalize(U);
                Obs(:,i) = 1;
                Xhat(:,i) = X(:,i);
                inds = unique(randi(d, [int32(d*p), 1]));
                Usub = U(inds,:);
                if all(Usub == 0),
                    Usuborth = zeros(size(Usub));
                    Usubinv = zeros(size(Usub,2), size(Usub,2));
                else,
                    Usuborth = orthogonalize(Usub);
                    Usubinv = (Usub'*Usub)^(-1);
                end;
            else,
                Xhat(:,i) = 0;
            end;
        else,
            if norm(xto - Usuborth*Usuborth'*xto) > tau,
                U = [U X(:,i)];
                U = orthogonalize(U);
                Obs(:,i) = 1;
                Xhat(:,i) = X(:,i);
                inds = unique(randi(d, [int32(d*p), 1]));
                Usub = U(inds,:);
                if all(Usub == 0),
                    Usuborth = zeros(size(Usub));
                    Usubinv = zeros(size(Usub,2), size(Usub, 2));
                else,
                    Usuborth = orthogonalize(Usub);
                    if size(Usuborth,2) ~= size(Usub,2),
                        return;
                    end;
                    Usubinv = (Usub'*Usub)^(-1);
                end;
            else,
                Xhat(:,i) = U*Usubinv*Usub'*xto;
            end;
        end;
    end;

function [U] = orthogonalize(U),
   [u,s,v] = svd(U);
   U = u(:,diag(s)>0.001);
