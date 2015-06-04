function [ W ,tensorW,Ls,train_time] = MLMTL_Convex( X, Y, indicators, beta, lambda, outerNiTPre, thresholdPre, groundW )
%MLMTL Summary of this function goes here
%   Detailed explanation goes here
global verbose
outerNiT=1000;
if nargin>5 && ~isempty(outerNiTPre)
    outerNiT=outerNiTPre;
end
threshold = 10e-5;
if nargin>6 && ~isempty(thresholdPre)
    threshold=thresholdPre;
end

tic;
nTotalTasks=length(Y);
nAttrs=getNAttrs(X);
nModes=length(indicators);

if nTotalTasks~=prod(indicators(2:end))
    [nTotalTasks prod(indicators(2:end))]
    error('There are inconsistencies between the indicators and the number of tasks.');
end

XX_plus_betaNI=cell(1,nTotalTasks);
XY=cell(1,nTotalTasks);
for t=1:nTotalTasks
    if isempty(X{t})
        XX_plus_betaNI{t}=beta*nModes*eye(nAttrs);
        XY{t}=0;
    else
        XX_plus_betaNI{t}=1/lambda*(X{t}*X{t}')+beta*nModes*eye(nAttrs);
        XY{t}=X{t}*Y{t};
    end
end

A=cell(1,nModes);
B=cell(1,nModes);
for n=1:nModes
    A{n}=tenzeros(indicators);
    B{n}=tenzeros(indicators);
end
sumA=tenzeros(indicators); 
sumB=tenzeros(indicators); 
oldW=Inf(nAttrs, nTotalTasks);
oit=0;

Ls = [];
oldL =0;
while true    
    oit=oit+1;
    % Optimizing over X
    matSumAux=tenmat(sumA, 1) + beta*tenmat(sumB, 1);
    matSumAux=matSumAux.data;
    matW=zeros(nAttrs, nTotalTasks);
    for t=1:nTotalTasks
        matW(:,t)=(XX_plus_betaNI{t})\(1/lambda*XY{t} + matSumAux(:,t));
    end
    W=tensor(reshape(matW, indicators));
   
    sumA=tenzeros(indicators); 
    sumB=tenzeros(indicators); 
    for n=1:nModes
%         [oit n]
        % Optimizing over B
        matW=tenmat(W, n);
        matW=matW.data;
        matA=tenmat(A{n},n);
        matA=matA.data;
        matB=shrink(matW-1/beta*matA, 1/beta);
        BTensor=tensor(reshape(matB, [indicators(n), indicators(1:n-1), indicators(n+1:end)]));
        B{n}=permute(BTensor, [2:n, 1, n+1:nModes]);
        sumB=sumB+B{n};
        
        % Optimizing over A
        A{n}=A{n}-beta*(W-B{n});
        sumA=sumA+A{n};
    end
    
    Wmat=tenmat(full(W), 1);
    Wmat=Wmat.data;
    if nargin>7 && ~isempty(groundW)
        disp(['RSE=' num2str(norm(Wmat(1:end)-groundW(1:end))/norm(groundW(1:end)))]);
    end
    if oit>outerNiT %norm(Wmat(1:end)-oldW(1:end))<threshold
        break
    end
    L = MLMTL_Objfunc(X,Y,W,B,A,lambda,beta);
    Ls = [Ls,L];
    
%     if( abs(L-oldL)< threshold)
    if(norm(Wmat(1:end)-oldW(1:end)) < threshold)
%         if verbose
%             fprintf('MLMTL_Convex:Converge after %d iteration \n', oit);
%         end
        break;
    end
    oldW=Wmat;
    oldL=L;
end
% plot(Ls);
% disp('L_Inf');
for i=1:nModes
    mat=tenmat(W, i);
    [u l v]=mySVD(mat.data);
    max(diag(l));
end
% norm(mat.data, 'fro')

tensorW=W;
W=tenmat(full(W), 1);
W=W.data;
train_time = toc;
end

