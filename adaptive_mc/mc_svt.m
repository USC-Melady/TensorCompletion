function [Xhat Obs] = mc_svt(X, p),
%% The SVT algorithm of Cai et al. 
  [d n] = size(X);
  tau = 5*max(d,n);
  kmax = 500;

  inds = unique(randi(n*d, [int32(n*d*p), 1]));
  obs = zeros(d*n,1);
  obs(inds) = 1;
  Obs = reshape(obs, [d,n]);

  delta = 1.2/p;
  Msub = X.*Obs;
  Xrep = reshape(X, [d*n,1]);
  
  if norm(Msub,2) ~= 0,
      k0 = tau/(delta*norm(Msub, 2));
  else,
      Xhat = zeros(d,n);
      return;
  end;
  Xk = zeros(d,n);
  Yk = k0*delta*Msub;
  for k=1:kmax,
      %% fprintf('k=%d\n', k);
      [u,s,v] = svd(Yk, 'econ');
      s = diag(max(diag(s)-tau,0));
      Xk = u*s*v';
      Xkrep = reshape(Xk, [d*n,1]);
      if norm(Xrep(inds) - Xkrep(inds), 'fro')/norm(Msub, 'fro') <= 0.0001,
          break;
      end;
      if norm(Xrep(inds) - Xkrep(inds), 'fro')/norm(Msub, 'fro') >= 1000,
          Xhat = zeros(d,n);
          return;
      end;
      %% if mod(k, 100) == 0,
          %% fprintf('%d %5.5f\n', k, norm(Xrep(inds) - Xkrep(inds), ...
          %% 'fro')/norm(Msub, 'fro'));
          %% end;
      Yk = Yk + delta*(reshape(Xrep - Xkrep,[d,n])).*Obs;
  end;
  Xhat = Xk;
  