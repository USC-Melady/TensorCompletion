
ns = 50:200:1000;
rhos = 0.01:0.05:0.2;
rs = 1:2:8;
ds = 2:6;
cs = [0.5 1 2];

nRuns = 50;
%%
results_pre1 = zeros(numel(ns),numel(rhos),numel(rs),numel(ds),numel(cs),nRuns);
results_rec1 = zeros(numel(ns),numel(rhos),numel(rs),numel(ds),numel(cs),nRuns);
results_pre2 = zeros(numel(ns),numel(rhos),numel(rs),numel(ds),numel(cs),nRuns);
results_rec2 = zeros(numel(ns),numel(rhos),numel(rs),numel(ds),numel(cs),nRuns);
%%
for di = 1:numel(ds)
    for rhoi = 1:numel(rhos)
        for ri = 1:numel(rs)
            for ci = 1:numel(cs)
                for ni = 1:numel(ns)
                    n = ns(ni);
                    r = rs(ri);
                    d = ds(di);
                    c = cs(ci)*log(n);
                    p = c* r /n;
                    rho = rhos(rhoi);                    
                    fprintf('n:%d\t r:%d\t d:%d\t p:%f\t rho:%f\n',n,r,d,p,rho);
                    
                    parfor irun = 1:nRuns
%                         fprintf('n:%d\t r:%d\t d:%d\t p:%f\t rho:%f irun:%d\n',n,r,d,p,rho,irun);
                        [ G, subIndexes ] = generateGraph( n,rho,p );                        
                        [   results_pre1(ni,rhoi,ri,di,ci,irun),...
                            results_rec1(ni,rhoi,ri,di,ci,irun),...
                            results_pre2(ni,rhoi,ri,di,ci,irun),...
                            results_rec2(ni,rhoi,ri,di,ci,irun)...
                            ] = getResult( G, d, r ,n, 1);
                    end
                    
                    
                end
            end
        end
    end
end

save('results','results_pre1','results_rec1','results_pre2','results_rec2');
