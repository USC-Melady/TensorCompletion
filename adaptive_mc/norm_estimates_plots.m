%% Generate plot of column norms and our estimates for them.

M = incoherent_approx_matrix(100, 100, 5, 1);
f = figure();
imagesc(M);
saveas(f, './figs/incoherent_approx_matrix.eps', 'epsc');


t_probs = diag(M'*M)/sum(diag(M'*M));

[Xhat, Obs1, e_probs1] = mc_approx_adapt(M, 0.05, 20);
[Xhat, Obs2, e_probs2] = mc_approx_adapt(M, 0.2, 20);

f = figure();
l1 = plot(t_probs, 'blue');
hold on;
l2 = plot(e_probs1, 'red');
l3 = plot(e_probs2, 'green');

legend([l1,l2,l3], 'Column Norms', 'Estimate p=0.05', ['Estimate ' ...
                    'p=0.2']);
saveas(f, './figs/norm_estimates.eps', 'epsc');

f = figure();
imagesc(Obs2);
saveas(f, './figs/sampling_pattern.eps', 'epsc');