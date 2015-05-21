%% Generate probabity of success plotsfor the adaptive mc algorithm
%% as we vary mu. 
load('./data/exact_mu_threshold.mat');

f = figure();
l1 = plot(ps, mu1, 'b--o', 'MarkerSize', 8, 'LineWidth', 2);
hold on;
l2 = plot(ps, mu2, 'g-o', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'g', 'LineWidth', 2);
l3 = plot(ps, mu3, 'r-d', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'r', 'LineWidth', 2);
l4 = plot(ps, mu4, 'k-s', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'k', 'LineWidth', 2);

xlim([0, 0.6]);
h = xlabel('Fraction of Samples/Column (p)');
set(h, 'FontSize', 20);
h = ylabel('Probability of Recovery');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);

legend([l1 l2 l3 l4], 'mu=1', 'mu=2', 'mu=3', 'mu=4', 4);

saveas(f, './figs/exact_mu_threshold_1.fig', 'fig');
saveas(f, './figs/exact_mu_threshold_1.eps', 'epsc');

f = figure();
l1 = plot(ps/1, mu1, 'b--o', 'MarkerSize', 8, 'LineWidth', 2);
hold on;
l2 = plot(ps/2, mu2, 'g-o', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'g', 'LineWidth', 2);
l3 = plot(ps/3, mu3, 'r-d', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'r', 'LineWidth', 2);
l4 = plot(ps/4, mu4, 'k-s', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'k', 'LineWidth', 2);

xlim([0, 0.2]);
h = xlabel('Rescaled Sampling Probability (p/mu)');
set(h, 'FontSize', 20);
h = ylabel('Probability of Recovery');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);

legend([l1 l2 l3 l4], 'mu=1', 'mu=2', 'mu=3', 'mu=4', 4);

saveas(f, './figs/exact_mu_threshold_2.fig', 'fig');
saveas(f, './figs/exact_mu_threshold_2.eps', 'epsc');
