%% Generate success probability plots for the adaptive mc algorithm
%% as we vary the problem size n.

load('./data/exact_p_threshold.mat');

f = figure();
l1 = plot(ps, m250, 'b--o', 'MarkerSize', 8, 'LineWidth', 2);
hold on;
l2 = plot(ps, m500, 'g-o', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'g', 'LineWidth', 2);
l3 = plot(ps, m750, 'r-d', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'r', 'LineWidth', 2);
l4 = plot(ps, m1000, 'k-s', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'k', 'LineWidth', 2);

h = xlabel('Fraction of Samples/Column (p)');
set(h, 'FontSize', 20);
h = ylabel('Probability of Recovery');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);

legend([l1 l2 l3 l4], 'n=250', 'n=500', 'n=750', 'n=1000', 4);

saveas(f, './figs/exact_p_threshold_1.fig', 'fig');
saveas(f, './figs/exact_p_threshold_1.eps', 'epsc');

f = figure();
l1 = plot(250*ps, m250, 'b--o', 'MarkerSize', 8, 'LineWidth', 2);
hold on;
l2 = plot(500*ps, m500, 'g-o', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'g', 'LineWidth', 2);
l3 = plot(750*ps, m750, 'r-d', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'r', 'LineWidth', 2);
l4 = plot(1000*ps, m1000, 'k-s', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'k', 'LineWidth', 2);

xlim([0, 150]);
h = xlabel('Number of Samples/Column (np)');
set(h, 'FontSize', 20);
h = ylabel('Probability of Recovery');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);

legend([l1 l2 l3 l4], 'n=250', 'n=500', 'n=750', 'n=1000', 4);

saveas(f, './figs/exact_p_threshold_2.fig', 'fig');
saveas(f, './figs/exact_p_threshold_2.eps', 'epsc');
