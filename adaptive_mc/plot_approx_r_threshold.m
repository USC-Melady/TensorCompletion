%% Generate error plots for adaptive matrix approximation algorithm
%% as we vary rank r. 

load('./data/approx_r_threshold.mat');

f = figure();
l1 = plot(ps, m10, 'b-o', 'MarkerSize', 4, 'LineWidth', 1);
hold on;
l2 = plot(ps, m20, 'g-s', 'MarkerSize', 4, 'MarkerFaceColor', ...
          'g', 'LineWidth', 1);
l3 = plot(ps, m30, 'r-d', 'MarkerSize', 4, 'MarkerFaceColor', ...
          'r', 'LineWidth', 1);

h = xlabel('Fraction of Samples/Column (p)');
set(h, 'FontSize', 20);
h = ylabel('Relative Error (epsilon)');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);

legend([l1 l2 l3], 'r=10', 'r=20', 'r=30', 1);

saveas(f, './figs/approx_r_threshold_1.fig', 'fig');
saveas(f, './figs/approx_r_threshold_1.eps', 'epsc');

f = figure();
l1 = plot(ps, m10/10^(1/2), 'b-o', 'MarkerSize', 4, 'LineWidth', 1);
hold on;
l2 = plot(ps, m20/20^(1/2), 'g-s', 'MarkerSize', 4, 'MarkerFaceColor', ...
          'g', 'LineWidth', 1);
l3 = plot(ps, m30/30^(1/2), 'r-d', 'MarkerSize', 4, 'MarkerFaceColor', ...
          'r', 'LineWidth', 1);

h = xlabel('Fraction of Samples/Column (p)');
set(h, 'FontSize', 20);
h = ylabel('Rescaled Relative Error (epsilon/sqrt(r))');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);

legend([l1 l2 l3], 'r=10', 'r=20', 'r=30', 1);

saveas(f, './figs/approx_r_threshold_2.fig', 'fig');
saveas(f, './figs/approx_r_threshold_2.eps', 'epsc');
