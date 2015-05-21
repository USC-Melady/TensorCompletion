%% Generate error plots for adaptive matrix approximation algorithm
%% as we vary problem size n

load('./data/approx_p_threshold.mat');

f = figure();
l1 = plot(ps, m250, 'b--o', 'MarkerSize', 4, 'LineWidth', 1);
hold on;
l2 = plot(ps, m500, 'g-o', 'MarkerSize', 4, 'MarkerFaceColor', ...
          'g', 'LineWidth', 1);
l3 = plot(ps, m750, 'r-d', 'MarkerSize', 4, 'MarkerFaceColor', ...
          'r', 'LineWidth', 1);
l4 = plot(ps, m1000, 'k-s', 'MarkerSize', 4, 'MarkerFaceColor', ...
          'k', 'LineWidth', 1);

h = xlabel('Fraction of Samples/Column (p)');
set(h, 'FontSize', 20);
h = ylabel('Relative Error (epsilon)');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);

legend([l1 l2 l3 l4], 'n=250', 'n=500', 'n=750', 'n=1000', 1);

saveas(f, './figs/approx_p_threshold_1.fig', 'fig');
saveas(f, './figs/approx_p_threshold_1.eps', 'epsc');

f = figure();
l1 = plot(250/log(250)*ps, m250, 'b--o', 'MarkerSize', 4, 'LineWidth', 1);
hold on;
l2 = plot(500/log(500)*ps, m500, 'g-o', 'MarkerSize', 4, 'MarkerFaceColor', ...
          'g', 'LineWidth', 1);
l3 = plot(750/log(750)*ps, m750, 'r-d', 'MarkerSize', 4, 'MarkerFaceColor', ...
          'r', 'LineWidth', 1);
l4 = plot(1000/log(1000)*ps, m1000, 'k-s', 'MarkerSize', 4, 'MarkerFaceColor', ...
          'k', 'LineWidth', 1);

h = xlabel('Rescaled Samples/Column (np/log(n))');
set(h, 'FontSize', 20);
h = ylabel('Relative Error (epsilon)');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);

legend([l1 l2 l3 l4], 'n=250', 'n=500', 'n=750', 'n=1000', 1);

saveas(f, './figs/approx_p_threshold_2.fig', 'fig');
saveas(f, './figs/approx_p_threshold_2.eps', 'epsc');

f = figure();
l1 = plot(ps, sqrt(ps).*m250, 'b--o', 'MarkerSize', 4, 'LineWidth', 1);
hold on;
l2 = plot(ps, sqrt(ps).*m500, 'g-o', 'MarkerSize', 4, 'MarkerFaceColor', ...
          'g', 'LineWidth', 1);
l3 = plot(ps, sqrt(ps).*m750, 'r-d', 'MarkerSize', 4, 'MarkerFaceColor', ...
          'r', 'LineWidth', 1);
l4 = plot(ps, sqrt(ps).*m1000, 'k-s', 'MarkerSize', 4, 'MarkerFaceColor', ...
          'k', 'LineWidth', 1);

h = xlabel('Fraction of Samples/Column (p)');
set(h, 'FontSize', 20);
h = ylabel('Rescaled Relative Error (sqrt(p)*epsilon)');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);

legend([l1 l2 l3 l4], 'n=250', 'n=500', 'n=750', 'n=1000', 1);

saveas(f, './figs/approx_p_threshold_3.fig', 'fig');
saveas(f, './figs/approx_p_threshold_3.eps', 'epsc');

f = figure();
l1 = plot(ps, ps.*(m250.^1), 'b--o', 'MarkerSize', 4, 'LineWidth', 1);
hold on;
l2 = plot(ps, ps.*(m500.^1), 'g-o', 'MarkerSize', 4, 'MarkerFaceColor', ...
          'g', 'LineWidth', 1);
l3 = plot(ps, ps.*(m750.^1), 'r-d', 'MarkerSize', 4, 'MarkerFaceColor', ...
          'r', 'LineWidth', 1);
l4 = plot(ps, ps.*(m1000.^1), 'k-s', 'MarkerSize', 4, 'MarkerFaceColor', ...
          'k', 'LineWidth', 1);

h = xlabel('Fraction of Samples/Column (p)');
set(h, 'FontSize', 20);
h = ylabel('Rescaled Relative Error (p*epsilon)');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);

legend([l1 l2 l3 l4], 'n=250', 'n=500', 'n=750', 'n=1000', 1);

saveas(f, './figs/approx_p_threshold_4.fig', 'fig');
saveas(f, './figs/approx_p_threshold_4.eps', 'epsc');
