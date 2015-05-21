%% Plot to compare soft thresholding and adaptive matrix
%% approximation on flat and spiky matrices. 

load('./data/approx_spiky_comparison.mat');

f = figure();

l1 = plot(ps, adapt, 'b-o', 'MarkerSize', 4, 'LineWidth', 1);
hold on;
l2 = plot(ps, soft, 'g-s', 'MarkerSize', 4, 'LineWidth', 1);

h = xlabel('Fraction of Samples/Column (p)');
set(h, 'FontSize', 20);
h = ylabel('Relative Error (epsilon)');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);

legend([l1 l2], 'Adapt', 'Passive');
saveas(f, './figs/approx_spiky_comparison.fig', 'fig');
saveas(f, './figs/approx_spiky_comparison.eps', 'epsc');

load('./data/approx_flat_comparison.mat');

f = figure();

l1 = plot(ps, adapt, 'b-o', 'MarkerSize', 4, 'LineWidth', 1);
hold on;
l2 = plot(ps, soft, 'g-s', 'MarkerSize', 4, 'LineWidth', 1);

h = xlabel('Fraction of Samples/Column (p)');
set(h, 'FontSize', 20);
h = ylabel('Relative Error (epsilon)');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);

legend([l1 l2], 'Adapt', 'Passive');
saveas(f, './figs/approx_flat_comparison.fig', 'fig');
saveas(f, './figs/approx_flat_comparison.eps', 'epsc');

