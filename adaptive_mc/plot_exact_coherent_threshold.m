%% Generate success probability plots for adaptive MC vs SVT on
%% cmatrices with coherent row space. 

load('./data/exact_coherent_threshold.mat')

f = figure();
l1 = plot(ps, a1, 'b--o', 'MarkerSize', 8, 'LineWidth', 2);
hold on;
l2 = plot(ps, a2, 'g-o', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'g', 'LineWidth', 2);
l3 = plot(ps, s1, 'r-d', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'r', 'LineWidth', 2);
l4 = plot(ps, s2, 'k-s', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'k', 'LineWidth', 2);

h = xlabel('Fraction of Samples/Column (p)');
set(h, 'FontSize', 20);
h = ylabel('Probability of Recovery');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);

legend([l1 l2 l3 l4], 'Adapt mu=1', 'Adapt mu=2', 'SVT mu=1', 'SVT mu=2', 4);

saveas(f, './figs/exact_coherent_threshold.fig', 'fig');
saveas(f, './figs/exact_coherent_threshold.eps', 'epsc');
