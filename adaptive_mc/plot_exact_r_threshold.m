%% Generate the plots of the success probability for the adaptive
%% MC algorithm as we vary the rank r. 

load('./data/exact_r_threshold.mat');

f = figure();
l1 = plot(ps, r10, 'b--o', 'MarkerSize', 8, 'LineWidth', 2);
hold on;
l2 = plot(ps, r20, 'g-o', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'g', 'LineWidth', 2);
l3 = plot(ps, r30, 'r-d', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'r', 'LineWidth', 2);
l4 = plot(ps, r40, 'k-s', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'k', 'LineWidth', 2);

xlim([0, 0.6]);
h = xlabel('Fraction of Samples/Column (p)');
set(h, 'FontSize', 20);
h = ylabel('Probability of Recovery');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);

legend([l1 l2 l3 l4], 'r=10', 'r=20', 'r=30', 'r=40', 4);

saveas(f, './figs/exact_r_threshold_1.fig', 'fig');
saveas(f, './figs/exact_r_threshold_1.eps', 'epsc');

f = figure();
l1 = plot(ps/(10*log(10)), r10, 'b--o', 'MarkerSize', 8, 'LineWidth', 2);
hold on;
l2 = plot(ps/(20*log(20)), r20, 'g-o', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'g', 'LineWidth', 2);
l3 = plot(ps/(30*log(30)), r30, 'r-d', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'r', 'LineWidth', 2);
l4 = plot(ps/(40*log(40)), r40, 'k-s', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'k', 'LineWidth', 2);

xlim([0, 0.015]);
h = xlabel('Rescaled Sampling Probability (p/(r log r))');
set(h, 'FontSize', 20);
h = ylabel('Probability of Recovery');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);

legend([l1 l2 l3 l4], 'r=10', 'r=20', 'r=30', 'r=40', 4);

saveas(f, './figs/exact_r_threshold_2.fig', 'fig');
saveas(f, './figs/exact_r_threshold_2.eps', 'epsc');
