%% Generate the plots of the success probability for the SVT
%% algorithm. 

load('./data/exact_svt_threshold.mat');

f = figure();
l1 = plot(100*ps/log(100), m100, 'b--o', 'MarkerSize', 8, 'LineWidth', 2);
hold on;
l2 = plot(150*ps/log(150), m150, 'g-s', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'g', 'LineWidth', 2);
l3 = plot(200*ps/log(200), m200, 'r-d', 'MarkerSize', 8, 'MarkerFaceColor', ...
          'r', 'LineWidth', 2);


%% xlim([0, 150]);
h = xlabel('Rescaled Sampling Probability (np/log(n))');
set(h, 'FontSize', 20);
h = ylabel('Probability of Recovery');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);

legend([l1 l2 l3], 'n=100', 'n=150', 'n=200', 4);

saveas(f, './figs/exact_svt_threshold_2.fig', 'fig');
saveas(f, './figs/exact_svt_threshold_2.eps', 'epsc');
