
%% get F1

results_f11 = 2.*results_pre1.*results_rec1./(results_pre1+results_rec1);
results_f12 = 2.*results_pre2.*results_rec2./(results_pre2+results_rec2);


%% average through runs

results_f11 = mean(results_f11,6);
results_f12 = mean(results_f12,6);


%%