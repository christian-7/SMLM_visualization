% Show fitting functions

load('K:\Christian\GitHub\SMLM_vis\exp_dist\dist_gap_fit_exp.mat')

x1 = 1:0.5:80;
y1 = distfit(x1);

x2 = 1:1:500;
y2 = gap_fit_exp(x2);

figure('Position',[100 200 800 300])

subplot(1,2,1)
scatter(x1,y1,4,'b','filled')
title('Localization precision');
xlabel('nm');
ylabel('probability');
box on
axis square


subplot(1,2,2)
scatter(x2,y2,10,'b','filled')
title('Dark time');
xlabel('frames');
ylabel('probability');
box on
axis square

%% Histogram of GT along the cross section

load('K:\Christian\GitHub\SMLM_vis\simulated_test_data\simulated_MT_3D_radius_20nm_030_GT.mat')

bins = -20:5:50;
y = hist(mol_list2(:,2),bins);

figure
scatter(bins-20, y); hold on;
plot(bins-20, y);
title('Cross section');
xlabel('nm');
ylabel('counts');
box on
axis square
