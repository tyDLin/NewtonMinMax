%%************************************************************************* 
%% Call different algorithms to solve AUC maximization problems
%%*************************************************************************

%%
clear;
clc; 
close all; 

addpath('solver')
addpath('subroutine')

ranseed = 1;
rng(ranseed, 'twister');

%% Problem setting.  
ndatasets = {'a9a', 'covtype', 'w8a'};
init      = [0.05 0.05 0.5];
x_epoch   = [20 20 20]; 
x_time    = [7 85 40]; 
yaxis     = [1 1 0.3];
T         = 200; 

for di=1:length(ndatasets) 
        
    %% generate data. 
    dataset_name = ndatasets{di};
    
    fprintf('\nProcessing (%d/%d) dataset: %s\n', di, length(ndatasets), dataset_name);
    data_path = './data/';
    load([data_path dataset_name '.mat']);
    
    A = full(samples);
    [n, d] = size(A); 
    b = labels;
    rho = 1/n; 
    eta = init(di); 
    
    %% generate optimal solution
    z_star = centroid_opt(A, b, rho);

    %% call SSNDE
    optsSSNDE.max_iters      = 10;
    optsSSNDE.display        = 0;
    optsSSNDE.displayfreq    = 1;
    optsSSNDE.checkfreq      = 1;
    optsSSNDE.savegaphist    = 1;
    optsSSNDE.savesamplehist = 1;
    optsSSNDE.savetimehist   = 1;
    [~, gaphist_SSNDE, iterhist_SSNDE, timehist_SSNDE] = centroid_SSNDE(z_star, A, b, rho, eta, optsSSNDE);

    %% call SEG
    optsSEG.max_iters     = T;
    optsSEG.display       = 0;
    optsSEG.displayfreq   = 100;
    options.checkfreq     = 100; 
    optsSEG.savegaphist   = 1;
    optsSEG.saveiterhist  = 1;
    optsSEG.savetimehist  = 1;
    [~, gaphist_SEG, iterhist_SEG, timehist_SEG] = centroid_SEG(z_star, A, b, rho, optsSEG);
 
    %% call SVREG
    optsSVREG.max_iters     = T;
    optsSVREG.display       = 0;
    optsSVREG.displayfreq   = 100;
    options.checkfreq       = 100; 
    optsSVREG.savegaphist   = 1;
    optsSVREG.saveiterhist  = 1;
    optsSVREG.savetimehist  = 1;
    [~, gaphist_SVREG, iterhist_SVREG, timehist_SVREG] = centroid_SVREG(z_star, A, b, rho, optsSVREG);
 
    %% plot the figures
    round_iter = 1:(T/10):T;

    figure; 
    plot(iterhist_SEG(round_iter), gaphist_SEG(round_iter), '-*', 'LineWidth', 3, 'MarkerSize', 15);
    hold on
    plot(iterhist_SVREG(round_iter), gaphist_SVREG(round_iter), '-d', 'LineWidth', 3, 'MarkerSize', 15);
    hold on
    plot(iterhist_SSNDE, gaphist_SSNDE, '-v', 'LineWidth', 3, 'MarkerSize', 15);
    hold off
    legend({'SEG', 'SVREG', 'Our Method'}, 'Location', 'northwest', ...
        'Orientation', 'vertical', 'NumColumns', 3);
        
    set(gca, 'YScale','log');
    set(gca, 'FontSize', 20);
    xlabel('Epoch Count');
    ylabel('Gap Function');
    xlim([0 x_epoch(di)])
    ylim([0 yaxis(di)])
    title([dataset_name ' (N=', num2str(n), ', n=', num2str(d), ')']);

    path = sprintf('../figs/auc_%s_epoch', dataset_name); 
    saveas(gcf, path, 'epsc');
    
    round_time = 1:(T/10):T;
 
    figure; 
    plot(timehist_SEG(round_time), gaphist_SEG(round_time), '-*', 'LineWidth', 3, 'MarkerSize', 15);
    hold on
    plot(timehist_SVREG(round_time), gaphist_SVREG(round_time), '-d', 'LineWidth', 3, 'MarkerSize', 15);
    hold on
    plot(timehist_SSNDE, gaphist_SSNDE, '-v', 'LineWidth', 3, 'MarkerSize', 15);
    hold off
    legend({'SEG', 'SVREG', 'Our Method'}, 'Location', 'northwest', ...
        'Orientation', 'vertical', 'NumColumns', 3);
         
    set(gca, 'YScale','log');
    set(gca, 'FontSize', 20);
    xlabel('Time (Seconds)');
    ylabel('Gap Function');
    xlim([0 x_time(di)])
    ylim([0 yaxis(di)])
    title([dataset_name ' (N=', num2str(n), ', n=', num2str(d), ')']);
 
    path = sprintf('../figs/auc_%s_time', dataset_name); 
    saveas(gcf, path, 'epsc');
end