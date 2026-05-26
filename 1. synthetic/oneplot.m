%%
clear;
clc; 
close all; 

%% Problem setting: generate A
ndim        = [20 50 100 200 500 1000];
T           = 50000;
x_iter      = [50 50 100 400 400 400]; 
x_time      = [0.06 0.3 2.5 30 150 300]; 
yaxis_upper = [1e+4 5e+4 1e+5 1e+6 1e+6 1e+6]; 
yaxis_lower = [1e-11 0 0 0 0 0];

for di=1:length(ndim)
        
    %% load results. 
    n = ndim(di);
    file_name_alg = sprintf('../results/cubic_%d.mat', n);
    load(file_name_alg);
    
    round_EG = 1:T/100:T; 
    round_OGDA = 1:T/50:2*T; 
    round_iter = 1:x_iter(di)/10:x_iter(di)+1;

    if n == 200
        round_EG = 1:T/20:T; 
        round_OGDA = 1:T/10:2*T; 
    end

    %% plot the figures
    figure; 
    plot(round_iter, gaphist_EG1(round_iter), '-p', 'LineWidth', 3, 'MarkerSize', 10);
    hold on
    plot(round_iter, gaphist_OGDA1(round_iter), '-d', 'LineWidth', 3, 'MarkerSize', 10);
    hold on
    plot(1:length(gaphist_NDE), gaphist_NDE, '-v', 'LineWidth', 3, 'MarkerSize', 10);
    hold on
    plot(counthist_EG2, gaphist_EG2, '-o', 'LineWidth', 3, 'MarkerSize', 10);
    hold on
    plot(counthist_OGDA2, gaphist_OGDA2, '-s', 'LineWidth', 3, 'MarkerSize', 10);
    hold off
    legend({'EG1', 'OGDA1', 'Our Method', 'EG2', 'OGDA2'}, 'Location', 'northwest', ...
        'Orientation', 'horizontal', 'NumColumns', 3);
    
    set(gca, 'YScale','log');
    set(gca, 'FontSize', 15);
    xlabel('Iteration Count');
    ylabel('Gap Function');
    xlim([0 x_iter(di)])
    ylim([0 yaxis_upper(di)])
    title(['n=', num2str(n)]);

    path = sprintf('../figs/cubic_%d_iteration', n); 
    saveas(gcf, path, 'epsc');

    figure; 
    plot(timehist_EG1(round_EG), gaphist_EG1(round_EG), '-p', 'LineWidth', 3, 'MarkerSize', 10);
    hold on
    plot(timehist_OGDA1(round_OGDA), gaphist_OGDA1(round_OGDA), '-d', 'LineWidth', 3, 'MarkerSize', 10);
    hold on
    plot(timehist_NDE, gaphist_NDE, '-v', 'LineWidth', 3, 'MarkerSize', 10);
    hold on
    plot(timehist_EG2, gaphist_EG2, '-o', 'LineWidth', 3, 'MarkerSize', 10);
    hold on
    plot(timehist_OGDA2, gaphist_OGDA2, '-s', 'LineWidth', 3, 'MarkerSize', 10);
    hold off
    legend({'EG1', 'OGDA1', 'Our Method', 'EG2', 'OGDA2'}, 'Location', 'northwest', ...
        'Orientation', 'horizontal', 'NumColumns', 3);
    
    set(gca, 'YScale','log');
    set(gca, 'FontSize', 15);
    xlabel('Time (Seconds)');
    ylabel('Gap Function');
    xlim([0 x_time(di)])
    ylim([yaxis_lower(di) yaxis_upper(di)])
    title(['n=', num2str(n)]);
    
    path = sprintf('../figs/cubic_%d_time', n); 
    saveas(gcf, path, 'epsc');

end