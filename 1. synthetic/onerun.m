%%
clear;
clc; 
close all; 

addpath('solver')
addpath('subroutine')

ranseed = 1;
rng(ranseed, 'twister');

%% Problem setting: generate A
ndim        = [500];          
T           = 50000; 

for di=1:length(ndim)
        
    %% generate data. 
    n = ndim(di);
    A = generate_A(n); 
    rho = 1/(20*n);
    
    %% generate data
    b = -1 + 2*rand(n, 1);
    R = centroid_opt(A, b, rho); 

    %% call EG1
    optsEG1.max_iters       = T;
    optsEG1.display         = 0;
    optsEG1.displayfreq     = 100;
    optsEG1.checkfreq       = 100;
    optsEG1.savegaphist     = 1;
    optsEG1.savetimehist    = 1;
    [~, gaphist_EG1, timehist_EG1] = centroid_EG1(A, b, rho, R, optsEG1);

    %% call OGDA1
    optsOGDA1.max_iters     = 2*T;
    optsOGDA1.display       = 0;
    optsOGDA1.displayfreq   = 100;
    optsOGDA1.checkfreq     = 100;
    optsOGDA1.savegaphist   = 1;
    optsOGDA1.savetimehist  = 1;
    [~, gaphist_OGDA1, timehist_OGDA1] = centroid_OGDA1(A, b, rho, R, optsOGDA1);

    %% call EG2
    optsEG2.max_iters       = 100;
    optsEG2.display         = 0;
    optsEG2.displayfreq     = 1;
    optsEG2.checkfreq       = 1;
    optsEG2.savecounthist   = 1; 
    optsEG2.savegaphist     = 1;
    optsEG2.savetimehist    = 1;
    [~, gaphist_EG2, counthist_EG2, timehist_EG2] = centroid_EG2(A, b, rho, R, optsEG2);

    %% call OGDA2
    optsOGDA2.max_iters     = 100;
    optsOGDA2.display       = 0;
    optsOGDA2.displayfreq   = 1;
    optsOGDA2.checkfreq     = 1;
    optsOGDA2.savecounthist = 1; 
    optsOGDA2.savegaphist   = 1;
    optsOGDA2.savetimehist  = 1;
    [~, gaphist_OGDA2, counthist_OGDA2, timehist_OGDA2] = centroid_OGDA2(A, b, rho, R, optsOGDA2);

    %% call NDE
    optsNDE.max_iters       = 100;
    optsNDE.display         = 0;
    optsNDE.displayfreq     = 1;
    optsNDE.checkfreq       = 1;
    optsNDE.savegaphist     = 1;
    optsNDE.savetimehist    = 1;
    [~, gaphist_NDE, timehist_NDE] = centroid_NDE(A, b, rho, R, optsNDE);

    filename = sprintf('../results/cubic_%d.mat', n);
    disp(filename);
    save(filename, 'gaphist_EG1', 'gaphist_OGDA1', 'gaphist_NDE', 'counthist_EG2', 'gaphist_EG2', 'counthist_OGDA2', 'gaphist_OGDA2', ...
        'timehist_EG1', 'timehist_OGDA1', 'timehist_NDE', 'timehist_EG2', 'timehist_OGDA2');
end