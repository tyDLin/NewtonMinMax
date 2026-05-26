%%************************************************************************
%% Call SEG to solve AUC maximization problems
function [z, gaphist, iterhist, timehist] = centroid_SVREG(z_star, A, b, rho, options)

[n, d] = size(A);

%% initialization
z = zeros(d+3, 1);
w = zeros(d+3, 1);
p = 0.0001; 
L = 0.2; 

nIter = 1000;
if isfield(options, 'max_iters')
    nIter = options.max_iters;
end

display = 1;         % option of displaying
displayfreq = 1;     % gap of display
checkfreq = 1;       % frequency of check
savegaphist = 0;     % save gap history
saveiterhist = 0;    % save epoch history
savetimehist = 0;    % save time history    
   
if isfield(options, 'display'),       display = options.display;            end    
if isfield(options, 'displayfreq'),   displayfreq = options.displayfreq;    end 
if isfield(options, 'checkfreq'),     checkfreq = options.checkfreq;        end  
if isfield(options, 'savegaphist'),   savegaphist = options.savegaphist;    end
if isfield(options, 'saveiterhist'),  saveiterhist = options.saveiterhist;  end
if isfield(options, 'savetimehist'),  savetimehist = options.savetimehist;  end

if display == 1
    fprintf('\n-------------- SVREG ---------------\n');
    fprintf('iter |   gap  |  epoch  |  time\n');
end

gap = auc_gap(z_star, A, b, rho, z); 
gaphist = [gap];
iterhist = [0]; 
timehist = [0];
num_samples = 0;
M = ceil(0.1*n); 

tstart = clock;
g_w = auc_grad(w, A, b, rho, n);
z_old = z; 
w_old = w; 

%% main loop
for iter = 1:nIter

    % set the intermediate iterate
    z_bar = (1-p)*z_old + p*w_old; 

    % compute the gradient at w
    if w ~= w_old; g_w = auc_grad(w, A, b, rho, n); end
 
    % compute the half iterate
    z_half = z_bar - L*g_w; 
    
    % compute the mini-batch gradient at z_half and w
    sg_half = auc_grad(z_half, A, b, rho, M);
    sg_w = auc_grad(w, A, b, rho, M);
    
    % compute the next iterate
    z = z_bar - L*(g_w + sg_half - sg_w); 
    z_old = z; 
    
    % compute the total number of samples processed
    num_samples = num_samples + M; 
    if w ~= w_old; num_samples = num_samples + n; end

    % set the auxillary iterate
    w_old = w; 
    r = rand(1); 
    if r < p; w = z; end

    if iter == 1 || mod(iter, checkfreq) == 0
        gap = auc_gap(z_star, A, b, rho, z);  
    end
    
    if savegaphist == 1; gaphist = [gaphist; gap]; end
    if saveiterhist == 1; iterhist = [iterhist; num_samples/n]; end
    if savetimehist == 1; timehist = [timehist; etime(clock, tstart)]; end
    if (display == 1) && ((mod(iter, displayfreq) == 0) && (mod(iter, checkfreq) == 0))
        fprintf('%5.0f|%0.3e|%3.2e|%3.2e\n', iter, gap, num_samples/n, etime(clock, tstart));
    end
end

end