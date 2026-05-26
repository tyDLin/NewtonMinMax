%%************************************************************************
%% Call OGDA1 to solve cubic regularized minimax problems
function [z, gaphist, timehist] = centroid_OGDA1(A, b, rho, R, options)

n = length(b);

%% initialization
z = ones(2*n, 1); 
z_old = z; 
L = 0.05; 

nIter = 20000;
if isfield(options, 'max_iters')
    nIter = options.max_iters;
end

display = 1;            % option of displaying
displayfreq = 1;        % gap of display
checkfreq = 1;          % frequency of check
savegaphist = 0;        % save gap history
savetimehist = 0;       % save time history
 
if isfield(options, 'display'),       display = options.display;            end    
if isfield(options, 'displayfreq'),   displayfreq = options.displayfreq;    end    
if isfield(options, 'checkfreq'),     checkfreq = options.checkfreq;        end   
if isfield(options, 'savegaphist'),   savegaphist = options.savegaphist;    end
if isfield(options, 'savetimehist'),  savetimehist = options.savetimehist;  end
 
if display == 1
    fprintf('\n-------------- OGDA1 ---------------\n');
    fprintf('iter |   gap  |   time\n');
end
 
r = gap(A, b, rho, R, z);
gaphist = [r]; 
timehist = [0];
[F_old, ~] = derivative(A, b, rho, z_old);

tstart = clock;
%% main loop
for iter = 1:nIter

    % compute F(z)
    [F, ~] = derivative(A, b, rho, z);

    % compute the next iterate
    z_next = z - L*(2*F - F_old); 

    % update 
    z = z_next;
    F_old = F; 
    
    if iter == 1 || mod(iter, checkfreq) == 0
        r  = gap(A, b, rho, R, z);
    end
     
    if savegaphist == 1, gaphist = [gaphist; r]; end
    if savetimehist == 1, timehist = [timehist; etime(clock, tstart)]; end
    if (display == 1) && ((mod(iter, displayfreq) == 0) && (mod(iter, checkfreq) == 0))
        fprintf('%5.0f|%0.3e|%3.2e\n', iter, r, etime(clock, tstart));
    end
end

end