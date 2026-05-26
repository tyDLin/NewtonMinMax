%%************************************************************************
%% Call SSNDE to solve AUC maximization problems
function [z, gaphist, samplehist, timehist] = centroid_SSNDE(z_star, A, b, rho, eta, options)

[n, d] = size(A);

%% initialization
z = zeros(d+3, 1);
z_hat = z; 

nIter = 100;
if isfield(options, 'max_iters')
    nIter = options.max_iters;
end

display = 1;         % option of displaying
displayfreq = 1;     % gap of display
checkfreq = 1;       % frequency of check
savegaphist = 0;     % save gap history
savesamplehist = 0;  % save epoch history
savetimehist = 0;    % save time history  
 
if isfield(options, 'display'),       display = options.display;            end    
if isfield(options, 'displayfreq'),   displayfreq = options.displayfreq;    end   
if isfield(options, 'checkfreq'),     checkfreq = options.checkfreq;        end  
if isfield(options, 'savegaphist'),   savegaphist = options.savegaphist;    end
if isfield(options, 'savesamplehist'),  savesamplehist = options.savesamplehist;  end
if isfield(options, 'savetimehist'),  savetimehist = options.savetimehist;  end
 
if display == 1
    fprintf('\n-------------- SSNDE ---------------\n');
    fprintf('iter |   gap  |  epoch  |  time\n');
end
 
r = auc_gap(z_star, A, b, rho, z); 
gaphist = [r];
samplehist = [0]; 
timehist = [0];
num_samples = 0; 

tstart = clock;
g = auc_grad(z, A, b, rho, n);
%% main loop
for iter = 1:nIter

    % compute the gradient and Hessian at z_hat
    g_hat = auc_grad(z_hat, A, b, rho, n);
    g_min = min(norm(g), norm(g_hat));
    if g_min > 1e-4
        M = ceil(min(n, max(eta*n, 8*log(d+3)/(g_min)^2)));
    else
        M = ceil(min(n, max(eta*n, 25*log(d+3)/(g_min)^2)));
    end
    J_hat = auc_jacb(z_hat, A, b, rho, M);  

    % subproblem solving
    z_delta = subproblem_NDE(g_hat, J_hat, rho);

    % compute the next iterate
    z = z_hat + z_delta;

    % compute the stepsize
    if M < 0.8*n
        lambda = 1/(1000*rho*norm(z_delta));
    else
        lambda = 1/(10*rho*norm(z_delta)); 
    end

    % compute the gradient at z
    g = auc_grad(z, A, b, rho, n);

    % update 
    z_hat = z_hat - lambda*g; 
    
    num_samples = num_samples+M; 

    r = auc_gap(z_star, A, b, rho, z);
    res_gap = abs(r)+1e-13; 
    
    if savegaphist == 1, gaphist = [gaphist; res_gap]; end
    if savesamplehist == 1; samplehist = [samplehist; num_samples/n]; end
    if savetimehist == 1, timehist = [timehist; etime(clock, tstart)]; end
    if (display == 1) && ((mod(iter, displayfreq) == 0) && (mod(iter, checkfreq) == 0))
        fprintf('%5.0f|%0.3e|%3.2e|%3.2e\n', iter, res_gap, num_samples/n, etime(clock, tstart));
    end

    if res_gap < 1e-6
        break; 
    end
end

end

function z_delta = subproblem_NDE(F, J, rho)

n = length(F);

[U, T] = schur(full(J));
U = sparse(U); 
T = sparse(T); 
tmp = -U'*F;

lambda = 1e-4;
tmp1 = T+lambda*speye(n); 
tmp2 = tmp1\tmp; 
z_delta = U*tmp2;
f = norm(z_delta) - lambda/(10*rho); 

while abs(f)> 1e-8
    tmp3 = U'*z_delta; 
    g = -tmp3'*(tmp1\tmp3)/norm(z_delta)-1/(10*rho); 
    lambda = lambda - f/g;

    tmp1 = T+lambda*speye(n); 
    tmp2 = tmp1\tmp; 
    z_delta = U*tmp2;
    f = norm(z_delta) - lambda/(10*rho); 
end

tmp1 = T+lambda*speye(n); 
tmp2 = tmp1\tmp; 
z_delta = U*tmp2;
end