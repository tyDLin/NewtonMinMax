%%************************************************************************
%% Call Newton dual extrapolation to solve cubic regularized minimax problems
function [z, gaphist, timehist] = centroid_NDE(A, b, rho, R, options)

n = length(b);

%% initialization
z = ones(2*n, 1);
z_hat = z; 

nIter = 100;
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
    fprintf('\n-------------- NDE ---------------\n');
    fprintf('iter |   gap  |   time\n');
end
 
r = gap(A, b, rho, R, z);
gaphist = [r]; 
timehist = [0];

tstart = clock;
%% main loop
for iter = 1:nIter

    % compute the gradient and Hessian at z_hat
    [F_hat, J_hat] = derivative(A, b, rho, z_hat);

    % subproblem solving
    z_delta = subproblem_NDE(F_hat, J_hat, rho);

    % compute the next iterate
    z = z_hat + z_delta;

    % compute the stepsize
    lambda = 1/(10*rho*norm(z_delta)); 

    % compute the gradient at z
    [F, ~] = derivative(A, b, rho, z);

    % update 
    z_hat = z_hat - lambda*F; 
    
    r  = gap(A, b, rho, R, z);
    res_gap = abs(r)+1e-13; 
    
    if savegaphist == 1, gaphist = [gaphist; res_gap]; end
    if savetimehist == 1, timehist = [timehist; etime(clock, tstart)]; end
    if (display == 1) && ((mod(iter, displayfreq) == 0) && (mod(iter, checkfreq) == 0))
        fprintf('%5.0f|%0.3e|%3.2e\n', iter, res_gap, etime(clock, tstart));
    end

    if res_gap < 1e-8
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