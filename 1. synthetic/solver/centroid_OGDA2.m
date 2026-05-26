%%************************************************************************
%% Call OGDA2 to solve cubic regularized minimax problems
function [z, gaphist, counthist, timehist] = centroid_OGDA2(A, b, rho, R, options)

n = length(b);

%% initialization
z = ones(2*n, 1); 
v = z; 

alpha = 0.5; 
beta = 0.8; 
eta_up = (100/32)^(3/2)/(4*rho); 
counter_LS = 0; 

nIter = 1000;
if isfield(options, 'max_iters')
    nIter = options.max_iters;
end

display = 1;         % option of displaying
displayfreq = 1;     % gap of display
checkfreq = 1;       % frequency of check
savecounthist = 0;   % save count history
savegaphist = 0;     % save gap history
savetimehist = 0;    % save time history
 
if isfield(options, 'display'),         display = options.display;              end    
if isfield(options, 'displayfreq'),     displayfreq = options.displayfreq;      end    
if isfield(options, 'checkfreq'),       checkfreq = options.checkfreq;          end
if isfield(options, 'savecounthist'),   savecounthist = options.savecounthist;  end
if isfield(options, 'savegaphist'),     savegaphist = options.savegaphist;      end
if isfield(options, 'savetimehist'),    savetimehist = options.savetimehist;    end
 
if display == 1
    fprintf('\n-------------- OGDA2 ---------------\n');
    fprintf('iter |   gap  |  count  |  time\n');
end
 
r = gap(A, b, rho, R, z);
gaphist = [r]; 
timehist = [0];
counthist = [0]; 
[F, J] = derivative(A, b, rho, z);

tstart = clock;
%% main loop
for iter = 1:nIter

    % warmstart: the initial trial stepsize is the inadmissible one in the last equation. 
    sigma = eta_up; 
    eta = sigma; 
    [z_next, F_next, res_next, flag_stepsize] = subproblem_OGDA2(z, F, J, v, A, b, rho, eta, alpha);
    counter_LS = counter_LS+1; 
    
    % backtracking and advancing subroutines
    if ~flag_stepsize
        % backtracking subroutine: the trial stepsize is inadmissible
        while ~flag_stepsize
            eta_up = eta; % record the inadmissible stepsize
            eta = beta*eta^2/sigma; % backtracking
            [z_next, F_next, res_next, flag_stepsize] = subproblem_OGDA2(z, F, J, v, A, b, rho, eta, alpha);
            counter_LS = counter_LS+1; 
        end

        % record the admissible stepsize
        eta_lo = eta; 
        z_lo = z_next; 
        F_lo = F_next; 
        v_lo = eta*res_next; 
    else
        % advancing subroutine: the trial stepsize is admissible
        while flag_stepsize && norm(F_next) >= 1e-16
            eta_lo = eta; 
            z_lo = z_next; 
            F_lo = F_next; 
            v_lo = eta*res_next; 
            eta = eta^2/beta/sigma; % advance
            [z_next, F_next, res_next, flag_stepsize] = subproblem_OGDA2(z, F, J, v, A, b, rho, eta, alpha);
            counter_LS = counter_LS+1; 
        end

        % record the admissible stepsize
        eta_up = eta;
    end

    % bisection subroutine
    while eta_up > eta_lo/beta
        eta = sqrt(eta_lo*eta_up);
        [z_next, F_next, res_next, flag_stepsize] = subproblem_OGDA2(z, F, J, v, A, b, rho, eta, alpha);
        counter_LS = counter_LS+1; 
        if flag_stepsize
            eta_lo = eta; 
            z_lo = z_next; 
            F_lo = F_next; 
            v_lo = eta*res_next; 
        else
            eta_up = eta;
        end
    end

    % update 
    z = z_lo;
    F = F_lo; 
    v = v_lo; 

    [~, J] = derivative(A, b, rho, z);
    
    r  = gap(A, b, rho, R, z);
    res_gap = abs(r)+1e-13; 
     
    if savecounthist == 1, counthist = [counthist; counter_LS]; end
    if savegaphist == 1, gaphist = [gaphist; res_gap]; end
    if savetimehist == 1, timehist = [timehist; etime(clock, tstart)]; end
    if (display == 1) && ((mod(iter, displayfreq) == 0) && (mod(iter, checkfreq) == 0))
        fprintf('%5.0f|%0.3e|%5.0f|%3.2e\n', iter, res_gap, counter_LS, etime(clock, tstart));
    end

    if res_gap < 1e-8
        break; 
    end
end

end

function [z_next, F_next, res_next, flag_stepsize] = subproblem_OGDA2(z, F, J, v, A, b, rho, eta, alpha)

n = length(F); 

[U, T] = schur(full(J));
U = sparse(U); 
T = sparse(T); 
tmp = -U'*(eta*F+v);
tmp1 = eta*T+speye(n); 
tmp2 = tmp1\tmp; 
z_delta = U*tmp2;
z_next = z+z_delta; 

[F_next, ~] = derivative(A, b, rho, z_next);
res_next = F_next-(F+J*z_delta); 
flag_stepsize = (eta*norm(res_next)<=0.5*alpha*norm(z_delta)); 

end