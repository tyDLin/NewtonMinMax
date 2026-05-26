%%************************************************************************
%% Call EG2 to solve cubic regularized minimax problems
function [z, gaphist, counthist, timehist] = centroid_EG2(A, b, rho, R, options)

n = length(b);

%% initialization
z = ones(2*n, 1);

nIter = 500;
if isfield(options, 'max_iters')
    nIter = options.max_iters;
end

counter_LS = 0; 

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
    fprintf('\n-------------- EG2 ---------------\n');
    fprintf('iter |   gap  |  count  |  time\n');
end
 
r = gap(A, b, rho, R, z);
gaphist = [r]; 
timehist = [0];
counthist = [0]; 

[F, J] = derivative(A, b, rho, z);
J_small = svds(J, 1, 'smallest'); 
F_norm = norm(F); 

tstart = clock;
%% main loop
for iter = 1:nIter

    % warmstart: the initial trial stepsize is the inadmissible one in the last equation. 
    eta_lo = J_small/(rho*F_norm); 
    eta_up = (100/32)^(3/2)/(4*rho); 
    [z_half, F_half] = subproblem_EG2(z, F, J, A, b, rho, eta_up);
    counter_LS = counter_LS+1; 
    
    % bisection subroutine
    if eta_up < 1/(8*rho*norm(z_half-z))
        % the upper trial stepsize is admissible
        eta = eta_up;
        z_next = z - eta*F_half; 
    elseif eta_lo >= eta_up
        % the lower trial stepsize is admissible
        eta = eta_lo; 
        [~, F_half] = subproblem_EG2(z, F, J, A, b, rho, eta);
        counter_LS = counter_LS+1; 
        z_next = z - eta*F_half; 
    else
        % bisection subroutine: the trial stepsize is inadmissible
        eta = 0.5*(eta_lo+eta_up);
        while abs(eta_up-eta_lo)/eta_up > 0.1
            [z_half, F_half] = subproblem_EG2(z, F, J, A, b, rho, eta);
            counter_LS = counter_LS+1; 
            D = 1/(12*rho*norm(z_half-z)); 
            if eta <= D
                eta_lo = eta; 
            else
                eta_up = eta; 
            end
            eta = 0.5*(eta_lo+eta_up);
        end
        z_next = z - eta*F_half; 
    end

    % update 
    z = z_next;
    
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
    
    [F, J] = derivative(A, b, rho, z);
    J_small = svds(J, 1, 'smallest'); 
    F_norm = norm(F); 
end

end

function [z_half, F_half] = subproblem_EG2(z, F, J, A, b, rho, eta)

n = length(F); 

[U, T] = schur(full(J));
U = sparse(U); 
T = sparse(T); 
tmp = -U'*(eta*F);
tmp1 = eta*T+speye(n); 
tmp2 = tmp1\tmp; 
z_delta = U*tmp2; 
z_half = z+z_delta; 

[F_half, ~] = derivative(A, b, rho, z_half);

end