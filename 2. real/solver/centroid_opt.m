%%************************************************************************
%% Call NLS to solve AUC maximization
%%************************************************************************
function z = centroid_opt(A, b, rho)

[n, d] = size(A);

%% problem parameter
p = sum(b > 0)/n;

A1 = A'*(b>0); 
A1 = -(2/n)*(1-p)*A1; 
A2 = A'*(b<0);
A2 = -(2/n)*p*A2; 

q = [A1-A2; zeros(2,1)]; 

H = auc_jacb(zeros(d+3), A, b, rho, n);
P = H(1:d+2,1:d+2); 
a = H(1:d+2,d+3); 
p = H(d+3,d+3); 

%% compute the optimal solution
fun = @(x)subproblem_opt(x, P, q, a, p, rho);
options = optimset('Display', 'off');
x = fsolve(fun, zeros(d+2, 1), options);
y = a'*x/p;

z = [x; y]; 

end

function f = subproblem_opt(x, P, q, a, p, rho)

f = q+a*(a'*x)/p+0.5*rho*norm(x)*x+P*x; 

end