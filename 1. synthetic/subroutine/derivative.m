function [grad, hess] = derivative(A, b, rho, z)  

%% return F and DF of cubic regularized blinear function: f(x,y)=(rho/6)x^3 + y^T(Ax-b)
n = length(b);
x = z(1:n); 
y = z(n+1:2*n); 

grad = [0.5*rho*norm(x)*x+A'*y; b-A*x]; 
if norm(x) == 0
    hess = [zeros(n) A'; -A zeros(n)]; 
else
    hess = [0.5*rho*norm(x)*eye(n)+0.5*rho*(x*x')/norm(x) A'; -A zeros(n)]; 
end

end
