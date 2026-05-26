function res = cubic_gap(A, b, rho, R, z)  

n = length(b); 
x = z(1:n); 
y = z(n+1:2*n); 
res = (rho/6)*norm(x)*norm(x)*norm(x) + R*norm(A*x-b) + (2/3)*sqrt(2/rho)*(norm(A'*y))^(3/2) + b'*y; 

end
