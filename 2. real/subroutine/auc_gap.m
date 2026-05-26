function res = auc_gap(z_star, A, b, rho, z)  

%% return the gap function value of objective function in AUC maximization
[n, d] = size(A);

%% problem parameter
p = sum(b > 0)/n;

A1 = A'*(b>0); 
A1 = -(2/n)*(1-p)*A1; 
A2 = A'*(b<0);
A2 = -(2/n)*p*A2; 

q = [A1-A2; zeros(3,1)]; 

H = auc_jacb(zeros(d+3), A, b, rho, n);
H(d+3,:) = - H(d+3,:); 

a_min = [z(1:d+2,1)-z_star(1:d+2,1); z_star(d+3)-z(d+3)];   
a_max = z+z_star; 

res = (1/2)*a_min'*H*a_max+q'*a_min+(rho/6)*norm(z(1:d+2,1))*norm(z(1:d+2,1))*norm(z(1:d+2,1)) ...
    -(rho/6)*norm(z_star(1:d+2,1))*norm(z_star(1:d+2,1))*norm(z_star(1:d+2,1));

end
