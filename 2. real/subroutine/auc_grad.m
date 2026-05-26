function g = auc_grad(z, A, b, rho, M)  

%% return the gradient of objective function in AUC maximization
[n, d] = size(A);

%% initialization
p = sum(b > 0)/n;
g = zeros(d+3, 1); 

theta = z(1:d); 
u = z(d+1);
v = z(d+2); 
y = z(d+3); 

%% full gradient
if M >= n
    tmp = A*theta;
    a1 = (tmp-u).*(b>0); 
    a1 = (2/n)*(1-p)*a1; 
    a2 = (tmp-v).*(b<0); 
    a2 = (2/n)*p*a2; 
    a3 = (2/n)*(p*ones(n,1)-(b>0)); 
end

%% subsampled gradient
if M < n
    I = randperm(n);
    A = A(I(1:M),:);
    b = b(I(1:M));
    tmp = A*theta;
    a1 = (tmp-u).*(b>0); 
    a1 = (2/M)*(1-p)*a1; 
    a2 = (tmp-v).*(b<0); 
    a2 = (2/M)*p*a2; 
    a3 = (2/M)*(p*ones(M,1)-(b>0));
end

g(1:d) = A'*(a1+a2+(1+y)*a3)+0.5*rho*norm(z(1:d+2))*theta; 
g(d+1) = -sum(a1)+0.5*rho*norm(z(1:d+2))*u; 
g(d+2) = -sum(a2)+0.5*rho*norm(z(1:d+2))*v;
g(d+3) = -theta'*A'*a3+2*y*p*(1-p); 

end
