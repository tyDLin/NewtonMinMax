function H = auc_jacb(z, A, b, rho, M)  

%% return the gradient of objective function in AUC maximization
[n, d] = size(A);

%% initialization
p = sum(b > 0)/n;
H = sparse(d+3,d+3); 

%% full Hessian
if M >= n
    tmp = (2/n)*(-b*p+(b>0));
    A1 = A'*(b>0); 
    A1 = -(2/n)*(1-p)*A1; 
    A2 = A'*(b<0);
    A2 = -(2/n)*p*A2; 
    M = n;
end

%% subsampled Hessian
if M < n
    I = randperm(n);
    A = A(I(1:M),:);
    b = b(I(1:M)); 
    tmp = (2/M)*(-b*p+(b>0));
    A1 = A'*(b>0); 
    A1 = -(2/M)*(1-p)*A1; 
    A2 = A'*(b<0);
    A2 = -(2/M)*p*A2; 
end

H(1:d,1:d) = A'*spdiags(tmp,0,M,M)*A;
H(1:d,d+1) = A1;
H(1:d,d+2) = A2; 
H(1:d,d+3) = A1-A2; 
H(d+1,1:d) = A1'; 
H(d+2,1:d) = A2'; 
H(d+3,1:d) = -A1'+A2'; 
H(d+1,d+1) = 2*p*(1-p); 
H(d+2,d+2) = 2*p*(1-p); 
H(d+3,d+3) = 2*p*(1-p);

% cubic regularization term
tmp = z(1:d+2);
if norm(tmp) > 0
    H(1:d+2,1:d+2) = H(1:d+2,1:d+2)+0.5*rho*norm(tmp)*eye(d+2)+0.5*rho*(tmp*tmp')/norm(tmp); 
end

end
