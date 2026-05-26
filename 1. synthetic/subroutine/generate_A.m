function A = generate_A(n)  

A = eye(n) + diag(-ones(n-1,1),1); 

A = sparse(A); 
end
