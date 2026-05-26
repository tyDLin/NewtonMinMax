function R = centroid_opt(A, b, rho)

x = A\b; 
y = -(rho/2)*norm(x)*(A'\x); 
R = norm(y); 

end