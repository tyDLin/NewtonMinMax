# Explicit Second-Order Min-Max Optimization

These codes provide implementations of solvers for solving cubic regularized bilinear min-max problems and AUC maximization problems using inexact regularized Newton-type methods. 

# About

We propose and analyze several inexact regularized Newton-type methods for finding a global saddle point of convex-concave unconstrained min-max optimization problems. Compared to first-order methods, our understanding of second-order methods for min-max optimization is relatively limited, as obtaining global rates of convergence with second-order information can be much more involved. 

In this paper, we examine how second-order information is used to speed up extra-gradient methods, even under inexactness. In particular, we show that the proposed methods generate iterates that remain within a bounded set and that the averaged iterates converge to an $\epsilon$-saddle point within $O(\epsilon^{-2/3})$ iterations in terms of a restricted gap function. We also provide a simple routine for solving the subproblem at each iteration, requiring a single Schur decomposition and $O(\log\log(1/\epsilon))$ calls to a linear system solver in a quasi-upper-triangular system. Thus, our method improves the existing line-search-based second-order min-max optimization methods by shaving off an $O(\log\log(1/\epsilon))$ factor in the required number of Schur decompositions. Finally, we evaluate our method on both synthetic benchmarks and a real-world application arising from AUC maximization on standard LIBSVM datasets, and find that the proposed second-order approach delivers stronger practical efficiency than representative first-order methods on these problems.

# Codes

The MATLAB Implementations on Synthetic and LIBSVM Data are provided.  

# References

T. Lin, P. Mertikopoulos and M. I. Jordan. Explicit Second-Order Min-Max Optimization: Practical Algorithms and Complexity Analysis. Transactions on Machine Learning Research (TMLR), 2026.
