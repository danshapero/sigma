This directory contains sources for linear solvers.
So far, we have implemented:

* conjugate gradient method (`cg_solvers.f90`): iterative Krylov subspace method (KSM) for symmetric positive definite matrices
* biconjugate gradient method (`bicg_solvers.f90`): KSM for general asymmetric matrices
* Jacobi's method (`jacobi_solvers.f90`): diagonal preconditioner, best for diagonally-dominant matrices
* incomplete LDU (`ldu_solvers.f90`): preconditioner for general matrices, reduces to incomplete Cholesky for symmetric definite matrices

Some solvers operate at the level of linear operators (CG, BiCG) while others (LDU) only work on a sparse matrix.

There is no distinction in the class hierarchy between solvers and preconditioners.
With this usage, a solver should be thought of as some procedure that, given an initial guess `x0` for the solution of `A*x = b`, can, for some matrices, give a better guess `x1`.
In this sense, even the Jacobi method can be thought of as a "solver", it just happens to be a very inefficient one.
This choice is intended to facilitate the implementation of multilevel methods such as multigrid or domain-decomposition, which can themselves be used as either solvers or preconditioners of outer-level KSMs.
Domain-decomposition methods can have rather deep nested solvers, i.e. the local sub-domain solves can be performed with ILDU-preconditioned Krylov subspace methods or even multigrid.
For such scenarios, enforcing a distinction between solvers and preconditioners is counter-productive.
