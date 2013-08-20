module fempack

! Use all the graph modules
use graphs
use ll_graphs
use coo_graphs
use cs_graphs

! Use all the matrix modules
use sparse_matrices
use csr_matrices
use coo_matrices

! Use all the block matrix modules
use block_sparse_matrices
use bcsr_matrices

! Use the solver and preconditioner modules
use iterative_solvers
use cg_solvers
use jacobi_preconditioners

! Use the C wrapper module
use wrapper

! Use other auxiliary modules
use conversions
use meshes


implicit none




end module fempack
