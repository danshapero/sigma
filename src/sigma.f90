module sigma

! Use all the graph modules
use graph_interfaces
use ll_graphs
use coo_graphs
use cs_graphs
use ellpack_graphs
use graphs

! Use all the linear operator modules
use linear_operator_interface
use linear_operator_sums
use linear_operator_products
use linear_operators

! Use all the matrix modules
use sparse_matrix_interfaces
use default_sparse_matrix_kernels
use default_matrices
use cs_matrices
use ellpack_matrices
use sparse_matrix_factory
use sparse_matrix_algebra
use sparse_matrices

! Use the solver and preconditioner modules
use cg_solvers
use bicgstab_solvers
use jacobi_solvers
use ldu_solvers

! Use the C wrapper module
!use wrapper

! Use other auxiliary modules
use eigensolver
use permutations
use vectors




end module sigma
