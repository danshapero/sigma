set(GRAPH_SRCS
graph/graphs.f90
graph/coo_graphs.f90
graph/cs_graphs.f90
graph/ellpack_graphs.f90
graph/ll_graphs.f90
)

set(MATRIX_SRCS
matrix/sparse_matrices.f90
matrix/block_sparse_matrices.f90
matrix/coo_matrices.f90
matrix/cs_matrices.f90
matrix/ellpack_matrices.f90
matrix/ll_matrices.f90
matrix/bcs_matrices.f90
)

set(SOLVER_SRCS
solver/iterative_solvers.f90
solver/bicgstab_solvers.f90
solver/cg_solvers.f90
solver/jacobi_preconditioners.f90
)

add_library( fempack SHARED types.f90 util.f90 vectors.f90 ${GRAPH_SRCS} ${MATRIX_SRCS} ${SOLVER_SRCS} permutations.f90 conversions.f90 meshes.f90 wrapper.f90 fempack.f90 )
target_link_libraries(fempack blas lapack)
