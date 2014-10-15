set(GRAPH_SRCS
graph/graph_interfaces.f90
graph/formats/coo_graphs.f90
graph/formats/cs_graphs.f90
graph/formats/ellpack_graphs.f90
graph/formats/ll_graphs.f90
graph/permutations.f90
graph/graphs.f90
)

set(MATRIX_SRCS
matrix/sparse_matrix_interfaces.f90
matrix/formats/default_sparse_matrix_kernels.f90
matrix/formats/default_matrices.f90
matrix/formats/cs_matrices.f90
matrix/formats/ellpack_matrices.f90
matrix/sparse_matrix_composites.f90
matrix/sparse_matrix_algebra.f90
matrix/sparse_matrix_factory.f90
matrix/sparse_matrices.f90
)

set(SOLVER_SRCS
solver/bicgstab_solvers.f90
solver/gcr_solvers.f90
solver/cg_solvers.f90
solver/jacobi_solvers.f90
solver/ldu_solvers.f90
)

set(LINEAR_OPERATOR_SRCS
linear_operator/linear_operator_interface.f90
linear_operator/linear_operator_sums.f90
linear_operator/linear_operator_products.f90
linear_operator/linear_operator_adjoints.f90
linear_operator/linear_operators.f90
)

#add_library( sigma SHARED types.f90 util.f90 vectors.f90 ${LINEAR_OPERATOR_SRCS} ${GRAPH_SRCS} ${MATRIX_SRCS} ${SOLVER_SRCS} eigensolver.f90 permutations.f90 wrapper.f90 sigma.f90 )
add_library( sigma SHARED types.f90 util.f90 vectors.f90 ${LINEAR_OPERATOR_SRCS} ${GRAPH_SRCS} ${MATRIX_SRCS} ${SOLVER_SRCS} eigensolver.f90 sigma.f90 )
target_link_libraries(sigma blas lapack)

