set(GRAPH_SRCS
graph/graphs.f90
graph/coo_graphs.f90
graph/cs_graphs.f90
graph/ellpack_graphs.f90
graph/ll_graphs.f90
)

set(MATRIX_SRCS
matrix/sparse_matrices.f90
)

set(SOLVER_SRCS
solver/iterative_solvers.f90
solver/bicgstab_solvers.f90
solver/cg_solvers.f90
solver/jacobi_preconditioners.f90
)

add_library( sigma SHARED types.f90 util.f90 vectors.f90 ${GRAPH_SRCS} ${MATRIX_SRCS} ${SOLVER_SRCS} permutations.f90 conversions.f90 meshes.f90 wrapper.f90 sigma.f90 )
target_link_libraries(sigma blas lapack)
