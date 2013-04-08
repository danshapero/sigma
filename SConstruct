env = Environment(tools=['default','gfortran'],F90='gfortran',LINK='gfortran',F90PATH=['.', 'src/linalg'])

linalg_sources = Split("""src/linalg/sparse_matrix_mod.f90 
        src/linalg/csr_matrix_mod.f90 src/linalg/iterative_solver_mod.f90
        src/linalg/cg_solver_mod.f90 src/linalg/nullpc_mod.f90
        src/linalg/jacobi_mod.f90""")

env.SharedLibrary('linalg',linalg_sources)
