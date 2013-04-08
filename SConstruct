env = Environment(tools=['default','gfortran'],F90='gfortran',LINK='gfortran',F90PATH=['.', 'src/linalg'])

env.SharedLibrary('linalg',['src/linalg/sparse_matrix_mod.f90','src/linalg/csr_matrix_mod.f90'])
