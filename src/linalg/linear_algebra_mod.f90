module linear_algebra_mod

    ! Sparse matrix format modules
    use csr_matrix_mod
    use sparse_matrix_mod

    ! Permutation module
    use permutation_mod

    ! Solver modules
    use iterative_solver_mod
    use cg_solver_mod

    ! Preconditioner modules
    use nullpc_mod
    use jacobi_mod

    implicit none

contains

! This module contains nothing. It is a wrapper which uses all of the other
! modules, so that we need only include a single module to get all of the
! linear algebra modules in an application.

end module linear_algebra_mod
