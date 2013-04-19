module linear_algebra_mod

    ! Sparse matrix format modules
    use sparse_matrix_mod
    use csr_matrix_mod

    ! Permutation module
    use permutation_mod

    ! Solver modules
    use iterative_solver_mod
    use cg_solver_mod

    ! Preconditioner modules
    use nullpc_mod
    use jacobi_mod
    use ilu_mod

    implicit none

contains


!--------------------------------------------------------------------------!
subroutine solver_setup(solver,pc,solver_name,tol_name,pc_name,pc_level)   !
!--------------------------------------------------------------------------!
    implicit none
    class(iterative_solver), intent(inout), allocatable, optional :: solver
    class(preconditioner), intent(inout), allocatable, optional :: pc
    character(len=32), intent(in), optional :: solver_name,tol_name, &
        & pc_name,pc_level
    ! local variables
    real(kind(1d0)) :: tolerance

    ! Pick which solver to use
    if (present(solver)) then
        select case(trim(solver_name))
            case("cg")
            allocate(cg_solver::solver)
        case default
            ! Need to change this when I actually get a new solver, e.g.
            ! make it default to GMRES or BiCG or something that will
            ! actually work for general systems.
            allocate(cg_solver::solver)
        end select

        ! Set the solver tolerance; default is 1E-6
        if (present(tol_name)) then
            read (tol_name,*) tolerance
        else
            tolerance = 1.d-6
        endif
    endif

    
end subroutine solver_setup


end module linear_algebra_mod
