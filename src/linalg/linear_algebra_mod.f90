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
    use jacobi_mod
    use bjacobi_mod
    use ilu_mod

    implicit none

contains


!--------------------------------------------------------------------------!
subroutine solver_setup(A,solver,pc,solver_name,pc_name,tolerance,level)   !
!--------------------------------------------------------------------------!
    implicit none
    class(sparse_matrix), intent(in) :: A
    class(iterative_solver), intent(inout), allocatable :: solver
    class(preconditioner), intent(inout), allocatable :: pc
    character(len=32), intent(in), optional :: solver_name,pc_name, &
        & tolerance,level
    ! local variables
    integer :: nn,lev
    real(kind(1d0)) :: tol

    ! Pick which solver to use.
    if (present(solver_name)) then
        select case(trim(solver_name))
            case("cg")
                allocate(cg_solver::solver)
            case default
                ! Need to change this when I actually get a new solver, e.g.
                ! make it default to GMRES or BiCG or something that will
                ! actually work for general systems.
                allocate(cg_solver::solver)
        end select
    else
        allocate(cg_solver::solver)
    endif

    ! Set the solver tolerance; default is 1E-6
    if (present(tolerance)) then
        read (tolerance,*) tol
    else
        tol = 1.d-8
    endif

    call solver%init(A%nrow,tol)

    ! Pick which preconditioner to use.
    if (present(pc_name)) then
        select case(trim(pc_name))
            case("jacobi")
                allocate(jacobi::pc)
                lev = 0
            case("bjacobi")
                allocate(bjacobi::pc)
                lev = 16
            case("ilu")
                allocate(ilu::pc)
                lev = 0
            case default
                allocate(nopc::pc)
                lev = 0
        end select
    else
        allocate(nopc::pc)
    endif

    if (present(level)) then
        read(level,*) lev
    endif

    call pc%init(A,lev)
    
end subroutine solver_setup


end module linear_algebra_mod
