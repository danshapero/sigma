module linalg

    ! Sparse matrix format modules
    use bsr
    use ellpack
    use csr
    use matrix

    ! Permutation module
    use permutation

    ! Solver modules
    use solver
    use cg
    use cgs

    ! Preconditioner modules
    use jacobi
    use bjacobi
    use ilu

    implicit none

contains



!--------------------------------------------------------------------------!
subroutine solver_setup(A,solver,pc,solver_name,pc_name,tolerance,level)   !
!--------------------------------------------------------------------------!
    implicit none
    class(sparse_matrix), intent(in) :: A
    class(iterative_solver), intent(inout), allocatable :: solver
    class(preconditioner), intent(inout), allocatable :: pc
    character(len=*), intent(in), optional :: solver_name,pc_name, &
        & tolerance,level
    ! local variables
    integer :: nn,lev
    real(kind(1d0)) :: tol

    ! Pick which solver to use.
    if (present(solver_name)) then
        select case(trim(solver_name))
            case("cg")
                allocate(cg_solver::solver)
            case("cgs")
                allocate(cgs_solver::solver)
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
                allocate(jacobi_preconditioner::pc)
                lev = 0
            case("bjacobi")
                allocate(bjacobi_preconditioner::pc)
                lev = 16
            case("ilu")
                allocate(ilu_preconditioner::pc)
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



end module linalg
