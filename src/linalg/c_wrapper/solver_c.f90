module solver_c

    use iso_c_binding
    use matrix_c
    use cg
    use cgs
    use jacobi
    use bjacobi
    use ilu

    implicit none


!--------------------------------------------------------------------------!
! C wrappers to Fortran solver and preconditioner data types               !
!--------------------------------------------------------------------------!
type, bind(c) :: iterative_solver_c
    type(c_ptr) :: p
    integer(c_int) :: nn,solver_type
end type iterative_solver_c


type, bind(c) :: preconditioner_c
    type(c_ptr) :: p
    integer(c_int) :: nn,level,pc_type
end type preconditioner_c



contains





!==========================================================================!
!==========================================================================!
!==== Subroutines to get Fortran pointers from C pointers to solvers   ====!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine get_iterative_solver_c(csolver,solver_type) bind(c)             !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type(iterative_solver_c), intent(inout) :: csolver
    integer(c_int), intent(in), value :: solver_type
    ! local variables
    type(cg_solver), pointer :: cg_sol
    type(cgs_solver), pointer :: cgs_sol

    csolver%solver_type = solver_type
    select case(solver_type)
        case(0)
            allocate(cg_sol)
            csolver%p = c_loc(cg_sol)
        case(1)
            allocate(cgs_sol)
            csolver%p = c_loc(cgs_sol)
    end select

end subroutine get_iterative_solver_c



!--------------------------------------------------------------------------!
subroutine iterative_solver_f(fsolver,csolver)                             !
!--------------------------------------------------------------------------!
    implicit none
    class(iterative_solver), pointer, intent(out) :: fsolver
    type(iterative_solver_c), intent(in) :: csolver
    ! local variables
    type(cg_solver), pointer :: cg_sol
    type(cgs_solver), pointer :: cgs_sol

    select case(csolver%solver_type)
        case(0)
            allocate(cg_sol)
            call c_f_pointer(csolver%p,cg_sol)
            fsolver => cg_sol
        case(1)
            allocate(cgs_sol)
            call c_f_pointer(csolver%p,cgs_sol)
            fsolver => cgs_sol
    end select

end subroutine iterative_solver_f



!--------------------------------------------------------------------------!
subroutine get_preconditioner_c(cpc,pc_type) bind(c)                       !
!--------------------------------------------------------------------------!
    implicit none
    type(preconditioner_c), intent(inout) :: cpc
    integer(c_int), intent(in), value :: pc_type
    ! local variables
    type(jacobi_preconditioner), pointer :: jacobi_pc
    type(bjacobi_preconditioner), pointer :: bjacobi_pc
    type(ilu_preconditioner), pointer :: ilu_pc

    cpc%pc_type = pc_type
    select case(pc_type)
        case(0)
            allocate(jacobi_pc)
            cpc%p = c_loc(jacobi_pc)
        case(1)
            allocate(bjacobi_pc)
            cpc%p = c_loc(bjacobi_pc)
        case(2)
            allocate(ilu_pc)
            cpc%p = c_loc(ilu_pc)
    end select

end subroutine get_preconditioner_c



!--------------------------------------------------------------------------!
subroutine preconditioner_f(fpc,cpc)                                       !
!--------------------------------------------------------------------------!
    implicit none
    class(preconditioner), pointer, intent(out) :: fpc
    type(preconditioner_c), intent(in) :: cpc
    ! local variables
    type(jacobi_preconditioner), pointer :: jacobi_pc
    type(bjacobi_preconditioner), pointer :: bjacobi_pc
    type(ilu_preconditioner), pointer :: ilu_pc

    select case(cpc%pc_type)
        case(0)
            allocate(jacobi_pc)
            call c_f_pointer(cpc%p,jacobi_pc)
            fpc => jacobi_pc
        case(1)
            allocate(bjacobi_pc)
            call c_f_pointer(cpc%p,bjacobi_pc)
            fpc => bjacobi_pc
        case(2)
            allocate(ilu_pc)
            call c_f_pointer(cpc%p,ilu_pc)
            fpc => ilu_pc
    end select

end subroutine preconditioner_f





!==========================================================================!
!==========================================================================!
!==== C wrapper subroutines to iterative solver type-bound procedures  ====!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine solver_init_c(csolver,nn,tolerance) bind(c)                     !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type(iterative_solver_c), intent(inout) :: csolver
    integer(c_int), intent(in), value :: nn
    real(c_double), intent(in), value :: tolerance
    ! local variables
    class(iterative_solver), pointer :: fsolver

    csolver%nn = nn
    call iterative_solver_f(fsolver,csolver)
    call fsolver%init(nn,tolerance)

end subroutine solver_init_c



!--------------------------------------------------------------------------!
subroutine solve_c(csolver,cmat,x,b,cpc,n) bind(c)                         !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type(iterative_solver_c), intent(in) :: csolver
    type(sparse_matrix_c), intent(in) :: cmat
    real(c_double), intent(in) :: b(n)
    real(c_double), intent(inout) :: x(n)
    type(preconditioner_c), intent(in) :: cpc
    integer(c_int), intent(in), value :: n
    ! local variables
    class(iterative_solver), pointer :: fsolver
    class(sparse_matrix), pointer :: fmat
    class(preconditioner), pointer :: fpc
    integer :: mask(0)

    call iterative_solver_f(fsolver,csolver)
    call sparse_matrix_f(fmat,cmat)
    call preconditioner_f(fpc,cpc)

    call fsolver%solve(fmat,x,b,fpc,mask)

end subroutine solve_c






!==========================================================================!
!==========================================================================!
!==== C wrapper subroutines to preconditioner type-bound procedures    ====!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine preconditioner_init_c(cpc,cmat,level) bind(c)                   !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type(preconditioner_c), intent(inout) :: cpc
    type(sparse_matrix_c), intent(in) :: cmat
    integer(c_int), intent(in), value :: level
    ! local variables
    class(preconditioner), pointer :: fpc
    class(sparse_matrix), pointer :: fmat

    call preconditioner_f(fpc,cpc)
    call sparse_matrix_f(fmat,cmat)

    call fpc%init(fmat,level)

end subroutine preconditioner_init_c



!--------------------------------------------------------------------------!
subroutine precondition_c(cpc,cmat,x,b,n) bind(c)                          !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type(preconditioner_c), intent(inout) :: cpc
    type(sparse_matrix_c), intent(in) :: cmat
    real(c_double), intent(inout) :: x(n)
    real(c_double), intent(in) :: b(n)
    integer(c_int), intent(in), value :: n
    ! local variables
    class(preconditioner), pointer :: fpc
    class(sparse_matrix), pointer :: fmat
    integer :: mask(0)

    call preconditioner_f(fpc,cpc)
    call sparse_matrix_f(fmat,cmat)

    call fpc%precondition(fmat,x,b,mask)

end subroutine precondition_c




end module solver_c
