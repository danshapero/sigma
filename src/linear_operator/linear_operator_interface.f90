!==========================================================================!
!==========================================================================!
module linear_operator_interface                                           !
!==========================================================================!
!==========================================================================!
!==== This module contains the definition of linear operator objects   ====!
!==== and linear solver objects.                                       ====!
!==========================================================================!
!==========================================================================!


use types, only: dp

implicit none




!--------------------------------------------------------------------------!
type, abstract :: linear_operator                                          !
!--------------------------------------------------------------------------!
    integer :: nrow, ncol
    class(linear_solver), pointer :: solver
contains
    procedure :: get_value => linear_operator_get_value
    procedure :: matvec => linear_operator_matvec
    procedure(opvec_add_ifc), deferred :: matvec_add
    procedure :: solve => linear_operator_solve
end type linear_operator



!--------------------------------------------------------------------------!
type :: linear_operator_pointer                                            !
!--------------------------------------------------------------------------!
    class(linear_operator), pointer :: ap
end type linear_operator_pointer



!--------------------------------------------------------------------------!
type, abstract :: linear_solver                                            !
!--------------------------------------------------------------------------!
    integer :: nn
    real(dp) :: tolerance
    class(linear_solver), pointer :: next
contains
    procedure(init_linear_solver_ifc), deferred :: init
    procedure(linear_solve_ifc), deferred :: linear_solve
    procedure :: linear_solve_pc
    generic :: solve => linear_solve, linear_solve_pc
end type linear_solver



!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    subroutine opvec_add_ifc(A,x,y,trans)
        import :: linear_operator, dp
        class(linear_operator), intent(in) :: A
        real(dp), intent(in) :: x(:)
        real(dp), intent(inout) :: y(:)
        logical, intent(in), optional :: trans
    end subroutine opvec_add_ifc
end interface



!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    subroutine init_linear_solver_ifc(solver,A)
        import :: linear_solver, linear_operator
        class(linear_solver), intent(inout) :: solver
        class(linear_operator), intent(in) :: A
    end subroutine init_linear_solver_ifc

    subroutine linear_solve_ifc(solver,A,x,b)
        import :: linear_solver, linear_operator, dp
        class(linear_solver), intent(inout) :: solver
        class(linear_operator), intent(in) :: A
        real(dp), intent(inout) :: x(:)
        real(dp), intent(in) :: b(:)
    end subroutine linear_solve_ifc

    subroutine linear_solve_pc_ifc(solver,A,x,b,pc)
        import :: linear_solver, linear_operator, dp
        class(linear_solver), intent(inout) :: solver
        class(linear_operator), intent(in) :: A
        real(dp), intent(inout) :: x(:)
        real(dp), intent(in) :: b(:)
        class(linear_solver), intent(inout) :: pc
    end subroutine linear_solve_pc_ifc
end interface



contains



!--------------------------------------------------------------------------!
function linear_operator_get_value(A,i,j) result(val)                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(linear_operator), intent(in) :: A
    integer, intent(in) :: i,j
    real(dp) :: val
    ! local variables
    real(dp) :: x(A%ncol), y(A%nrow)

    x(j) = 1.0_dp
    call A%matvec(x,y)
    val = y(i)

end function linear_operator_get_value



!--------------------------------------------------------------------------!
subroutine linear_operator_matvec(A,x,y,trans)                             !
!--------------------------------------------------------------------------!
    class(linear_operator), intent(in) :: A
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    logical, intent(in), optional :: trans

    y = 0.0_dp
    call A%matvec_add(x,y,trans)

end subroutine linear_operator_matvec



!--------------------------------------------------------------------------!
subroutine linear_operator_solve(A,x,b)                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(linear_operator), intent(in) :: A
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: b(:)
    ! local variables
    class(linear_solver), pointer :: solver

    solver => A%solver

    ! This subroutine is a facade for more complex operations that occur
    ! in a dedicated solver object contained in the operator itself
    call solver%solve(A,x,b)

end subroutine linear_operator_solve



!--------------------------------------------------------------------------!
subroutine linear_solve_pc(solver,A,x,b,pc)                                !
!--------------------------------------------------------------------------!
    class(linear_solver), intent(inout) :: solver
    class(linear_operator), intent(in)  :: A
    real(dp), intent(inout)             :: x(:)
    real(dp), intent(in)                :: b(:)
    class(linear_solver), intent(inout) :: pc

    call solver%solve(A,x,b)

end subroutine linear_solve_pc




end module linear_operator_interface
