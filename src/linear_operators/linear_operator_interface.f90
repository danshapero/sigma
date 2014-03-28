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
end type linear_operator



!--------------------------------------------------------------------------!
type :: linear_operator_pointer                                            !
!--------------------------------------------------------------------------!
    class(linear_operator), pointer :: ap
end type linear_operator_pointer



!--------------------------------------------------------------------------!
type, abstract :: linear_solver                                            !
!--------------------------------------------------------------------------!
    class(linear_solver), pointer :: next
contains
    procedure(linear_solve_ifc), deferred :: solve
end type linear_solver



!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    subroutine opvec_add_ifc(A,x,y,trans)
        import :: linear_operator
        class(linear_operator), intent(in) :: A
        real(dp), intent(in) :: x(:)
        real(dp), intent(inout) :: y(:)
        logical, intent(in), optional :: trans
    end subroutine opvec_add_ifc
end interface



!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    subroutine linear_solve_ifc(solver,A,x,b)
        import :: linear_solver
        class(linear_solver), intent(inout) :: solver
        class(linear_operator), intent(in) :: A
        real(dp), intent(inout) :: x(:)
        real(dp), intent(in) :: b(:)
    end subroutine linear_solve_ifc
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

end subroutine linear_operator_get_value



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




end module linear_operator_interface
