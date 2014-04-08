module jacobi_solvers

use types, only: dp
use linear_operator_interface

implicit none


!--------------------------------------------------------------------------!
type, extends(linear_solver) :: jacobi_solver                              !
!--------------------------------------------------------------------------!
    real(dp), allocatable :: diag(:)
contains
    procedure :: init => jacobi_init
    procedure :: linear_solve => jacobi_solve
    procedure :: free => jacobi_free
end type jacobi_solver


contains




!--------------------------------------------------------------------------!
subroutine jacobi_init(solver,A)                                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(jacobi_solver), intent(inout) :: solver
    class(linear_operator), intent(in)  :: A
    ! local variables
    integer :: i

    solver%nn = A%nrow

    allocate(solver%diag(solver%nn))

    do i=1,A%nrow
        solver%diag(i) = A%get_value(i,i)
    enddo

end subroutine jacobi_init



!--------------------------------------------------------------------------!
subroutine jacobi_solve(solver,A,x,b)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(jacobi_solver), intent(inout) :: solver
    class(linear_operator), intent(in)  :: A
    real(dp), intent(inout)             :: x(:)
    real(dp), intent(in)                :: b(:)
    ! local variables
    integer :: k

    associate(diag => solver%diag)

    x = b/diag

    end associate

end subroutine jacobi_solve



!--------------------------------------------------------------------------!
subroutine jacobi_free(solver)                                             !
!--------------------------------------------------------------------------!
    class(jacobi_solver), intent(inout) :: solver

    solver%nn = 0

    deallocate(solver%diag)

end subroutine jacobi_free





end module jacobi_solvers
