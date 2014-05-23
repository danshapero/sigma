module jacobi_solvers

use types, only: dp
use linear_operator_interface

implicit none


!--------------------------------------------------------------------------!
type, extends(linear_solver) :: jacobi_solver                              !
!--------------------------------------------------------------------------!
    ! inverse of diagonal entries of matrix
    real(dp), allocatable :: idiag(:)
contains
    procedure :: setup => jacobi_setup
    procedure :: linear_solve => jacobi_solve
    procedure :: destroy => jacobi_destroy
end type jacobi_solver


contains



!--------------------------------------------------------------------------!
function jacobi()                                                          !
!--------------------------------------------------------------------------!
    class(linear_solver), pointer :: jacobi

    allocate(jacobi_solver::jacobi)

end function jacobi



!--------------------------------------------------------------------------!
subroutine jacobi_setup(solver,A)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(jacobi_solver), intent(inout) :: solver
    class(linear_operator), intent(in)  :: A
    ! local variables
    integer :: i

    solver%nn = A%nrow

    if (A%ncol/=A%nrow) then
        print *, 'Cannot make a Jacobi solver for a non-square matrix'
        print *, 'Terminating.'
        call exit(1)
    endif

    if (.not.solver%initialized) then
        allocate(solver%idiag(solver%nn))

        solver%initialized = .true.
    endif

    do i=1,A%nrow
        solver%idiag(i) = 1.0_dp/A%get_value(i,i)
    enddo

end subroutine jacobi_setup



!--------------------------------------------------------------------------!
subroutine jacobi_solve(solver,A,x,b)                                      !
!--------------------------------------------------------------------------!
    class(jacobi_solver), intent(inout) :: solver
    class(linear_operator), intent(in)  :: A
    real(dp), intent(inout)             :: x(:)
    real(dp), intent(in)                :: b(:)

    associate(idiag => solver%idiag)

    x = idiag*b

    end associate

end subroutine jacobi_solve



!--------------------------------------------------------------------------!
subroutine jacobi_destroy(solver)                                          !
!--------------------------------------------------------------------------!
    class(jacobi_solver), intent(inout) :: solver

    solver%nn = 0
    solver%initialized = .false.

    deallocate(solver%idiag)

end subroutine jacobi_destroy





end module jacobi_solvers
