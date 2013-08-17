module jacobi_preconditioners

use types
use sparse_matrices
use iterative_solvers

implicit none


!--------------------------------------------------------------------------!
type, extends(preconditioner) :: jacobi_preconditioner                     !
!--------------------------------------------------------------------------!
    real(dp), allocatable :: diag(:)
contains
    procedure :: init => jacobi_init
    procedure :: precondition => jacobi_precondition
    procedure :: free => jacobi_free
end type jacobi_preconditioner


contains




!--------------------------------------------------------------------------!
subroutine jacobi_init(pc,A,level)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(jacobi_preconditioner), intent(inout) :: pc
    class(sparse_matrix), intent(in)            :: A
    integer, intent(in)                         :: level
    ! local variables
    integer :: i

    pc%level = level
    pc%nn = A%nrow

    allocate(pc%diag(pc%nn))

    do i=1,A%nrow
        pc%diag(i) = A%get_value(i,i)
    enddo

end subroutine jacobi_init



!--------------------------------------------------------------------------!
subroutine jacobi_precondition(pc,A,x,b)                                   !
!--------------------------------------------------------------------------!
    class(jacobi_preconditioner), intent(inout) :: pc
    class(sparse_matrix), intent(in)            :: A
    real(dp), intent(inout)                     :: x(:)
    real(dp), intent(in)                        :: b(:)

    x = b/pc%diag

end subroutine jacobi_precondition



!--------------------------------------------------------------------------!
subroutine jacobi_free(pc)                                                 !
!--------------------------------------------------------------------------!
    class(jacobi_preconditioner), intent(inout) :: pc

    pc%level = 0
    pc%nn = 0

    deallocate(pc%diag)

end subroutine jacobi_free





end module jacobi_preconditioners



