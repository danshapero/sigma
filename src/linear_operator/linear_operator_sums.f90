module linear_operator_sums

use types, only: dp
use linear_operator_interface

implicit none



!--------------------------------------------------------------------------!
type, extends(linear_operator) :: operator_sum                             !
!--------------------------------------------------------------------------!
    integer :: num_summands
    type(linear_operator_pointer), allocatable :: summands(:)
contains
    procedure :: get_value => operator_sum_get_value
    procedure :: matvec_add => operator_sum_matvec_add
end type operator_sum



contains



!--------------------------------------------------------------------------!
function operator_sum_get_value(A,i,j) result(val)                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(operator_sum), intent(in) :: A
    integer, intent(in) :: i,j
    real(dp) :: val
    ! local variables
    integer :: k

    val = 0.0_dp

    do k=1,A%num_summands
        val = val + A%summands(k)%ap%get_value(i,j)
    enddo

end function operator_sum_get_value



!--------------------------------------------------------------------------!
subroutine operator_sum_matvec_add(A,x,y,trans)                            !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(operator_sum), intent(in) :: A
    real(dp), intent(in) :: x(:)
    real(dp), intent(inout) :: y(:)
    logical, intent(in), optional :: trans
    ! local variables
    integer :: k

    do k=1,A%num_summands
        call A%summands(k)%ap%matvec_add(x,y,trans)
    enddo

end subroutine operator_sum_matvec_add




end module linear_operator_sums
