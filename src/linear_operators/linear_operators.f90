module linear_operators

use linear_operator_interface
use linear_operator_sums

implicit none




!--------------------------------------------------------------------------!
function add_operators(A,B) result(C)                                      !
!--------------------------------------------------------------------------!
    class(linear_operator), pointer, intent(in) :: A, B
    class(linear_operator), pointer, intent(out) :: C

    allocate(operator_sum::C)
    allocate(C%summands(2))

    C%summands(1)%ap => A
    C%summands(2)%ap => B

end function add_operators




end module linear_operators
