module linear_operator_products

use types, only: dp
use linear_operator_interface

implicit none



!--------------------------------------------------------------------------!
type, extends(linear_operator) :: operator_product                         !
!--------------------------------------------------------------------------!
    integer :: num_products, temp_vec_size
    type(linear_operator_pointer), allocatable :: products(:)
    real(dp), pointer :: z1(:), z2(:)
contains
    procedure :: matvec_add => operator_product_matvec_add
end type operator_product



contains



!--------------------------------------------------------------------------!
subroutine operator_product_matvec_add(A,x,y,trans)                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(operator_product), intent(in) :: A
    real(dp), intent(in) :: x(:)
    real(dp), intent(inout) :: y(:)
    logical, intent(in), optional :: trans
    ! local variables
    integer :: k
    real(dp), pointer :: z1(:), z2(:)

    z1 => A%z1
    z2 => A%z2

    z1(1:A%ncol) = x(1:A%ncol)
    z2(:) = 0.0_dp
    do k=A%num_products,1,-1
        call A%products(k)%ap%matvec(z1,z2,trans)
        z1(:) = z2(:)
    enddo
    y = y+z2

    z1(:) = 0.0_dp
    z2(:) = 0.0_dp

end subroutine operator_product_matvec_add



end module linear_operator_products
