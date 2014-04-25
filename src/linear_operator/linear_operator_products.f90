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



!--------------------------------------------------------------------------!
interface operator(*)                                                      !
!--------------------------------------------------------------------------!
    module procedure multiply_operators
end interface



contains




!--------------------------------------------------------------------------!
function multiply_operators(A,B) result(C)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(linear_operator), target, intent(in) :: A, B
    class(linear_operator), pointer :: C
    ! local variables
    integer :: d

    ! Do some error checking
    if (A%ncol/=B%nrow) then
        print *, 'Dimensions of operators to be multiplied are inconsistent'
        call exit(1)
    endif

    ! Make a pointer to an operator_product
    allocate(operator_product::C)

    ! Set the dimension of C
    C%nrow = A%nrow
    C%ncol = B%ncol

    ! Make the factors of C point to A and B
    select type(C)
        type is(operator_product)
            C%num_products = 2
            allocate(C%products(2))

            ! Make space for temporary vectors in the operator product
            C%temp_vec_size = maxval([A%nrow, A%ncol, B%nrow, B%ncol])
            allocate(C%z1(C%temp_vec_size),C%z2(C%temp_vec_size))

            ! Make the operator_product point to its two factors
            C%products(1)%ap => A
            C%products(2)%ap => B
    end select

end function multiply_operators



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
