module linear_operators

use linear_operator_interface
use linear_operator_sums
use linear_operator_products
use linear_operator_adjoints

implicit none


interface assignment(=)
    module procedure assign_operators
end interface

interface operator(+)
    module procedure add_operators
end interface

interface operator(*)
    module procedure multiply_operators
end interface



contains



!--------------------------------------------------------------------------!
subroutine assign_operators(A,B)                                           !
!--------------------------------------------------------------------------!
    class(linear_operator), pointer, intent(out) :: A
    class(linear_operator), target, intent(in) :: B

    A => B

end subroutine assign_operators



!--------------------------------------------------------------------------!
function add_operators(A,B) result(C)                                      !
!--------------------------------------------------------------------------!
    class(linear_operator), target, intent(in) :: A, B
    class(linear_operator), pointer :: C

    ! Do some error checking
    if (A%nrow/=B%nrow .or. A%ncol/=B%ncol) then
        print *, 'Dimensions of operators to be summed are not consistent'
        call exit(1)
    endif

    ! Make a pointer to an operator_sum
    allocate(operator_sum::C)

    ! Set the dimension of C
    C%nrow = A%nrow
    C%ncol = A%ncol

    ! Make the summands of C point to A and B
    select type(C)
        type is(operator_sum)
            C%num_summands = 2
            allocate(C%summands(2))

            ! Make the operator_sum point to its two summands
            C%summands(1)%ap => A
            C%summands(2)%ap => B
    end select

end function add_operators



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
function adjoint(A) result(B)                                              !
!--------------------------------------------------------------------------!
    class(linear_operator), target, intent(in) :: A
    class(linear_operator), pointer :: B

    allocate(operator_adjoint::B)

    B%nrow = A%ncol
    B%ncol = A%nrow

    select type(B)
        type is(operator_adjoint)
            B%op => A
    end select

end function adjoint



end module linear_operators
