module linear_operator_adjoints

use types, only: dp
use linear_operator_interface

implicit none



!--------------------------------------------------------------------------!
type, extends(linear_operator) :: operator_adjoint                         !
!--------------------------------------------------------------------------!
    class(linear_operator), pointer :: op => null()
contains
    procedure :: get_value => operator_adjoint_get_value
    procedure :: matvec_add => operator_adjoint_matvec_add
    procedure :: matvec_t_add => operator_adjoint_matvec_t_add
    procedure :: destroy => operator_adjoint_destroy
end type operator_adjoint



contains




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
            call B%op%add_reference()
    end select

end function adjoint



!--------------------------------------------------------------------------!
function operator_adjoint_get_value(A,i,j) result(z)                       !
!--------------------------------------------------------------------------!
    class(operator_adjoint), intent(in) :: A
    integer, intent(in) :: i,j
    real(dp) :: z

    z = A%op%get_value(j,i)

end function operator_adjoint_get_value



!--------------------------------------------------------------------------!
subroutine operator_adjoint_matvec_add(A,x,y)                              !
!--------------------------------------------------------------------------!
    class(operator_adjoint), intent(in) :: A
    real(dp), intent(in) :: x(:)
    real(dp), intent(inout) :: y(:)

    call A%op%matvec_t_add(x,y)

end subroutine operator_adjoint_matvec_add



!--------------------------------------------------------------------------!
subroutine operator_adjoint_matvec_t_add(A,x,y)                            !
!--------------------------------------------------------------------------!
    class(operator_adjoint), intent(in) :: A
    real(dp), intent(in) :: x(:)
    real(dp), intent(inout) :: y(:)

    call A%op%matvec_add(x,y)

end subroutine operator_adjoint_matvec_t_add



!--------------------------------------------------------------------------!
subroutine operator_adjoint_destroy(A)                                     !
!--------------------------------------------------------------------------!
    class(operator_adjoint), intent(inout) :: A

    call A%op%remove_reference()

    if (A%op%reference_count <= 0) then
        call A%op%destroy()
        deallocate(A%op)
    endif

    nullify(A%op)

    A%reference_count = 0

end subroutine operator_adjoint_destroy



end module linear_operator_adjoints

