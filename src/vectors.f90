module vectors

use types

implicit none



!--------------------------------------------------------------------------!
type :: vector                                                             !
!--------------------------------------------------------------------------!
    real(dp), allocatable :: val(:)
    integer :: nn, nc
    integer, allocatable :: ptr(:)
contains
    procedure :: init_vector
    procedure :: init_multi_vector
    procedure :: vec_get_value
    procedure :: vec_set_value
    procedure :: vec_add_value
    procedure :: vec_get_value_multi_index
    procedure :: vec_set_value_multi_index
    procedure :: vec_add_value_multi_index
    generic :: init => init_vector, init_multi_vector
    generic :: get_value => vec_get_value, vec_get_value_multi_index
    generic :: set_value => vec_set_value, vec_set_value_multi_index
    generic :: add_value => vec_add_value, vec_add_value_multi_index
    procedure :: zero => vec_zero
end type vector


contains



!--------------------------------------------------------------------------!
subroutine init_vector(v,nn)                                               !
!--------------------------------------------------------------------------!
    class(vector), intent(inout) :: v
    integer, intent(in) :: nn

    v%nn = nn
    v%nc = 1

    allocate(v%val(nn),v%ptr(2))

    v%val = 0.0_dp
    v%ptr = [1, nn+1]

end subroutine init_vector



!--------------------------------------------------------------------------!
subroutine init_multi_vector(v,nn_per_field)                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(vector), intent(inout) :: v
    integer, intent(in) :: nn_per_field(:)
    ! local variables
    integer :: i

    v%nn = sum(nn_per_field)
    v%nc = size(nn_per_field)

    allocate(v%val(v%nn),v%ptr(v%nc+1))
    
    v%val = 0.0_dp
    v%ptr(1) = 1
    do i=1,v%nc
        v%ptr(i+1) = v%ptr(i)+nn_per_field(i)
    enddo

end subroutine init_multi_vector



!--------------------------------------------------------------------------!
function vec_get_value(v,i)                                                !
!--------------------------------------------------------------------------!
    class(vector), intent(in) :: v
    integer, intent(in) :: i
    real(dp) :: vec_get_value

    vec_get_value = v%val(i)

end function vec_get_value



!--------------------------------------------------------------------------!
function vec_get_value_multi_index(v,multi_index)                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(vector), intent(in) :: v
    integer, intent(in) :: multi_index(2)
    real(dp) :: vec_get_value_multi_index
    ! local variables
    integer :: i,k

    i = multi_index(1)
    k = multi_index(2)

    vec_get_value_multi_index = v%val( v%ptr(k)+i-1 )

end function vec_get_value_multi_index



!--------------------------------------------------------------------------!
subroutine vec_set_value(v,i,val)                                          !
!--------------------------------------------------------------------------!
    class(vector), intent(inout) :: v
    integer, intent(in) :: i
    real(dp), intent(in) :: val

    v%val(i) = val

end subroutine vec_set_value



!--------------------------------------------------------------------------!
subroutine vec_set_value_multi_index(v,multi_index,val)                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(vector), intent(inout) :: v
    integer, intent(in) :: multi_index(2)
    real(dp), intent(in) :: val
    ! local variables
    integer :: i,k

    i = multi_index(1)
    k = multi_index(2)

    v%val( v%ptr(k)+i-1 ) = val

end subroutine vec_set_value_multi_index



!--------------------------------------------------------------------------!
subroutine vec_add_value(v,i,val)                                          !
!--------------------------------------------------------------------------!
    class(vector), intent(inout) :: v
    integer, intent(in) :: i
    real(dp), intent(in) :: val

    v%val(i) = v%val(i)+val

end subroutine vec_add_value



!--------------------------------------------------------------------------!
subroutine vec_add_value_multi_index(v,multi_index,val)                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(vector), intent(inout) :: v
    integer, intent(in) :: multi_index(2)
    real(dp), intent(in) :: val
    ! local variables
    integer :: i,k

    i = multi_index(1)
    k = multi_index(2)

    v%val( v%ptr(k)+i-1 ) = v%val( v%ptr(k)+i-1 )+val

end subroutine vec_add_value_multi_index



!--------------------------------------------------------------------------!
subroutine vec_zero(v)                                                     !
!--------------------------------------------------------------------------!
    class(vector), intent(inout) :: v

    v%val = 0.0_dp

end subroutine vec_zero


end module vectors
