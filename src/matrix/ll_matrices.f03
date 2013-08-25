module ll_matrices

use sparse_matrices
use ll_graphs
use types

implicit none


!--------------------------------------------------------------------------!
type, extends(sparse_matrix) :: ll_matrix                                  !
!--------------------------------------------------------------------------!
    real(dp), allocatable :: val(:)
    type(linked_list), allocatable :: ptrs(:)
    integer :: last
    class(ll_graph), pointer :: g
    ! procedure pointers to implementations of matrix operations
    ! Should these be private?
    procedure(ll_find_entry_ifc), pointer :: find_entry
    procedure(ll_matvec_ifc), pointer :: matvec_impl, matvec_t_impl
    procedure(ll_permute_ifc), pointer :: left_permute_impl, &
                                        & right_permute_impl
contains
    ! front-ends to matrix operations
    procedure :: init => ll_matrix_init
    procedure :: assemble => ll_assemble
    procedure :: neighbors => ll_matrix_neighbors
    procedure :: get_value => ll_get_value
    procedure :: set_value => ll_set_value, add_value => ll_add_value
    procedure :: sub_matrix_add => ll_sub_matrix_add
    procedure :: left_permute => ll_left_permute, &
                & right_permute => ll_right_permute
    procedure :: matvec => ll_matvec, matvec_t => ll_matvec_t
    procedure, private :: ll_set_value_not_preallocated
end type ll_matrix



!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    function ll_find_entry_ifc(A,i,j)
        import :: ll_matrix
        class(ll_matrix), intent(in) :: A
        integer, intent(in) :: i,j
        integer :: ll_find_entry_ifc
    end function ll_find_entry_ifc

    subroutine ll_matvec_ifc(A,x,y)
        import :: ll_matrix, dp
        class(ll_matrix), intent(in) :: A
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
    end subroutine

    subroutine ll_permute_ifc(A,p)
        import :: ll_matrix
        class(ll_matrix), intent(inout) :: A
        integer, intent(in) :: p(:)
    end subroutine ll_permute_ifc

end interface




contains



!==========================================================================!
!==========================================================================!
!== Front-ends to matrix operations                                      ==!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine ll_matrix_init(A,nrow,ncol)                                     !
!--------------------------------------------------------------------------!
    class(ll_matrix), intent(inout) :: A
    integer, intent(in) :: nrow, ncol

    A%nrow = nrow
    A%ncol = ncol
    A%max_degree = 0

end subroutine ll_matrix_init



!--------------------------------------------------------------------------!
subroutine ll_assemble(A,g)                                                !
!--------------------------------------------------------------------------!
    class(ll_matrix), intent(inout) :: A
    class(ll_graph), pointer, intent(in) :: g
    ! local variables
    integer :: i,k

    A%g => g

    if (A%orientation=='col') then
        A%ncol = g%n
        A%nrow = g%m

        A%find_entry => llc_find_entry

        A%matvec_impl   => llc_matvec
        A%matvec_t_impl => llr_matvec

        A%left_permute_impl => ll_matrix_permute_vals
        A%right_permute_impl => ll_matrix_permute_ptrs
    else
        A%nrow = g%n
        A%ncol = g%m

        A%find_entry => llr_find_entry

        A%matvec_impl   => llr_matvec
        A%matvec_t_impl => llc_matvec

        A%left_permute_impl => ll_matrix_permute_ptrs
        A%right_permute_impl => ll_matrix_permute_vals
    endif

    A%last = 1
    allocate(A%ptrs(g%n))
    do i=1,g%n
        do k=1,g%lists(i)%length
            call A%ptrs(i)%prepend(A%last)
            A%last = A%last+1
        enddo
    enddo

    A%nnz = g%ne
    A%max_degree = g%max_degree

    allocate(A%val(A%nnz))
    A%val = 0.0_dp

end subroutine ll_assemble



!--------------------------------------------------------------------------!
subroutine ll_matrix_neighbors(A,i,nbrs)                                   !
!--------------------------------------------------------------------------!
    class(ll_matrix), intent(in) :: A
    integer, intent(in) :: i
    integer, intent(out) :: nbrs(:)

    call A%g%neighbors(i,nbrs)

end subroutine ll_matrix_neighbors



!--------------------------------------------------------------------------!
function ll_get_value(A,i,j)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    real(dp) :: ll_get_value
    ! local variables
    integer :: k

    ll_get_value = 0.0_dp
    k = A%find_entry(i,j)
    if (k/=-1) ll_get_value = A%val(k)

end function ll_get_value



!--------------------------------------------------------------------------!
subroutine ll_set_value(A,i,j,val)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k

    k = A%find_entry(i,j)
    if (k/=-1) then
        A%val(k) = val
    else
        call A%ll_set_value_not_preallocated(i,j,val)
    endif

end subroutine ll_set_value



!--------------------------------------------------------------------------!
subroutine ll_add_value(A,i,j,val)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k

    k = A%find_entry(i,j)
    if (k/=-1) then
        A%val(k) = A%val(k)+val
    else
        call A%ll_set_value_not_preallocated(i,j,val)
    endif

end subroutine ll_add_value



!--------------------------------------------------------------------------!
subroutine ll_sub_matrix_add(A,B)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_matrix), intent(inout) :: A
    class(ll_matrix), intent(in)    :: B
    ! local variables
    integer :: i,j,k
    real(dp) :: Bij

    do i=1,B%g%n
        do k=1,B%g%lists(i)%length
            j = B%g%lists(i)%get_value(k)
            Bij = B%val( B%ptrs(i)%get_value(k) )
            call A%add_value(i,j,Bij)
        enddo
    enddo

end subroutine ll_sub_matrix_add



!--------------------------------------------------------------------------!
subroutine ll_left_permute(A,p)                                            !
!--------------------------------------------------------------------------!
    class(ll_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%left_permute_impl(p)

end subroutine ll_left_permute



!--------------------------------------------------------------------------!
subroutine ll_right_permute(A,p)                                           !
!--------------------------------------------------------------------------!
    class(ll_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%right_permute_impl(p)

end subroutine ll_right_permute



!--------------------------------------------------------------------------!
subroutine ll_matvec(A,x,y)                                                !
!--------------------------------------------------------------------------!
    class(ll_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)

    call A%matvec_impl(x,y)

end subroutine ll_matvec



!--------------------------------------------------------------------------!
subroutine ll_matvec_t(A,x,y)                                              !
!--------------------------------------------------------------------------!
    class(ll_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)

    call A%matvec_t_impl(x,y)

end subroutine ll_matvec_t



!--------------------------------------------------------------------------!
subroutine ll_set_value_not_preallocated(A,i,j,val)                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k,l
    real(dp), allocatable :: val_tmp(:)

    k = A%last
    if (k>size(A%val)) then
        allocate(val_tmp(max(2*A%nnz,10)))
        val_tmp(1:A%nnz) = A%val
        call move_alloc(from=val_tmp, to=A%val)
    endif

    if (A%orientation=='col') then
        call A%g%add_edge(j,i)
        call A%ptrs(j)%prepend(k)
    else
        call A%g%add_edge(i,j)
        call A%ptrs(i)%prepend(k)
    endif

    A%val(k) = val

    A%max_degree = A%g%max_degree
    A%nnz = A%nnz+1
    A%last = A%last+1

end subroutine ll_set_value_not_preallocated





!==========================================================================!
!==========================================================================!
!== Implementations of matrix operations                                 ==!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
function llr_find_entry(A,i,j)                                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    integer :: llr_find_entry
    ! local variables
    integer :: k

    llr_find_entry = -1
    k = A%g%find_edge(i,j)
    if (k/=-1) llr_find_entry = A%ptrs(i)%get_value(k)

end function llr_find_entry



!--------------------------------------------------------------------------!
function llc_find_entry(A,i,j)                                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    integer :: llc_find_entry
    ! local variables
    integer :: k

    llc_find_entry = -1
    k = A%g%find_edge(j,i)
    if (k/=-1) llc_find_entry = A%ptrs(j)%get_value(k)

end function llc_find_entry



!--------------------------------------------------------------------------!
subroutine llr_matvec(A,x,y)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    ! local variables
    integer :: i,j,k
    real(dp) :: z,Aij

    do i=1,A%g%n
        z = 0.0_dp
        do k=1,A%ptrs(i)%length
            j = A%g%lists(i)%get_value(k)
            Aij = A%val( A%ptrs(i)%get_value(k) )
            z = z+Aij*x(j)
        enddo
        y(i) = z
    enddo

end subroutine llr_matvec



!--------------------------------------------------------------------------!
subroutine llc_matvec(A,x,y)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    ! local variables
    integer :: i,j,k
    real(dp) :: z,Aij

    y = 0.0_dp
    do j=1,A%g%n
        z = x(j)
        do k=1,A%ptrs(j)%length
            i = A%g%lists(j)%get_value(k)
            Aij = A%val( A%ptrs(j)%get_value(k) )
            y(i) = y(i)+Aij*z
        enddo
    enddo

end subroutine llc_matvec



!--------------------------------------------------------------------------!
subroutine ll_matrix_permute_vals(A,p)                                     !
!--------------------------------------------------------------------------!
    class(ll_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%g%right_permute(p)

end subroutine ll_matrix_permute_vals



!--------------------------------------------------------------------------!
subroutine ll_matrix_permute_ptrs(A,p)                                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i,j,k
    type(linked_list) :: ptrs(A%g%n)

    do i=1,A%g%n
        do k=1,A%ptrs(i)%length
            j = A%ptrs(i)%get_value(k)
            call ptrs(i)%append(j)
        enddo
        call A%ptrs(i)%free()
    enddo

    do i=1,A%g%n
        do k=1,ptrs(i)%length
            j = ptrs(i)%get_value(k)
            call A%ptrs(p(i))%append(j)
        enddo
        call ptrs(i)%free()
    enddo

    call A%g%left_permute(p)

end subroutine ll_matrix_permute_ptrs





end module ll_matrices
